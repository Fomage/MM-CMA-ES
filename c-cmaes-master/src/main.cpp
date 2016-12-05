#ifdef __cplusplus
    #include <cstdlib>
#else
    #include <stdlib.h>
#endif
#include "example_restart.h"
#include "cmaes_interface.h"
#include "cmaes.h"

#include "math.h"
#include "time.h"
#include "stdio.h"

int callSDL(int argc, char** argv);
int testParameter(int nFun, double theParameter);

int main ( int argc, char** argv )
{
    srand(time(NULL));

    printf("Starting the console.\n");

    int i,j;
    /*sscanf(argv[1],"%d",&i);
    int nFun=i%4;
    switch(nFun){
    case 0:
        break;
    case 1:
        break;
    case 2:
        nFun=5;
        break;
    case 3:
        nFun=14;
        break;
    }
    double theParameter = pow(10,(i+2)/4);
    printf("Testing function %d with parameter value %2.5e\n",nFun,theParameter);*/

    /*double theParameter;
    for(i=0;i<100;i++){
        printf("Run %d\n",i);
        for(j=2;j<10;j++){
            theParameter = pow(10,j);
            testParameter(14,theParameter);
        }
    }*/

    return callSDL(argc, argv);
}

#include <SDL/SDL.h>

#define XSIZE 600
#define YSIZE 600

typedef struct{
    double xoffset,yoffset,zoom,ffac;
} scope_t;

void screenToPoint(int x, int y, double *xp, double *yp, scope_t s){
    *xp=((double)x/(double)XSIZE -.5+s.xoffset)*s.zoom;
    *yp=((double)y/(double)YSIZE -.5+s.yoffset)*s.zoom;
}

void pointToScreen(double xp, double yp, int *x, int *y, scope_t s){
    *x=((xp)/s.zoom+.5-s.xoffset)*XSIZE;
    *y=((yp)/s.zoom+.5-s.yoffset)*YSIZE;
}

int callSDL(int argc,char**argv){

    int i=0;

    /* range the benchmark functions */
    typedef double (*pfun_t)(double const *);
	pfun_t rgpFun[15];
    rgpFun[0] = f_sphere;
	rgpFun[1] = f_elli;
	rgpFun[2] = f_cigar;
	rgpFun[3] = f_cigtab;
	rgpFun[4] = f_tablet;
	rgpFun[5] = f_rosenbrock;
	rgpFun[6] = f_parabR;
	rgpFun[7] = f_sharpR;
	rgpFun[8] = f_diffpow;
	rgpFun[9] = f_kugelmin1;
	rgpFun[10] = f_ellinumtest;
	rgpFun[11] = f_elli100;
	rgpFun[12] = f_gleichsys5;
	rgpFun[13] = f_stepsphere;
	rgpFun[14] = f_rastrigin;
	int nFun=14;

    /* parameters */
    int N=2;
    parameters p;
    p.divisionThreshold = 0.7*log(0.7*N)+0.4;
	p.fusionThreshold = .7;//0.000346*N*N*N+0.4;
	p.cdiv = .5;
	p.recovery = 1;
	p.tooYoung = 10;
	p.mud = 1+floor(.75*log(N));
	p.max_villages = 20;//N+2;
	p.N = N;
	double* xinit=(double*)malloc(N*sizeof(double));for(i=0;i<N;i++)xinit[i]=.5;
	double* sinit=(double*)malloc(N*sizeof(double));for(i=0;i<N;i++)sinit[i]=.3;

	/* mm-cma-es init */
	result res;
	mm_cmaes_t evo;
    mm_cmaes_init(&evo, p.max_villages,
                  p.recovery,
                  p.tooYoung,
                  p.fusionThreshold,
                  0,//fusion factor
                  p.N,
                  xinit, sinit, 0, 0, p.divisionThreshold, p.cdiv, p.mud, NULL);
    evo.fbestever = (*rgpFun[nFun])(evo.villages[0]->rgxmean);

    // initialize SDL video
    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 )
    {
        printf( "Unable to init SDL: %s\n", SDL_GetError() );
        return 1;
    }

    // make sure SDL cleans up before exit
    atexit(SDL_Quit);

    //enable key repeat
    SDL_EnableKeyRepeat(500,200);

    // create a new window
    SDL_Surface* screen = SDL_SetVideoMode(XSIZE, YSIZE, 16,
                                           SDL_HWSURFACE|SDL_DOUBLEBUF);
    if ( !screen )
    {
        printf("Unable to set %dx%d video: %s\n",XSIZE,YSIZE,SDL_GetError());
        return 1;
    }
    scope_t scope;
    scope.xoffset=scope.yoffset=0;
    scope.ffac=scope.zoom=1;

    //init colours
    Uint32 *colours = (Uint32*)malloc(evo.max_villages*sizeof(Uint32));
    for(int i=0;i<evo.max_villages;i++){
        double temp=((double)i/(double)evo.max_villages)*512;
        if(temp<256)
            colours[i]=SDL_MapRGB(screen->format,255-temp,temp,0);
        else
            colours[i]=SDL_MapRGB(screen->format,0,255-temp,temp);
    }
    char window_title[100]="Visual CMAES";
    SDL_WM_SetCaption(window_title,NULL);

    // program main loop
    bool done = false;
    bool redraw = true;
    bool fastforward = false;
    while (!done)
    {
        if(fastforward){
            if(!evo.stahp){
                mm_cmaes_run(&evo,rgpFun[nFun],0);
                sprintf(window_title,"Running : %f",evo.fbestever);
                SDL_WM_SetCaption(window_title,NULL);
                redraw=true;
                printf("%d\n",evo.countevals);
            }else{
                sprintf(window_title,"Finished : %f",evo.fbestever);
                SDL_WM_SetCaption(window_title,NULL);
                fastforward=false;
            }
        }

        // message processing loop
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            // check for messages
            switch (event.type)
            {
                // exit if the window is closed
            case SDL_QUIT:
                done = true;
                break;

                // check for keypresses
            case SDL_KEYDOWN:
                {
                    // exit if ESCAPE is pressed
                    if (event.key.keysym.sym == SDLK_ESCAPE)
                        done = true;
                    else if (event.key.keysym.sym == SDLK_SPACE){
                        fastforward = !fastforward;
                    } else if (event.key.keysym.sym == SDLK_RETURN){
                        if(!evo.stahp){
                            mm_cmaes_run(&evo,rgpFun[nFun],0);
                            sprintf(window_title,"Running : %f",evo.fbestever);
                            SDL_WM_SetCaption(window_title,NULL);
                            redraw=true;
                            printf("%d\n",evo.countevals);
                        }else{
                            sprintf(window_title,"Finished : %f",evo.fbestever);
                            SDL_WM_SetCaption(window_title,NULL);
                        }
                    } else if (event.key.keysym.sym == SDLK_LEFT){
                        scope.xoffset-=scope.zoom/8;
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_RIGHT){
                        scope.xoffset+=scope.zoom/8;
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_UP){
                        scope.yoffset-=scope.zoom/8;
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_DOWN){
                        scope.yoffset+=scope.zoom/8;
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_q){
                        scope.zoom*=2;
                        printf("Zoom out");
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_w){
                        scope.zoom/=2;
                        printf("Zoom in");
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_a){
                        scope.ffac*=2;
                        printf("Heighten");
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_s){
                        scope.ffac/=2;
                        printf("Flatten");
                        redraw=true;
                    } else if (event.key.keysym.sym == SDLK_z){
                        scope.ffac=scope.zoom=1;
                        scope.xoffset=scope.yoffset=0;
                        redraw=true;
                    }
                    break;

                }
            } // end switch
        } // end of message processing

        // DRAWING STARTS HERE

        if(redraw){
            redraw=false;
            //printf("Redrawing...\n");
            double *point=(double*)malloc(3*sizeof(double));
            point[0]=2;
            point++;

            // draw background
            SDL_Rect pixel;
            for(int x=0;x<XSIZE;x++){
                for(int y=0;y<YSIZE;y++){
                    screenToPoint(x,y,&point[0],&point[1],scope);
                    double eval=(*rgpFun[nFun])(point);
                    //eval!=0?printf("%f ",eval):NULL;
                    eval*=255.*scope.ffac;
                    pixel.h = pixel.w = 1;
                    pixel.x=x;pixel.y=y;
                    SDL_FillRect(screen,&pixel,SDL_MapRGB(screen->format,(int)eval,(int)eval,(int)eval));
                }
            }

            // draw points
            for(int v=0;v<evo.max_villages;v++){
                if(evo.villages[v]){
                    pixel.h = pixel.w = 2;
                    for(int i=0;i<evo.villages[v]->sp.lambda;i++){
                        int x,y;
                        pointToScreen(evo.villages[v]->rgrgx[i][0],evo.villages[v]->rgrgx[i][1],&x,&y,scope);
                        pixel.x=x;pixel.y=y;
                        SDL_FillRect(screen,&pixel,colours[v]);
                    }
                }
            }

            //draw ellipses
            for(int v=0;v<evo.max_villages;v++){
                if(evo.villages[v]){
                    pixel.h = pixel.w = 1;
                    for(double t=0;t<360;t++){
                        int x,y;
                        double myPoint[4];
                        myPoint[0]=evo.villages[v]->rgD[0] * evo.villages[v]->sigma * cos(t);
                        myPoint[1]=evo.villages[v]->rgD[1] * evo.villages[v]->sigma * sin(t);
                        myPoint[2]=myPoint[0]*evo.villages[v]->B[0][0] + myPoint[1]*evo.villages[v]->B[0][1];
                        myPoint[3]=myPoint[0]*evo.villages[v]->B[1][0] + myPoint[1]*evo.villages[v]->B[1][1];
                        myPoint[2]+=evo.villages[v]->rgxmean[0];
                        myPoint[3]+=evo.villages[v]->rgxmean[1];
                        pointToScreen(myPoint[2],myPoint[3],&x,&y,scope);
                        pixel.x=x;pixel.y=y;
                        SDL_FillRect(screen,&pixel,colours[v]);
                    }
                }
            }

            // DRAWING ENDS HERE

            // finally, update the screen :)
            SDL_Flip(screen);
        }

    } // end main loop

    res.countevals = evo.countevals;
    res.fitness = evo.fbestever;
    res.nbSplits = evo.nbSplits;
    res.nbMerges = evo.nbMerges;

    mm_cmaes_exit(&evo);

    free(xinit);
    free(sinit);
    free(colours);

    return 0;
}

int testParameter(int nFun, double theParameter)
{
    int i=0;

    /* range the benchmark functions */
    typedef double (*pfun_t)(double const *);
	pfun_t rgpFun[15];
    rgpFun[0] = f_sphere;
	rgpFun[1] = f_elli;
	rgpFun[2] = f_cigar;
	rgpFun[3] = f_cigtab;
	rgpFun[4] = f_tablet;
	rgpFun[5] = f_rosenbrock;
	rgpFun[6] = f_parabR;
	rgpFun[7] = f_sharpR;
	rgpFun[8] = f_diffpow;
	rgpFun[9] = f_kugelmin1;
	rgpFun[10] = f_ellinumtest;
	rgpFun[11] = f_elli100;
	rgpFun[12] = f_gleichsys5;
	rgpFun[13] = f_stepsphere;
	rgpFun[14] = f_rastrigin;

    /* parameters */
    int N=2;
    parameters p;
    p.divisionThreshold = 0.7*log(0.7*N)+0.4;
	p.fusionThreshold = .5;//0.000346*N*N*N+0.4;
	p.cdiv = .5;
	p.recovery = 1;
	p.tooYoung = 10;
	p.mud = 1+floor(.75*log(N));
	p.max_villages = 20;//N+2;
	p.N = N;
	double* xinit=(double*)malloc(N*sizeof(double));for(i=0;i<N;i++)xinit[i]=.5;
	double* sinit=(double*)malloc(N*sizeof(double));for(i=0;i<N;i++)sinit[i]=.3;

	/* mm-cma-es init */
	mm_cmaes_t evo;
    mm_cmaes_init(&evo, p.max_villages,
                  p.recovery,
                  p.tooYoung,
                  p.fusionThreshold,
                  0,//fusion factor
                  p.N,
                  xinit, sinit, 0, 0, p.divisionThreshold, p.cdiv, p.mud, NULL);
    evo.poison_threshold = theParameter;
    evo.maxevals = 1e5;
    evo.villages[0]->sp.sobol=0;
    evo.fbestever = (*rgpFun[nFun])(evo.villages[0]->rgxmean);

    while(!evo.stahp){
        mm_cmaes_run(&evo,rgpFun[nFun],0);
    }

    FILE *fp;
    fp = fopen("output.csv", "a");
    fprintf(fp, "%d,%f,%d,%f\n",nFun,theParameter,evo.countevals,evo.fbestever);
    fclose(fp);

    mm_cmaes_exit(&evo);

    free(xinit);
    free(sinit);

    return 0;
}

/*-------------------------------------------------------------------*/
/* benchmarking functions */

static double SQR(double d)
{
	return (d*d);
}

/* ----------------------------------------------------------------------- */
double f_stepsphere( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);
	for (i = 0; i < DIM; ++i)
		sum += floor(x[i]*x[i]);
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_sphere( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);
	if(DIM<=0)
        printf("Warning f_sphere : DIM=%d\n",DIM);
	for (i = 0; i < DIM; ++i)
		sum += x[i]*x[i];
    //printf("TEST f_sphere : return=%f\n",sum);
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_cigar( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	for (i = 1; i < DIM; ++i)
		sum += x[i]*x[i];
	sum *= 1e6;
	sum += x[0]*x[0];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_cigtab( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	sum = x[0]*x[0] + 1e8*x[DIM-1]*x[DIM-1];
	for (i = 1; i < DIM-1; ++i)
		sum += 1e4*x[i]*x[i];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_tablet( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	sum = 1e6*x[0]*x[0];
	for (i = 1; i < DIM; ++i)
		sum += x[i]*x[i];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_elli( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	if (DIM == 1)
		return x[0] * x[0];
	for (i = 0; i < DIM; ++i)
		sum += exp(log(1000.) * 2. * (double)(i)/(DIM-1)) * x[i]*x[i];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_elli100( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	if (DIM == 1)
		return x[0] * x[0];
	for (i = 0; i < DIM; ++i)
		sum += exp(log(100.) * 2. * (double)(i)/(DIM-1)) * x[i]*x[i];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_diffpow( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	if (DIM == 1)
		return x[0] * x[0];
	for (i = 0; i < DIM; ++i)
		sum += pow(fabs(x[i]), 2.+10*(double)(i)/(DIM-1));
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_kugelmin1( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);

	for (i = 1; i < DIM; ++i)
		sum += x[i]*x[i];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_rosenbrock( double const *x)
/*
	Rosenbrock's Function, generalized.
 */
{
	double qualitaet;
	int i;
	int DIM = (int)(x[-1]);
	qualitaet = 0.0;

	for( i = DIM-2; i >= 0; --i)
		qualitaet += 100.*SQR(SQR(x[i])-x[i+1]) + SQR(1.-x[i]);
	return ( qualitaet);
} /* f_rosenbrock() */

/* ----------------------------------------------------------------------- */
double f_parabR( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);
	for (i = 1; i < DIM; ++i)
		sum += x[i]*x[i];
	return -x[0] + 100.*sum;
}

/* ----------------------------------------------------------------------- */
double f_sharpR( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);
	for (i = 1; i < DIM; ++i)
		sum += x[i]*x[i];
	return -x[0] + 100*sqrt(sum);
}

/* ----------------------------------------------------------------------- */
double f_ellinumtest( double const *x)
{
	int i;
	double sum = 0.;
	int DIM = (int)(x[-1]);
	static double maxVerhaeltnis = 0.;
	if (maxVerhaeltnis == 0.)
	{
		for (maxVerhaeltnis = 1.;
				maxVerhaeltnis < 1e99 && maxVerhaeltnis < 2. * maxVerhaeltnis;
				maxVerhaeltnis *= 2.)
			if (maxVerhaeltnis == maxVerhaeltnis + 1.)
				break;
		maxVerhaeltnis *= 10.;
		maxVerhaeltnis = sqrt (maxVerhaeltnis);
	}
	if (DIM < 3)
		return x[0] * x[0];
	for (i = 1; i < DIM; ++i)
		sum += exp(log(maxVerhaeltnis) * 2. * (double)(i-1)/(DIM-2)) * x[i]*x[i];
	return sum;
}

/* ----------------------------------------------------------------------- */
double f_gleichsys5( double const *x)
/*
	Gleichungssystem 5-dimensional von Juergen Bremer
	Fuer jede Zeile soll gelten:
	 c_1*x[1] + c_2*x[2] + c_3*x[3] + c_4*x[4] + c_5*x[5] + c_0 = 0
	 Deswegen geht das Quadrat der linken Seite in die
	 Qualitaetsfunktion ein.
 */
{
	double qualitaet = 0.0;

#if 1
	static double koeff[5][6] =
	{/* c_1,   c_2,  c_3,   c_4,  c_5,   c_0 */
			{   4,   191,   27,   199,   21,   172},
			{ 191, 10883, 1413,  5402,  684, -8622},
			{  27,  1413,  191,  1032,  118,   -94},
			{ 199,  5402, 1032, 29203, 2331, 78172},
			{  21,   684,  118,  2331,  199,  5648}
	};
	int i, j;
	double sum;

	for( i = 0; i < 5; ++i)
	{
		sum = koeff[i][5];
		for ( j = 0; j < 5; ++j)
		{
			sum += koeff[i][j] * x[j];
		}
		qualitaet += sum * sum;
	}
#endif
return qualitaet;
} /* f_gleichsys5() */

/* ----------------------------------------------------------------------- */
double f_rastrigin( double const *x)
{
    int i;
	int DIM = (int)(x[-1]);
    double sum = 10.*DIM;
    for(i=0;i<DIM;i++)
        sum+=x[i]*x[i]-10*cos(2*M_PI*x[i]);
    return sum;
}


