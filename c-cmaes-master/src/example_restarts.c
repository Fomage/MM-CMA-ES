#include <stdio.h>
#include <string.h> /* strncmp */
#include <math.h>
#include <stdlib.h>
#include <stddef.h> /* NULL */
#include <pthread.h>
#include "cmaes_interface.h"
/* result type for performances analysis */
typedef struct
{
    int countevals;
    double fitness;
    int nbSplits;
    int nbMerges;
} result;

/* parameters type for performance analysis */
typedef struct
{
    int N;
    int max_villages;
    double divisionThreshold;
    double fusionThreshold;
} parameters;

/* argument for the thread functions */
#define MAX_THREADS 8
typedef struct
{
	double (*pfun)(double const *);
	char *filename;
	char *used;
	result *m_result;
	parameters *m_parameters;
	pthread_mutex_t *used_mutex,*console_mutex;
} threadArgument;

/*___________________________________________________________________________
 *
 * Function Declarations
 *___________________________________________________________________________
 */
double **OrthogonalBasis(int DIM);
double f_rosenbrock( double const *x);
double f_rand( double const *x);
double f_constant( double const *x);
double f_kugelmin1( double const *x);
double f_sphere( double const *x);
double f_stepsphere( double const *x);
double f_cigar( double const *x);
double f_cigtab( double const *x);
double f_tablet( double const *x);
double f_elli( double const *x);
double f_ellirot( double const *x);
double f_elli100( double const *x);
double f_ellinumtest( double const *x);
double f_parabR( double const *x);
double f_sharpR( double const *x);
double f_diffpow( double const *x);
double f_diffpowrot( double const *x);
double f_gleichsys5( double const *x);
double f_rastrigin( double const *x);

result optimize(double(*pFun)(double const *),
        parameters p,
		char *input_parameter_filename);
void writeResults();

void *threadMain(void* arg);

extern void   cmaes_random_init( cmaes_random_t *, long unsigned seed /*=0=clock*/);
extern void   cmaes_random_exit( cmaes_random_t *);
extern double cmaes_random_Gauss( cmaes_random_t *); /* (0,1)-normally distributed */


#define N_MAX 40
#define N_STEP 1
#define SUCCESS_THRESHOLD 1e-13

/*___________________________________________________________________________
//___________________________________________________________________________
//
// reads from file "cmaes_initials.par" here and in cmaes_init()
//___________________________________________________________________________
 */
int main(int argn, char **args)
{
	typedef double (*pfun_t)(double const *);
	pfun_t rgpFun[99];  /* array (range) of pointer to objective function */
	char *filename = "cmaes_initials.par"; /* input parameter file */
	const char *s = "testOutput.dat";
    FILE *fp;

	int i,N;

	/* Put together objective functions */
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
	rgpFun[18] = f_gleichsys5;
	rgpFun[19] = f_rand;
	rgpFun[20] = f_constant;
	rgpFun[21] = f_stepsphere;
	rgpFun[22] = f_ellirot;
	rgpFun[23] = f_diffpowrot;
	rgpFun[24] = f_rastrigin;
	//int maxnb = 24;

    int nb=0;

	/* specify jobs to do, i.e every optimize() call to do */
	#define NB_JOBS 8
	parameters *jobsParameters = (parameters*) calloc(NB_JOBS,sizeof(parameters));
	result *jobsResults = (result*) calloc(NB_JOBS,sizeof(result));
	char* used = (char*) calloc(NB_JOBS,sizeof(char));

	for(i=0;i<NB_JOBS;i++){
		jobsParameters[i].N=i+2;
		jobsParameters[i].max_villages=i+2;
		jobsParameters[i].fusionThreshold=0;
	}

	pthread_mutex_t usedMutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_init(&usedMutex,NULL);
	pthread_mutex_t consoleMutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_init(&consoleMutex,NULL);

    /* launch the threads */
	pthread_t thr[MAX_THREADS];
    for(i=1;i<MAX_THREADS;i++){
		threadArgument theArg;
		theArg.pfun = rgpFun[nb];
		theArg.filename = filename;
		theArg.used = used;
		theArg.m_result = jobsResults;
		theArg.m_parameters = jobsParameters;
		theArg.used_mutex = &usedMutex;
		theArg.console_mutex = &consoleMutex;
        if(pthread_create(&thr[i],NULL,threadMain,&theArg)){
			printf("Error creating thread %d\n",i);
			return -1;
		}
	}

	/* do the job on this thread too */
	threadArgument theArg;
	theArg.pfun = rgpFun[nb];
	theArg.filename = filename;
	theArg.used = used;
	theArg.m_result = jobsResults;
	theArg.m_parameters = jobsParameters;
	theArg.used_mutex = &usedMutex;
	theArg.console_mutex = &consoleMutex;
	threadMain(&theArg);

	/* wait for other threads */
	for(i=1;i<MAX_THREADS;i++){
		pthread_join(thr[i],NULL);
	}

	/* write results */
	fp = fopen( s, "a");
    if(fp == NULL) {
        printf("cmaes_WriteToFile(): could not open '%s' with flag 'a'\n", s);
        return -2;
    }
	for(i=0;i<NB_JOBS;i++){
        fprintf(fp, "Dimension %d : counteval %d, divisionThreshold %.10e\n",
			jobsParameters[i].N,jobsResults[i].countevals,jobsParameters[i].divisionThreshold);
    }
    fclose(fp);

	/* clean up */
	free(jobsParameters);
	free(jobsResults);
	free(used);
	pthread_mutex_destroy(&usedMutex);
	pthread_mutex_destroy(&consoleMutex);

	return 0;

} /* main() */

void *threadMain(void* arg)
{
	threadArgument* a = (threadArgument*) arg;
	double dt,bestdt;
	result r;
	int i,job;

	/* search a job to do */
	for(job=0;job<NB_JOBS;job++){
		pthread_mutex_lock(a->used_mutex);
		if(a->used[job])
			pthread_mutex_unlock(a->used_mutex);
		else{
			a->used[job]=1;
			pthread_mutex_unlock(a->used_mutex);
			pthread_mutex_lock(a->console_mutex);
			printf("Begining job %d\n",job);
			pthread_mutex_unlock(a->console_mutex);
			/* begin parameter search */
			a->m_result[job].countevals=0;
			for(dt=.4;dt<4.;dt*=1.4){
				a->m_parameters[job].divisionThreshold=dt;
				r = optimize(a->pfun,a->m_parameters[job],a->filename);
				if((r.fitness < SUCCESS_THRESHOLD) &&
					((r.countevals < a->m_result[job].countevals) || (a->m_result[job].countevals == 0))){
					bestdt=dt;
					a->m_result[job]=r;
				}
			}
			a->m_parameters[job].divisionThreshold = bestdt;
		}
	}
}


/*___________________________________________________________________________
//___________________________________________________________________________
//
// Somewhat extended interface for optimizing pFun with cmaes_t
// implementing a restart procedure with increasing population size
//___________________________________________________________________________
 */

result optimize(double(*pFun)(double const *), parameters p, char * filename)
{
	result res;
	mm_cmaes_t evo;
    mm_cmaes_init(&evo, p.max_villages,
                  0,
                  0,
                  p.fusionThreshold,
                  0,
                  p.N, NULL, NULL, 0, 0, p.divisionThreshold, filename);

    //printf("Dimension %d, ft %.3e\n",p.N,p.fusionThreshold);
    //printf("Supplemented : %d, dt : %f\n",evo.villages[0]->sp.flgsupplemented, evo.villages[0]->sp.divisionThreshold);

    while(!evo.stahp){
	  mm_cmaes_run(&evo,pFun,0);
    }

    /*printf("MM CMA ES terminated in %d function evaluations, for a final fitness of %.5e.\n",
        evo.countevals, evo.fbestever);*/

    res.countevals = evo.countevals;
    res.fitness = evo.fbestever;
    res.nbSplits = evo.nbSplits;
    res.nbMerges = evo.nbMerges;

    mm_cmaes_exit(&evo);

	return res;
}

#if 1
/*___________________________________________________________________________
//___________________________________________________________________________
 */
double f_rand( double const *x)
{
	double d = (double)rand() / RAND_MAX;
	while (d == 0.)
		d = (double)rand() / RAND_MAX;
	return d;
}
double f_constant( double const *x)
{
	return 1;
}
#endif

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
/* a hack, memory is never released */
double **OrthogonalBasis(int DIM) {
	static int b_dim;
	static double **b;
	double sp;
	int i,j,k;
	cmaes_random_t R;

	if(b_dim != 0) { /* Initialization was done */

		if (b_dim != DIM) {
			printf("function OrthogonalBasis cannot change dimensionality in file example2.c");
			exit(0);
		}

		return b;
	}

	/* Otherwise initialize basis b */
	cmaes_random_init(&R, 2); /* TODO: choose not always the same basis? */

	/* allocate b */
	b = (double **) calloc((unsigned) DIM, sizeof(double*));
	if (!b) {
		printf("calloc failed in function OrthogonalBasis in file example2.c");
		exit(0);
	}
	for (i = 0; i < DIM; ++i) {
		b[i] = (double *) calloc((unsigned) DIM, sizeof(double));
		if (!b[i]) {
			printf("calloc failed in function Orthogonalbasis in file example2.c");
			exit(0);
		}
	}
	b_dim = DIM;

	/* generate orthogonal basis */
	for (i = 0; i < DIM; ++i) {
		/* sample components gaussian */
		for (j = 0; j < DIM; ++j)
			b[i][j] = cmaes_random_Gauss(&R);
		/* substract projection of previous vectors */
		for (j = i-1; j >= 0; --j) {
			for (sp = 0., k = 0; k < DIM; ++k)
				sp += b[i][k]*b[j][k]; /* scalar product */
			for (k = 0; k < DIM; ++k)
				b[i][k] -= sp * b[j][k]; /* substract */
		}
		/* normalize */
		for (sp = 0., k = 0; k < DIM; ++k)
			sp += b[i][k]*b[i][k]; /* squared norm */
		for (k = 0; k < DIM; ++k)
			b[i][k] /= sqrt(sp);
	}
	cmaes_random_exit(&R);

	return b;

} /* OrthogonalBasis(int DIM) */

/* ----------------------------------------------------------------------- */
double f_ellirot( double const *x)
{
	int i, k;
	double sum = 0., y;
	int DIM = (int)(x[-1]);
	double **B = OrthogonalBasis(DIM);

	if (DIM == 1)
		return x[0] * x[0];
	for (i = 0; i < DIM; ++i) {
		for (y = 0., k = 0; k < DIM; ++k)
			y += B[i][k] * x[k];
		sum += exp(log(1e3) * 2. * (double)(i)/(DIM-1)) * y*y;
	}
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
double f_diffpowrot( double const *x)
{
	int i, k;
	double sum = 0., y;
	int DIM = (int)(x[-1]);
	double **B = OrthogonalBasis(DIM);

	if (DIM == 1)
		return x[0] * x[0];
	for (i = 0; i < DIM; ++i) {
		for (y = 0., k = 0; k < DIM; ++k)
			y += B[i][k] * x[k];
		sum += pow(fabs(y), 2.+10*(double)(i)/(DIM-1));
	}
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

/*
  05/10/05: revised buggy comment on handling constraints by resampling
 */
