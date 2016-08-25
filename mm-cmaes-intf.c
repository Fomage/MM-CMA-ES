#include <math.h>
#include <stdbool.h>
#include <emscripten.h>
#include <cmaes_interface.h>

mm_cmaes_t evo;

void cmaesInitialize(void) {
  double inxstart[] = {0, 0};
  double inrgstddev[] = {0.5, 0.5};
  long int inseed = 0;
  int lambda = 0;
  mm_cmaes_init(&evo, 4, 2, inxstart, inrgstddev, inseed, lambda, NULL);
}

double fitness;

void cmaesSetFitness(double newFitness) {
  fitness = newFitness;
}

double cmaesTargetFunc(double const *a) {
  EM_ASM_ARGS({
        _cmaesSetFitness(cmaesFunction($0, $1));
      }, a[0], a[1]);
  return fitness;
}

bool cmaesStep() {
  //const char *termination;
  if (evo.stahp) {
	//termination = cmaes_TestForTermination(&evo);
    /*EM_ASM_ARGS({
        cmaesSetTerminationMessage(Pointer_stringify($0))
    }, termination);*/
    return false;
  }

  const double *xmean;
  const double *rgD;
  int i,ivillage;

  mm_cmaes_run(&evo,&cmaesTargetFunc);
  
  EM_ASM(cmaesXMean = []);
  EM_ASM(cmaesEigenVectors = []);
  EM_ASM(cmaesAxisLengths = []);
  EM_ASM(cmaesSigma = []);
  EM_ASM(cmaesPoints = []);
  
  for(ivillage=0;ivillage<evo.nb_villages;ivillage++){
	  xmean = cmaes_GetPtr(evo.villages[ivillage], "xmean");
	  EM_ASM_ARGS({
		  cmaesXMean[$0] = ([$1, $2]);
		}, ivillage, xmean[0], xmean[1]);

	  EM_ASM_ARGS({
		  cmaesEigenVectors[$0] = ([[$1, $2], [$3, $4]]);
		}, ivillage,
		evo.villages[ivillage]->B[0][0],
		evo.villages[ivillage]->B[1][0],
		evo.villages[ivillage]->B[0][1], 
		evo.villages[ivillage]->B[1][1]);

	  rgD = cmaes_GetPtr(evo.villages[ivillage], "diag(D)");
	  EM_ASM_ARGS({
		  cmaesAxisLengths[$0] = ([$1, $2]);
		}, ivillage, rgD[0], rgD[1]);
	  
	  EM_ASM_ARGS({
		  cmaesSigma[$0] = ($1);
		}, ivillage, evo.villages[ivillage]->sigma);
		
	  /*EM_ASM_ARGS({
		  cmaesDivisionThreshold = ($0);
		}, evo.divisionThreshold);*/
		
	  EM_ASM_ARGS({
		  cmaesPoints[$0] = [];
		}, ivillage);
	  for (i = 0; i < cmaes_Get(evo.villages[ivillage], "popsize"); i++) {
		/* FIXME: workaround since EM_ASM_DOUBLE seems to truncate the result. */
		EM_ASM_ARGS({
			cmaesPoints[$0].push([$1, $2, $3, $4]);
		  }, ivillage,
			evo.villages[ivillage]->rgrgx[i][0],
			evo.villages[ivillage]->rgrgx[i][1],
			evo.villages[ivillage]->causeDivision[evo.villages[ivillage]->index[i]],
			evo.villages[ivillage]->clusters[evo.villages[ivillage]->index[i]]);
	  }
	  /*EM_ASM_ARGS({
		  cmaesShouldSplit[$0] = ($1);
		}, ivillage, evo.villages[ivillage]->shouldSplit);*/
  }
  
  EM_ASM(cmaesPushStep());
  return true;
}
