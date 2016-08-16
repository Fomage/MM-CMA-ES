#include <math.h>
#include <stdbool.h>
#include <emscripten.h>
#include <cmaes_interface.h>

cmaes_t evo;

void cmaesInitialize(void) {
  double inxstart[] = {0, 0};
  double inrgstddev[] = {0.5, 0.5};
  long int inseed = 0;
  int lambda = 0;
  cmaes_init(&evo, 2, inxstart, inrgstddev, inseed, lambda, NULL);
}

double fitness;

void cmaesSetFitness(double newFitness) {
  fitness = newFitness;
}

bool cmaesStep() {
  const char *termination = cmaes_TestForTermination(&evo);
  if (termination) {
    EM_ASM_ARGS({
        cmaesSetTerminationMessage(Pointer_stringify($0))
    }, termination);
    return false;
  }
  double *const*pop;
  const double *xmean;
  const double *rgD;
  int i;

  pop = cmaes_SamplePopulation(&evo);

  xmean = cmaes_GetPtr(&evo, "xmean");
  EM_ASM_ARGS({
      cmaesXMean = ([$0, $1]);
    }, xmean[0], xmean[1]);

  EM_ASM_ARGS({
      cmaesEigenVectors = ([[$0, $1], [$2, $3]]);
    }, evo.B[0][0], evo.B[1][0], evo.B[0][1], evo.B[1][1]);

  rgD = cmaes_GetPtr(&evo, "diag(D)");
  EM_ASM_ARGS({
      cmaesAxisLengths = ([$0, $1]);
    }, rgD[0], rgD[1]);
  
  EM_ASM_ARGS({
      cmaesSigma = ($0);
    }, evo.sigma);
	
  EM_ASM_ARGS({
      cmaesDivisionThreshold = ($0);
    }, evo.divisionThreshold);

  EM_ASM(cmaesPoints = []);
  for (i = 0; i < cmaes_Get(&evo, "popsize"); i++) {
    /* FIXME: workaround since EM_ASM_DOUBLE seems to truncate the result. */
    EM_ASM_ARGS({
        _cmaesSetFitness(cmaesFunction($0, $1));
      }, pop[i][0], pop[i][1]);
    evo.publicFitness[i] = fitness;
  }
  cmaes_UpdateDistribution(&evo, evo.publicFitness);
  for (i = 0; i < cmaes_Get(&evo, "popsize"); i++) {
    /* FIXME: workaround since EM_ASM_DOUBLE seems to truncate the result. */
    EM_ASM_ARGS({
        cmaesPoints.push([$0, $1, $2, $3]);
      }, pop[evo.index[i]][0], pop[evo.index[i]][1],evo.causeDivision[evo.index[i]], evo.clusters[evo.index[i]]);
  }
  EM_ASM_ARGS({
      cmaesShouldSplit = ($0);
    }, evo.shouldSplit);
  EM_ASM(cmaesPushStep());
  return true;
}
