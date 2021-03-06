/* --------------------------------------------------------- */
/* --- File: cmaes.h ----------- Author: Nikolaus Hansen --- */
/* ---------------------- last modified: IX 2010         --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*
     CMA-ES for non-linear function minimization.

     Copyright (C) 1996, 2003-2010  Nikolaus Hansen.
     e-mail: nikolaus.hansen (you know what) inria.fr

     License: see file cmaes.c

*/
#ifndef NH_cmaes_h /* only include ones */
#define NH_cmaes_h

#include <time.h>

typedef struct
/* cmaes_random_t
 * sets up a pseudo random number generator instance
 */
{
  /* Variables for Uniform() */
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int *rgrand;

  /* Variables for Gauss() */
  short flgstored;
  double hold;

} cmaes_random_t;

typedef struct
/* cmaes_timings_t
 * time measurement, used to time eigendecomposition
 */
{
  /* for outside use */
  double totaltime; /* zeroed by calling re-calling cmaes_timings_start */
  double totaltotaltime;
  double tictoctime;
  double lasttictoctime;

  /* local fields */
  clock_t lastclock;
  time_t lasttime;
  clock_t ticclock;
  time_t tictime;
  short istic;
  short isstarted;

  double lastdiff;
  double tictoczwischensumme;
} cmaes_timings_t;

typedef struct
/* cmaes_readpara_t
 * collects all parameters, in particular those that are read from
 * a file before to start. This should split in future?
 */
{
  char * filename;  /* keep record of the file that was taken to read parameters */
  short flgsupplemented;

  /* input parameters */
  int N; /* problem dimension, must stay constant, should be unsigned or long? */
  unsigned int seed;
  double * xstart;
  double * typicalX;
  int typicalXcase;
  double * rgInitialStds;
  double * rgDiffMinChange;

  /* termination parameters */
  double stopMaxFunEvals;
  double facmaxeval;
  double stopMaxIter;
  struct { int flg; double val; } stStopFitness;
  double stopTolFun;
  double stopTolFunHist;
  double stopTolX;
  double stopTolUpXFactor;

  /* internal evolution strategy parameters */
  short flgNoRandom;   /* experimental */
  int lambda;          /* -> mu, <- N */
  int mu;              /* -> weights, (lambda) */
  double mucov, mueff; /* <- weights */
  double *weights;     /* <- mu, -> mueff, mucov, ccov */
  double damps;        /* <- cs, maxeval, lambda */
  double cs;           /* -> damps, <- N */
  double ccumcov;      /* <- N */
  double ccov;         /* <- mucov, <- N */
  double diagonalCov;  /* number of initial iterations */
  double divisionThreshold;
  struct { int flgalways; double modulo; double maxtime; } updateCmode;
  double facupdateCmode;

  /* supplementary variables */

  char *weigkey;
  char resumefile[99];
  const char **rgsformat;
  void **rgpadr;
  const char **rgskeyar;
  double ***rgp2adr;
  int n1para, n1outpara;
  int n2para;
} cmaes_readpara_t;

struct cmaes_type
/* cmaes_t
 * CMA-ES "object"
 */
{
  const char *version;
  /* char *signalsFilename; */
  cmaes_readpara_t sp;
  cmaes_random_t rand; /* random number generator */

  double sigma;  /* step size */

  double *rgxmean;  /* mean x vector, "parent" */
  double *rgxbestever;
  double **rgrgx;   /* range of x-vectors, lambda offspring */
  int *index;       /* sorting index of sample pop. */
  double *arFuncValueHist;
  int *causeDivision;
  int *clusters;

  short flgIniphase; /* not really in use anymore */
  short flgStop;

  double chiN;
  double **C;  /* lower triangular matrix: i>=j for C[i][j] */
  double **B;  /* matrix with normalize eigenvectors in columns */
  double *rgD; /* axis lengths */

  double *rgpc;
  double *rgps;
  double *rgxold;
  double *rgout;
  double *rgBDz;   /* for B*D*z */
  double *rgdTmp;  /* temporary (random) vector used in different places */
  double *rgFuncValue;
  double *publicFitness; /* returned by cmaes_init() */
  double **grid; /* used when flgNoRandom is set */

  double gen; /* Generation number */
  double countevals;
  double state; /* 1 == sampled, 2 == not in use anymore, 3 == updated */

  double maxdiagC; /* repeatedly used for output */
  double mindiagC;
  double maxEW;
  double minEW;

  char sOutString[330]; /* 4x80 */

  short flgEigensysIsUptodate;
  short flgCheckEigen; /* control via cmaes_signals.par */
  double genOfEigensysUpdate;
  cmaes_timings_t eigenTimings;

  double dMaxSignifKond;
  double dLastMinEWgroesserNull;

  short flgresumedone;

  time_t printtime;
  time_t writetime; /* ideally should keep track for each output file */
  time_t firstwritetime;
  time_t firstprinttime;

  short shouldSplit;
  short canSplit;

  struct cmaes_type* other;/* used in case of split */
  int splitGen;
};
typedef struct cmaes_type cmaes_t;


typedef struct
/* the MM-CMA-ES object, managing the different cmaes instances */
{
  /* parameters */
  int max_villages;
  int recoveryTimeAfterSplit;
  int tooYoungToMerge;
  double fusionThreshold;
  double fusionFactor;

  /* measures */
  int countevals;
  double fbestever;
  double* xbestever;
  int nbSplits;
  int nbMerges;

  /*internal variables */
  int nb_villages;

  cmaes_t** villages;
  double*const** pop;

  char stahp;
  char allowSplit;
} mm_cmaes_t;

#endif
