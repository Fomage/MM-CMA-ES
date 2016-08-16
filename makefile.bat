
emcc -o mm-cmaes-intf.js -Ic-cmaes-master/src/ mm-cmaes-intf.c c-cmaes-master/src/cmaes.c -s EXPORTED_FUNCTIONS="['_cmaesInitialize', '_cmaesStep', '_cmaesSetFitness']" 

