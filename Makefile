all: mm-cmaes-intf.js

clean:
	rm -f mm-cmaes-intf.js

mm-cmaes-intf.js: mm-cmaes-intf.c
	emcc -o $@ -Ic-cmaes-master/src/ $< c-cmaes-master/src/cmaes.c -s EXPORTED_FUNCTIONS="['_cmaesInitialize', '_cmaesStep', '_cmaesSetFitness']"
