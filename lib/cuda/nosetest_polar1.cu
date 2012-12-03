/*
This is a very simple comparision of nosecuda_polar1.cu and nosecpu_polar1.cu. 
USAGE: 
1) compile it with nvcc:

nvcc -lm nosecuda_polar1.cu nosecpu_polar1.cu nosetest_polar1.cu -o nosetest_polar1
2) run it with an arbitrary trigger parameter. Trigger means that if the real or the imaginary part of S returned from cpu and gpu method differ for more than trigger a new line is printed.

./nosetest_polar1 [trigger]
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern float *nosecuda_polar_1(int* N, float *rU, float *iU, float *d, float *rho0);
extern float *nosecpu_polar_1(int* N, float *rU, float *iU, float *d, float *rho0);

float random01(){//returns random number between 0 and 1
	return rand()/(float)RAND_MAX;
}

int main(int argc, char **argv){
	//set dimensios of parameters
	int N[3]={4, 4, 300}; 
	float d[N[1]*N[0]]; 
	float rho0[N[0]]; 
	float rU[N[0]*N[1]*N[0]*N[1]*N[2]]; 
	float iU[N[0]*N[1]*N[0]*N[1]*N[2]]; 
	int i, j, a, b, t;

	//set parameters to random values
	for (i=0;i<N[0];i++){//set rho0, d
		rho0[i]= random01()*10.0;
		for (a=0;a<N[1];a++){
			d[i+N[0]*a]= random01()*10;//float between 0 and 100
		}
	}
	for (t=0;t<N[2];t++){//set rU, iU
		for (i=0;i<N[0];i++){//set d
			for (a=0;a<N[1];a++){
				for (j=0;j<N[0];j++){
					for (b=0;b<N[1];b++){
						rU[i+N[0]*(a+N[1]*(j+N[0]*(b+N[1]*t)))]= random01();
						iU[i+N[0]*(a+N[1]*(j+N[0]*(b+N[1]*t)))]= random01();
					}
				}
			}
		}
	}
	//now everything is set already and can be used
	printf("starting comparision\n");

	float trigger= 0.0;
	if (argc>1) trigger= atof(argv[1]);
	printf("trigger= %f\n", trigger);

	//execute both cpu and gpu methods
	float *SC, *SG;
	printf("executing cpu method\n");
	SC= (float*)nosecpu_polar_1(N, rU, iU, d, rho0);
	printf("executing gpu method\n");
	SG= (float*)nosecuda_polar_1(N, rU, iU, d, rho0);
	printf("comparing methods\n");
	
	//compare all values and print if differ for more than trigger
	for (t=0;t<N[2];t++){
		if ((fabs(SG[2*t]-SC[2*t])>=trigger)||(fabs(SG[2*t+1]-SC[2*t+1])>= trigger)) printf("for t=%i: SG[t]=%f+i%f, SC[t]=%f+i%f\n", t, SG[2*t], SG[2*t+1], SC[2*t], SC[2*t+1]);
	}
	
	//free memory
	free(SC);
	free(SG);

}
