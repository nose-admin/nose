/*
	This file contains method nosecuda_polar_1(...) which returns response function S(t) in time. Theory can be found in the white paper Absorption Spectrum at NOSE web.
*/

#include <stdio.h>
#include <stdlib.h>

/*
	Computes response functions S(t) in time using CPU.
	Parameters: 
		int N[3]: N[1] = number of levels in the ground state band
		          N[2] = number of levels in the excited state band
			  N[3] = number of the time steps
		float rU[N[0]*N[1]*N[0]*N[1]*N[2]] = real part of the coherence evolution superoperator in time. The superoperator has 4 indeces, the 5th dimension represents discrete time.

		float iU[N[0]*N[1]*N[0]*N[1]*N[2]] = imaginary part of the coherence evolution superoperator. The superoperator has 4 indeces, the 5th dimension represents discrete time.
		float d[N[0]*N[1]] = transition dipole moment elements
		float rho0[N[0]] = diagonal elements of the equilibrium ground state density matrix
		float S[2*N[0]] = response function in time (both real and imaginary parts)
*/

float *nosecpu_polar_1(int* N, float* rU, float* iU, float* d, float* rho0){
	/*creates S array*/
	float *S;
	S= (float*)malloc(N[2]*2*sizeof(float));
	
	int t;
	int i, j, a, b;
	for (t=0;t<N[2];t++){ 
		//puts initial zeros into S
		S[2*t]= 0;
		S[2*t+1]= 0;

		for (i=0;i<N[0];i++){
			for (j=0;j<N[0];j++){
				for (a=0;a<N[1];a++){
					for (b=0;b<N[1];b++){
						/*In rU (resp. iU) the relation between rU[n] and indeces i, j, a, b, t is n= i+ N[0]*(a+ N[1]*(j+ N[0]*(b+ N[1]*t))). 
						In d array the relation is n=i+ N[0]* a.*/

						/*real part*/
						S[2*t]+= d[b*N[0]+j]* d[a*N[0]+i]* rU[i+N[0]*(a+N[1]*(j+N[0]*(b+N[1]*t)))]*rho0[i];
						/*imaginary part*/
						S[2*t+1]+= d[b*N[0]+j]* d[a*N[0]+i]* iU[i+N[0]*(a+N[1]*(j+N[0]*(b+N[1]*t)))]*rho0[i];
					}
				}
			}
		}
	}

	/*returns S(t)*/
	return S;
}
/*
int main(){//tests nosecuda_polar1 with simple parameters
	float *S;
	//initialization of some parameters / these parameters should make S[0]=1+0i
	int N[3]={2,2,1};
	float rU[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	float iU[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float d[4]={1,1,1,1};
	float rho0[2]={1,1};

	S= (float*)nosecuda_polar_1(N, rU, iU, d, rho0);//computes S
	int i;
	for (i=0;i<1;i++)//writes S array
		printf("S(%i)=%f+i%f\n", i, S[2*i], S[2*i+1]);
	free(S);
}*/
