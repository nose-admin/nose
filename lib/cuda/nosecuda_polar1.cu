/*
	This file contains method nosecuda_polar_1(...) which returns response function S(t) in time. Theory can be found in the white paper Absorption Spectrum at NOSE web.

	ATTENTION: This version is not yet optimized!
*/

#include <stdio.h>
#include <stdlib.h>

int BLOCKDIM=16;//size of thread blocks
/*
	Computes response functions S(t) in time using GPU.
	Parameters: 
		int N[3]: N[1] = number of lovels in the ground state band
		          N[2] = number of levels in the excited state band
			  N[3] = number of the time steps
		float rU[N[0]*N[1]*N[0]*N[1]*N[2]] = real part of the coherence evolution superoperator in time. The superoperator has 4 indeces, the 5th dimension represents discrete time.

		float iU[N[0]*N[1]*N[0]*N[1]*N[2]] = imaginary part of the coherence evolution superoperator. The superoperator has 4 indeces, the 5th dimension represents discrete time.
		float d[N[0]*N[1]] = transition dipole moment elements
		float rho0[N[0]] = diagonal elements of the equilibrium ground state density matrix
		float S[2*N[0]] = response function in time (both real and imaginary parts)
*/
__global__ void _evolveS(int* N, float* rU, float* iU, float* d, float* rho0, float *S, int BLOCKDIM){ 
	int t= (gridDim.x*blockIdx.y+ blockIdx.x)*BLOCKDIM*BLOCKDIM+ blockDim.x*threadIdx.y+threadIdx.x;//calculates time of current thread
	if (t>=N[2]) return;//because 16x16 blocks are alocated, some threads at the end are not necessary
	//initialized temporary variables
	float im= 0;
	float re= 0;

	int i, j, a, b;
	for (i=0;i<N[0];i++){
		for (j=0;j<N[0];j++){
			for (a=0;a<N[1];a++){
				for (b=0;b<N[1];b++){
					/*In rU (resp. iU) the relation between rU[n] and indeces i, j, a, b, t is n= i+ N[0]*(a+ N[1]*(j+ N[0]*(b+ N[1]*t))). 
					In d array the relation is n=i+ N[0]* a.*/

					/*real part*/
					re+= d[b*N[0]+j]* d[a*N[0]+i]* rU[i+N[0]*(a+N[1]*(j+N[0]*(b+N[1]*t)))]*rho0[i];
					/*imaginary part*/
					im+= d[b*N[0]+j]* d[a*N[0]+i]* iU[i+N[0]*(a+N[1]*(j+N[0]*(b+N[1]*t)))]*rho0[i];
				}
			}
		}
	}
	//stores resulsts into S
	S[2*t]= re;
	S[2*t+1]= im;
}

float *nosecuda_polar_1(int* N, float* rU, float* iU, float* d, float* rho0){
	/*creates S array*/
	float *S;
	S= (float*)malloc(2*N[2]*sizeof(float));

	/*initializes device variables*/
	float *_rU, *_iU, *_d, *_rho0, *_S;
	int *_N;
	/*allocate N array*/
	cudaMalloc((void**) &_N, 3*sizeof(int));/*length of N is 3*/
	cudaMemcpy(_N, N, 3*sizeof(float), cudaMemcpyHostToDevice);
	/*allocate rU array*/
	cudaMalloc((void**) &_rU, N[0]*N[1]*N[0]*N[1]*N[2]*sizeof(float)); 
	cudaMemcpy(_rU, rU, N[0]*N[1]*N[0]*N[1]*N[2]*sizeof(float), cudaMemcpyHostToDevice);
	/*allocate iU array*/
	cudaMalloc((void**) &_iU, N[0]*N[1]*N[0]*N[1]*N[2]*sizeof(float)); 
	cudaMemcpy(_iU, iU, N[0]*N[1]*N[0]*N[1]*N[2]*sizeof(float), cudaMemcpyHostToDevice);
	/*allocate d array*/
	cudaMalloc((void**) &_d, N[0]*N[1]*sizeof(float));
	cudaMemcpy(_d, d, N[0]*N[1]*sizeof(float), cudaMemcpyHostToDevice);
	/*allocate rho0 array*/
	cudaMalloc((void**) &_rho0, N[0]*sizeof(float));
	cudaMemcpy(_rho0, rho0, N[0]*sizeof(float), cudaMemcpyHostToDevice);
	/*allocate S array*/
	cudaMalloc((void**) &_S, 2*N[2]*sizeof(float));//2*N[2] contains real and imaginarypart

	//calculates grid dimension
	dim3 dimBlock(BLOCKDIM, BLOCKDIM);

	int b= N[2]/(BLOCKDIM*BLOCKDIM)+1;//number of blocks
	int i=1;//grid dimension in blocks
	while ((i+1)*(i+1)<b) i+= 1;//find greatest square grid whose area is lesser than b
	dim3 dimGrid(i, (b-i*i)/i+ 1);

	/*kernel execution*/
	_evolveS<<<dimGrid, dimBlock>>>(_N, _rU, _iU, _d, _rho0, _S, BLOCKDIM);

	/*copy result from the device to the host*/
	cudaMemcpy(S, _S, N[2]*2*sizeof(float), cudaMemcpyDeviceToHost);

	/*free used memory*/
	cudaFree(_N);
	cudaFree(_rU);
	cudaFree(_iU);
	cudaFree(_d);
	cudaFree(_rho0);
	cudaFree(_S);
	
	/*returns S(t)*/
	return S;//the function which calls this method should free S
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

	S= nosecuda_polar_1(N, rU, iU, d, rho0);//computes S
	int i;
	for (i=0;i<1;i++)//writes S array
		printf("S(%i)=%f+i%f\n", i, S[2*i], S[2*i+1]);
	free(S);
}*/
