/*
	This file contains method nosecuda_secular_1(...) which returns response function S(t) in time. The secular approximation U_{bjai}(t)=\delta_ab \delta_ij U_{aiai}(t). Theory can be found in the white paper Absorption Spectrum at NOSE web.

	ATTENTION: This version is not yet optimized!

	At this state only N[0]==1 works as it should. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int BLOCKDIM=4;//size of thread blocks
/*
	Computes response functions S(t) in time using GPU.
	Parameters: 
		int N[3]: N[1] = number of levels in the ground state band
		          N[2] = number of levels in the excited state band
			  N[3] = number of the time steps
		float rU[N[0]*N[1]*N[2]] = real part of the coherence evolution superoperator in secular approximation. The superoperator has 4 indeces, the 5th dimension represents discrete time. 1st and 3rd indeces and 2nd and 4th indeces are the same.

		float iU[N[0]*N[1]*N[2]] = imaginary part of the coherence evolution superoperator in secular approximation. The superoperator has 4 indeces, the 5th dimension represents discrete time.
		float d[N[0]*N[1]] = transition dipole moment elements
		float rho0[N[0]] = diagonal elements of the equilibrium ground state density matrix
		float S[2*N[0]] = response function in time (both real and imaginary parts)
*/

void checkCUDAError(const char* msg);

__global__ void _evolveS(int N0, int N1, int N2, float *data, float *rS, float *iS){ //rU, iU don't contain a!=b || i!=j elements yent
	int t= (gridDim.x*blockIdx.y+ blockIdx.x)*blockDim.x*blockDim.y+ blockDim.x*threadIdx.y+threadIdx.x;//calculates time of current thread
	//if (t>=N2) return;//because 4x4 blocks are alocated, some threads at the end are not necessary
	int bIdx= gridDim.x*blockIdx.y+blockIdx.x;
	int tIdx= blockDim.x*threadIdx.y+threadIdx.x;

	//initialized temporary variables
	float im= 0;
	float re= 0;
	//int BLOCKSIZE= blockDim.x*blockDim.y;
	int size= N0*N1*N2/(gridDim.x*gridDim.y);//size of rU (resp iU) region supplied to one block (in floats)
	int _rUOffset= bIdx*size;//offset of _rU according to first address of rU
	int _iUOffset= bIdx*size;//offset of _iU according to first address of iU

	//shared memory duplicates of global memory variables
	extern __shared__ float _sData[];//sData[];//all dynamically alocated arrays in shared memory start at the same address. So I have to create offsets. Structure of sData is: rU[size],iU[size],d[N0*N1],rho0[N0]
	float *_rU, *_iU, *_d, *_rho0;
	float *rU, *iU, *d, *rho0;
	rU= &data[_rUOffset];//1*100*16
	iU= &data[N0*N1*N2+ _iUOffset];
	d= &data[2*N0*N1*N2];
	rho0= &data[2*N0*N1*N2+ N0*N1];

	_rU= _sData;
	_iU= &_rU[size];
	_d= &_iU[size];
	_rho0= &_d[N0*N1];
	
	//copy rU, iU regions and d and rho0 into shared memory / threads should be coalesced
	int i, a;//indeces i, a 
	for (i=tIdx;i<size;i+= 16){
		_rU[i]= rU[i];
		_iU[i]= iU[i];
		if (i<N0*N1) {
			_d[i]= d[i];
			if (i<N0){
				_rho0[i]= rho0[i];	
			}
		}
	}

	__syncthreads();//all values copied into the shared memory (by all threads)

	//calculates S in secular approximation
	for (i=0;i<N0;i++){
		for (a=0;a<N1;a++){
			/*real part*/
			re+= _d[a*N0+i]* _d[a*N0+i]* _rU[i+N0*(a+N1*tIdx)]*_rho0[i];
			/*imaginary part*/
			im+= _d[a*N0+i]* _d[a*N0+i]* _iU[i+N0*(a+N1*tIdx)]*_rho0[i];
		}
	}
	__syncthreads();//all threads have accomplished their computation

	//stores resulsts into S / threads are coalesced because t=threadIdx
	rS[t]= re;
	iS[t]= im;
}

float *nosecuda_secular_1(int *N, float *rU, float *iU, float *d, float *rho0){
	//input data serialization
	float *data, *_data;//this array will by passed to the kernel. It containes serialized rU,iU,d,rho0
	int size= N[0]*N[1]*N[2];//size of rU array in floats
	data= (float*)malloc(sizeof(float)*(2*size+N[0]*N[1]+N[0]));
	int i, t;
	for (i=0;i<size;i++) data[i]= rU[i];
	for (i=0;i<size;i++) data[i+size]= iU[i];
	for (i=0;i<N[0]*N[1];i++) data[i+2*size]= d[i]; 
	for (i=0;i<N[0];i++) data[i+ 2*size+ N[0]*N[1]]= rho0[i];

	/*creates S array and arrays to store data from GPU*/
	float *S, *iS, *rS, *_iS, *_rS;
	S= (float*)malloc(2*N[2]*sizeof(float));
	iS= (float*)malloc(N[2]*sizeof(float)); 
	rS= (float*)malloc(N[2]*sizeof(float)); 

	/*initializes device variables*/
	/*allocate _data array*/
	cudaMalloc((void**) &_data, (2*size+N[0]*N[1]+N[0])*sizeof(float)); 
	cudaMemcpy(_data, data, (2*size+N[0]*N[1]+N[0])*sizeof(float), cudaMemcpyHostToDevice);
	/*allocate S array*/
	cudaMalloc((void**) &_iS, N[2]*sizeof(float));//N[2] contains real 
	cudaMalloc((void**) &_rS, N[2]*sizeof(float));//N[2] contains real

	//calculates grid dimension
	dim3 dimBlock(BLOCKDIM, BLOCKDIM);

	int b= ceil(N[2]/(float)(BLOCKDIM*BLOCKDIM));//number of blocks needed
	i=1;//grid dimension in blocks
	while ((i+1)*(i+1)<=b) i+= 1;//find greatest square grid whose area is lesser than b
	dim3 dimGrid(i, i+ceil((b-i*i)/(float)i));

	//initialize testing f array
	/*float *_f, *f;
	cudaMalloc((void**) &_f, sizeof(float));*/

	/*kernel execution*/
	_evolveS<<<dimGrid, dimBlock, 16000>>>(N[0], N[1], N[2], _data, _rS, _iS);
	checkCUDAError("kernel invocation");

	/*copy result from the device to the host*/
	cudaMemcpy(rS, _rS, N[2]*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(iS, _iS, N[2]*sizeof(float), cudaMemcpyDeviceToHost);
	checkCUDAError("cpy");

	//copy data from rS and iS to S
	for (t=0;t<N[2];t++) {
		S[2*t]= rS[t];
		S[2*t+1]= iS[t];
	}
	/*f= (float*)malloc(sizeof(float));
	cudaMemcpy(f, _f, sizeof(float), cudaMemcpyDeviceToHost);
	printf("f=%f\n", *f);
	cudaFree(_f);*/
	

	/*free used memory*/
	cudaFree(_data);
	cudaFree(_iS);
	cudaFree(_rS);
	free(data);
	free(rS);
	free(iS);
	
	/*returns S(t)*/
	return S;//the function which calls this method should free S
}

float random01(){//returns random number between 0 and 1
	return rand()/(float)RAND_MAX;
}

int main(int argc, char **argv){//tests nosecuda_secular1 function
	//set dimensios of parameters
	int N[3]={1, 64, 64};
	float *d, *rU, *iU, *rho0;
	d= (float*)malloc(N[1]*N[0]*sizeof(float));
	rho0= (float*)malloc(N[0]*sizeof(float));
	rU= (float*)malloc(N[0]*N[1]*N[0]*N[1]*N[2]*sizeof(float));
	iU= (float*)malloc(N[0]*N[1]*N[0]*N[1]*N[2]*sizeof(float));
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
	//perform cpu computation
	float cpuS[2*N[2]];
	for (t=0;t<N[2];t++){//set rU, iU
		cpuS[2*t]= 0;
		cpuS[2*t+1]= 0;
		for (i=0;i<N[0];i++){//set d
			for (a=0;a<N[1];a++){
				cpuS[2*t]+= d[i+a*N[0]]* d[i+a*N[0]]* rU[i+N[0]*(a+N[1]*t)]* rho0[i];
				cpuS[2*t+1]+= d[i+a*N[0]]* d[i+a*N[0]]* iU[i+N[0]*(a+N[1]*t)]* rho0[i];
			}
		}
	}
	//perform gpu computation
	float *S= nosecuda_secular_1(N, rU, iU, d, rho0);
	//compare both cpu and gpu methods and print results
	for (t=0;t<N[2];t++){
		printf("S[%i]=%f+i%f;cpu: %f+i%f\n", t, S[2*t], S[2*t+1], cpuS[2*t], cpuS[2*t+1]);
	}
	//free memory
	free(d);
	free(rho0);
	free(rU);
	free(iU);
	free(S);
}

void checkCUDAError(const char *msg){//checks for cudaErrors
	cudaError_t err = cudaGetLastError();
	if( cudaSuccess != err) {   
		fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
//		exit(EXIT_FAILURE);
	}    
}

