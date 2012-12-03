/*
	You can call cudaTest() from a C code. This function:
	returns 1 if cuda works well
		0 if no GPU is present
		-1 if an error in computation occured
*/

#include <stdio.h>
const float t = 2.5;//multiplicative parameter
const int N= 1024;//matrixes are NxN
const int blockSize= 16;//blocks are 16x16


/*
	cuda kernel
	Each thread multiplies one element of the matrix by the t parameter
*/
__global__ void _multiplyArray(float* a, float t, int N){
    int x = blockIdx.x* blockDim.x+ threadIdx.x;
    int y = blockIdx.y* blockDim.y+ threadIdx.y;
    
    int i= x+ N*y;
    
    a[i]= t*a[i];
}

int cudaTest(){
    /*tests count of present cuda devices
    This doesn't work yet. I dont't know how to determine there is certainly the emulation mode in progress.*/
    int i;
    /*determines whether a cuda capable device is used*/
    cudaGetDevice(&i);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    if ((prop.major< 1)||(prop.major> 1)) {/*in Feb 09 all cuda enabled devices have prop.major==1*/
    	cudaGetDeviceCount(&i);
	if (i< 2) {/*no other device in the system*/
		return 0;/*no GPU*/
	} else {
		printf("You are not using the cuda capable device!");
		return 0;
	}
    }

    /*host array initialization*/
    float f[N*N];
    for (i=0;i<N*N;i++) f[i]= 1;


    const int size= N*N*sizeof(float);
    
    /*initializes the device array and fills it with the content of f*/
    float *_f;
    cudaMalloc((void**) &_f, size);
    cudaMemcpy(_f, f, size, cudaMemcpyHostToDevice);
    
    /*executes the kernel*/
    dim3 dimBlock(blockSize, blockSize);
    dim3 dimGrid(N/dimBlock.x, N/dimBlock.y);
    _multiplyArray<<<dimGrid, dimBlock>>>(_f, t, N);
    
    /*copies data back to the host array f and cleans the memory up*/
    cudaMemcpy(f, _f, size, cudaMemcpyDeviceToHost);
    cudaFree(_f);
    
    /*tests the result*/
    for (i=0;i<N*N;i++) 
    	if (f[i]!= t) return -1; /*computation failed*/
    
    return 1; /*everything all right*/
}

