#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// include <cutil.h>
// include <cuda.h>
#include <builtin_types.h>

#define INDX4(n, m) n + 4*m
#define INDX(a, b, c, d) a + b 

// device configuration

#define MULTIPROC_NR        4
#define MULTIPROC_MULTIPLE  2

#define BLOCK_SIZE   MULTIPROC_MULTIPLE*MULTIPROC_NR



/*

     Callable routines 

*/
extern "C" void primitive_gpu_(int*, int*, float* , float* , float* , float* , float* , float* );


extern "C" void test_int_(int* );
extern "C" void test_float_(float*);
extern "C" void test_ptr_(float*, int*);

/*

     Internal declarations

*/
__device__ void diffs(float*, float*, float*, float*, float*, float*, int, int);

/*

     This function performs the calculation

*/
__global__ void global_primitive(float* x, float* y, float* w,
                                 float* d, float* pr, float* pi,
                                 int Nx, int Nw) {

    // here the diffs must be summed in some sensible way

	diffs(x, y, w, d, pr, pi, Nx, Nw);


}



/*

    Solution of the 4x4 system of equations with Gauss elimination

*/
__device__ void gaussel(float* A, float* b, float* x) {

    int m;
    int n;

    __shared__ float a1[3][3];
    __shared__ float a2[2][2];
    __shared__ float b1[3];
    __shared__ float b2[2];

    for (m = 0; m < 3; m++) {
        for (n = 0; n < 3; n++) {
            a1[m][n] = A[INDX4(m,n)] - A[INDX4(m,0)]*A[INDX4(1,n)]/A[0];
        }
        b1[m] = b[m] - A[INDX4(m,0)]*b[0]/A[0];
    }
    
    for (m = 0; m < 2; m++) {
        for (n = 0; n < 2; n++) {
            a2[m][n] = a1[m][n] - a1[m][0]*a1[0][n]/a1[0][0];
        }
        b2[m] = b1[m] - a1[m][0]*b1[0]/a1[0][0];
    }

    x[3] = (a2[0][0]*b2[1] - a2[1][0]*b2[0])/(a2[0][0]*a2[1][1]-a2[1][0]*a2[0][1]);
    x[2] = (b2[0]/a2[0][0]) - (a2[0][1]/a2[0][0])*x[3];
    x[1] = (b1[0]/a1[0][0]) - (a1[0][1]/a1[0][0])*x[2] - (a1[0][2]/a1[0][0])*x[3];
    x[0] = (b[0] - A[INDX4(0,1)]*x[1] - A[INDX4(0,2)]*x[2] - A[INDX4(0,3)]*x[3])/A[0];

}


/*

    Subroutine for spline integration

*/
__device__ void get_abcd(float x1, float x2, float y1, float y2,
                         float d1, float d2, float* v) {

    __shared__ float M[16];
    __shared__ float b[4];

    M[0] = 1.0f;
    M[INDX4(1,0)] = 1.0f;
    M[INDX4(2,1)] = 1.0f;
    M[INDX4(3,1)] = 1.0f;
    M[INDX4(0,1)] = x1;
    M[INDX4(0,2)] = M[INDX4(0,1)]*x1;
    M[INDX4(0,3)] = M[INDX4(0,2)]*x1;
    M[INDX4(1,1)] = x2;
    M[INDX4(1,2)] = M[INDX4(1,1)]*x2;
    M[INDX4(1,3)] = M[INDX4(1,2)]*x2;
    M[INDX4(2,2)] = 2.0f*x1;
    M[INDX4(2,3)] = 3.0f*x1*x1;
    M[INDX4(3,2)] = 2.0f*x2;
    M[INDX4(3,3)] = 3.0f*x2*x2;

    b[0] = y1;
    b[1] = y2;
    b[2] = d1;
    b[3] = d2;

    gaussel(M, b, v);

} 


/*

    Multiplication of two complex numbers with and without constructor

*/
__device__ float2 m2c(float2 a, float2 b) {

    return make_float2(a.x*b.x - a.y*b.y,a.x*b.y + a.y*b.x);

}


/*

    Real by complex multiplication with and without constructor

*/
__device__ float2 mrc(float r, float2 a) {

    return make_float2(r*a.x,r*a.y);

}


/*

    Complex addition without constructor

*/
__device__ float2 d2c(float2 a, float2 b) {

    return make_float2(a.x - b.x,a.y - b.y);

}

/*

    Calculates integrals between interpolation points

*/
__device__ void diffs(float* x, float* y, float* w,
                      float* d, float* pr, float* pi,
                      int Nx, int Nw) {


    float x1, x2, y1, y2, d1, d2;
    float2 BeB, AeA, B2eB, A2eA, dI0, dI1, dI2, dI3;
    float2 I0, I1, I2, I3, eA, eB, iom; 
    float2 aux;
    
    __shared__ float A[4];

    int i;
    int k;

    // chose i and k according to the thread id

    i = 0;
    k = 0;

    // here I load from global to shared memory
    x1 = x[i];
    x2 = x[i+1];
    y1 = y[i];
    y2 = y[i+1];
    d1 = d[i];
    d2 = d[i+1];

    // calculation

    get_abcd(x1,x2,y1,y2,d1,d2,A);

    eA  = make_float2(cosf(w[k]*x1),sinf(w[k]*x1));
    eB  = make_float2(cosf(w[k]*x2),sinf(w[k]*x2));
    iom = make_float2(0.0f, -1.0f/w[k]);
 
    BeB = mrc(x2,eB);
    AeA = mrc(x1,eA);
    B2eB = mrc(x2,BeB);
    A2eA = mrc(x1,AeA);
    dI0 = make_float2(eB.x-eA.x,eB.y-eA.y);
    dI1 = make_float2(BeB.x-AeA.x,BeB.y-AeA.y);
    dI2 = make_float2(B2eB.x-A2eA.x,B2eB.y-A2eA.y);
    aux = mrc(x2,B2eB);
    dI3 = make_float2(x2*B2eB.x,x2*B2eB.y);
    dI3.x = dI3.x - aux.x;
    dI3.y = dI3.y - aux.y;

    I0 = m2c(iom,dI0);
    aux = d2c(dI1,I0);
    I1 = m2c(iom,aux);
    aux = mrc(2.0f,I1);
    aux = d2c(dI2,aux);
    I2 = m2c(iom,aux);
    aux = mrc(3.0f,I2);
    aux = d2c(dI2,aux);
    I3 = m2c(iom,aux);

    // returning to global memory

    pr[INDX(i+1,k,Nk,Nw)] = A[0]*I0.x + A[1]*I1.x + A[2]*I2.x + A[3]*I3.x;
    pi[INDX(i+1,k,Nk,Nw)] = A[0]*I0.y + A[1]*I1.y + A[2]*I2.y + A[3]*I3.y;

}


/*

   Function callable from host 

*/
void primitive_gpu_(int* Nxi, int* Nwi, float* x, float* y, float* w, float* d, float* pr, float* pi) {

   int nThread;
   int gridSize;
   int Nx, Nw, i;
   struct cudaDeviceProp dprop;


   Nx = *Nxi;
   Nw = *Nwi;

   cudaGetDeviceCount(&i);
   cudaGetDeviceProperties(&dprop,0);
   

   printf("primitive_gpu ...%i - %i \n",i,dprop.multiProcessorCount);

   printf("   %i %i \n",Nx,Nw);

   /* allocate device memory */

   float* d_x;
   float* d_y;
   float* d_w;
   float* d_d;
   float* d_pr;
   float* d_pi;


   // number of threads and the size of the grid
   nThread = Nx*Nw; 
   gridSize = int(nThread/ int(BLOCK_SIZE));

   if (gridSize*BLOCK_SIZE != nThread) gridSize++;

   printf("   %i %i \n",Nx,Nw);
   printf("   block size = %i \n",BLOCK_SIZE);
   printf("   grid size = %i \n",gridSize);
   printf("   number of threads = %i \n", nThread);
   
   printf("%f \n",x[0]);
   printf("%f \n",x[1]);


   // copy into device memory
   cudaMalloc((void**)&d_x,Nx);
   cudaMemcpy(d_x,x,Nx,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&d_y,Nx);
   cudaMemcpy(d_y,y,Nx,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&d_d,Nx);
   cudaMemcpy(d_d,d,Nx,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&d_pr,Nx*Nx);
   cudaMalloc((void**)&d_pi,Nx*Nx);
   cudaMalloc((void**)&d_w,Nw);
   cudaMemcpy(d_w,w,Nw,cudaMemcpyHostToDevice);


   // call GPU
   dim3 dimBlock(BLOCK_SIZE);
   dim3 dimGrid(gridSize);
   global_primitive<<<dimBlock, dimGrid>>>(d_x, d_y, d_w, d_d,
                                           d_pr, d_pi, Nx, Nw);


   cudaMemcpy(pi, d_pi,Nx*Nw,cudaMemcpyDeviceToHost);
   cudaMemcpy(pr, d_pr,Nx*Nw,cudaMemcpyDeviceToHost);

   // copy data to the host

   cudaFree(d_x);
   cudaFree(d_y);
   cudaFree(d_d);
   cudaFree(d_w);
   cudaFree(d_pr);
   cudaFree(d_pi);

   printf("... primitive_gpu finished\n");

}

void test_int_(int* N) {

  int Ni;

  Ni = *N;

  printf("%i \n", Ni);
  
  *N = 2;

}

void test_float_(float* r) {

  printf("%f \n",r[0]);

}

void test_ptr_(float* x, int* N) {
	
	int Ni;
	Ni = *N;
	x[Ni-1] = 10.001;
	
}
