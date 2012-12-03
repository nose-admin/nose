noinst_LIBRARIES = libnosecuda.a
if CUDA 
  noinst_PROGRAMS = primtest.x ctest.x
  ctest_x_SOURCES = callCudaTest.cu cudaTest.cu
else
  noinst_PROGRAMS = primtest.x
endif

primtest_x_SOURCES = prim.F90 

if CUDA 
libnosecuda_a_SOURCES = primitive.cu
else 
libnosecuda_a_SOURCES = 
endif

primtest_x_LDADD = -L./ -lnosecuda $(CUDA_LIBS) $(CUDA_CFLAGS)
 
# suffixes
.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS)
	
.F90.o:
	${FC} ${DEFS} ${FCFLAGS} ${INCLUDES}  -c $<

if CUDA 	
CUDA_FFLAGS = -lGL -lGLU
CNDTN = libnosecuda.a prim.o
CMPLR = ${FC} -o primtest.x${EXEEXT} $(NVCFLAGS) prim.o -lnosecuda -L./ $(CUDA_LIBS) ${CUDA_CFLAGS} ${CUDA_FFLAGS} 
CTEST_O = callCudaTest.o cudaTest.o 
CTEST_CL = ${NVCC} -o ctest.x${EXEEXT} ${CTEST_O}
else
CNDTN = prim.o
CMPLR = ${FC} -o primtest.x${EXEEXT}  prim.o
CTEST_O =
CTEST_CL = @echo 'No ctest'
endif

callCudaTest.o: cudaTest.o

#primtest.x${EXEEXT}: ${CNDTN}
#	${CMPLR}

ctest.x${EXEEXT}: ${CTEST_O}
	${CTEST_CL}


	