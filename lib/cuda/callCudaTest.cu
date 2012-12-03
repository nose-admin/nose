/*
	Calls cudaTest function from cudaTest.cu file and returns it's return value
*/
#include <stdio.h>

extern int cudaTest();

int main(){
	int i;
	i = cudaTest();
	
	printf("%i \n", i);
	
	return i;
}

