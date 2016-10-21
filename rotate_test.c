
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double real;

void vector_project(real *input, real angle, real *output){
	
	angle = angle*3.1415956/180.;
	output[0]=input[0]*cos(angle)+input[1]*sin(angle);
	output[1]=-input[0]*sin(angle)+input[1]*cos(angle);
	
}

int main(){
	
	real output[2];
	real input[2] = {1,1};
	
	vector_project(input,180.,output);
	printf("%f %f\n",output[0], output[1]);

	return 0;
}