#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
	double b = 12.45577823784;
	double* integer = malloc(sizeof(double));
	double c = modf(b, integer);
	printf("partie decimale de %lf = %lf\n", b, c);
	free(integer);
}
