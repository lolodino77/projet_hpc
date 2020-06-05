#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
	double b = 12.455;
	double a = 2.35;
	double c = fmod(b,a);
	printf("mod b/a = %lf\n", c);
}
