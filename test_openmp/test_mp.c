#include <stdlib.h>
#include <stdio.h>
#include <time.h>               /* chronometrage */
#include <string.h>             /* pour memset */
#include <math.h>
#include <sys/time.h>
#define NB_TIMES 10

#ifdef _OPENMP
#include <omp.h>
#endif

double my_gettimeofday()
{
        struct timeval tmp_time;
        gettimeofday(&tmp_time, NULL);
        return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int main(int argc, char **argv){
        double debut = my_gettimeofday();
	int N;
	if(argc == 2) {
		N = atoi(argv[1]);
        }
	//printf("N = %d\n", N);
	int T[N];
	int i;
	#pragma omp parallel for
	for(i = 0;i<N;i++){
		T[i] = 3;
	}
	//for(i = 0;i<N;i++){
	//	printf("%d ",T[i]);
	//}
	//printf("\n ");
	double fin = my_gettimeofday();

        fprintf(stdout, "For N=%d: total computation time (with gettimeofday()) : %g s\n", N, (fin - debut) / NB_TIMES);
}
