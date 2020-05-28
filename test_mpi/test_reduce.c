#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//mpicc -o test_mpi test_mpi.c 
//mpirun -n 5 -hostfile hostfile --map-by node ./test_mpi

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);

	printf("size = %d\n", p);

    // Get the rank of the process
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, my_rank, p);

    int root = 0;


    int *rz = malloc(sizeof(int));
	int *rz_part = malloc(sizeof(int));
    *rz_part = 10;
    *rz = 0;
	// printf("rz = %d\n", *rz);
 //    MPI_Reduce(&rz_part, &rz, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
 //    printf("after reduce\n");
 //    if(my_rank == 0){
 //        printf("rz = %d\n", *rz);
 //    }

	MPI_Finalize();
}
