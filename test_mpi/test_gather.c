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
    int sendcount = 5;
    int* A = malloc(p * sendcount * sizeof(int)); // p * sendcount = 4 * 5 = 20
    int* A_part = malloc(sendcount * sizeof(int)); 
	// int A_part[sendcount];
    for(int i = 0;i < sendcount;i++){
            // printf("my rank = %d\n", my_rank);
            A_part[i] = my_rank;
            // printf("%d ", A_part[i]);
    }
    printf("\n");

    int recvcount = sendcount;
    printf("..recvcount = %d\n", recvcount);
    printf("debut gather\n");        
    MPI_Gather(A_part, sendcount, MPI_INT, A, recvcount, MPI_INT, root, MPI_COMM_WORLD);    
    printf("fin gather\n");

    if(my_rank == 0){
        printf("my_rank = %d\n", my_rank);
        for(int i = 0;i < p*recvcount;i++){
            printf("%d ", A[i]);
        }   
    }

    MPI_Finalize();

}
