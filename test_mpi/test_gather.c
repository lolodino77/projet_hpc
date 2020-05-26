#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//mpicc -o test_mpi test_mpi.c 
//mpirun -n 5 -hostfile hostfile --map-by node ./test_mpi

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, my_rank, world_size);

    int root = 0;
    int p = 4;
    int sendcount = 5;
    int* A = malloc(p * sendcount * sizeof(int)); // p * sendcount = 4 * 5 = 20
    int* A_part = malloc(sendcount * sizeof(int)); 

    if(my_rank != 0){
        for(int i = 0;i < sendcount;i++){
            // printf("my rank = %d\n", my_rank);
            A_part[i] = my_rank;
            printf("%d ", A_part[i]);
        }
        printf("\n");
    }
    if(my_rank == 0){
        printf("my_rank = %d\n", my_rank);
        printf("debut gather\n");        
        int recvcount = p * sendcount;
        printf("recvcount = %d\n", recvcount);
        MPI_Gather(A_part, sendcount, MPI_INT, A, recvcount, MPI_INT, root, MPI_COMM_WORLD);
        printf("fin gather\n");
        for(int i = 0;i < recvcount;i++){
            printf("%d ", A[i]);
        }   
    }

    MPI_Finalize();

}
