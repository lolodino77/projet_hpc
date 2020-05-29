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
    int quotient = 5;
    int reste = 2;
    int sendcount = quotient;
    if(my_rank == 0){
        sendcount = reste;
    }
    int* A = malloc((p * quotient + 2) * sizeof(int)); // p * sendcount = 4 * 5 = 20
    int* A_part = malloc(sendcount * sizeof(int)); 
	// int A_part[sendcount];
    for(int i = 0;i < sendcount;i++){
            // printf("my rank = %d\n", my_rank);
            A_part[i] = my_rank;
            // printf("%d ", A_part[i]);
    }
    printf("\n");

    int recvcounts[p];
    for(int i = 1;i<p;i++){
        recvcounts[i] = sendcount;
    }

    int displs[p];
    displs[0] = 0;
    for(int i = 0;i<(p - 2);i++){
        displs[i + 2] = (i * sendcount) + 2;
    }
    printf("displs : \n");
    for(int i = 0; i < p; i ++){
        printf("%d\n", displs[i]);
    }

    printf("recvcounts = %d\n", recvcounts);
    printf("debut gather\n");        
    // MPI_Gatherv(A_part, sendcount, MPI_INT, A, recvcounts, MPI_INT, 0, MPI_COMM_WORLD)
    MPI_Allgatherv(A_part, sendcount, MPI_INT, A, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);    
    printf("fin gather\n");


    printf("my_rank = %d\n", my_rank);
    for(int i = 0;i < (p * sendcount + 2);i++){
        printf("%d ", A[i]);
    }
    printf("\n");   
    
    MPI_Finalize();

}
