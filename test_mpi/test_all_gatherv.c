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
    int n_part = quotient;
    if(my_rank == 0){
        n_part = reste;
    }
    int n = (p - 1) * n_part + 2;

    int* A = malloc(((p - 1) * quotient + reste) * sizeof(int)); // p * n_part = 4 * 5 = 20
    int* A_part = malloc(n_part * sizeof(int)); 
	// int A_part[n_part];
    printf("my_rank = %d\n", my_rank);
    printf("%d \n", (p - 1) * quotient + reste);
    for(int i = 0;i < n_part;i++){
            // printf("my rank = %d\n", my_rank);
            A_part[i] = my_rank;
            printf("%d ", A_part[i]);
    }
    printf("\n");

    int recvcounts[p];
    printf("recvcounts, p = %d :\n", p);
    recvcounts[0] = reste;
    printf("%d\n", recvcounts[0]);
    for(int i = 1;i<p;i++){
        recvcounts[i] = quotient;
        printf("%d\n", recvcounts[i]);
    }

    int displs[p - 1];
    displs[0] = 0;
    displs[1] = reste;
    for(int i = 2;i < p;i++){
        displs[i] = (i - 1) * quotient + reste; // displs[2] = displs[i] = 7 = 1*5 + 2 = (i-1)*5 + 2
    }
    printf("displs : \n");
    for(int i = 0; i < (p - 1); i ++){
        printf("%d ", displs[i]);
    }
    printf("\n");

    printf("debut gather\n");        
    // MPI_Gatherv(A_part, n_part, MPI_INT, A, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD)
    MPI_Allgatherv(A_part, n_part, MPI_INT, A, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);    
    printf("fin gather\n");


    printf("my_rank = %d\n", my_rank);
    for(int i = 0;i < n;i++){
        printf("%d ", A[i]);
    }
    printf("\n");   
    
    MPI_Finalize();

}
