#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

struct csr_matrix_t {
    int n;          // dimension
    int nz;         // number of non-zero entries
    int *Ap;        // row pointers
    int *Aj;        // column indices
    double *Ax;     // actual coefficient
};

int main(int argc, char** argv) {
    int n = 0;
    // MPI_Aint displacements[2] = {};
    // int block_lengths[2] = {};
    // MPI_Datatype types[2] = {};
    // MPI_Datatype MPI_csr_matrix_t;    
    // struct csr_matrix_t *A = malloc(sizeof(*A)); 

    // Initialize the MPI environment
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

    // printf("avant broadcast : A->n = %d, A->nz = %d\n", A->n, A->nz);

    if(my_rank == 0){
        n = 10000;
        A->n = 9;
        A->nz = 10;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A, 1, MPI_csr_matrix_t, 0, MPI_COMM_WORLD);

    // printf("apres broadcast : A->n = %d, A->nz = %d\n", A->n, A->nz);
    printf("n = %d\n", n);
    // Finalize the MPI environment.
    MPI_Finalize();
}
