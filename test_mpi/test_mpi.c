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
    // int n = 0;
    // int count = 5;
    // MPI_Aint displacements[5] = {};
    // int array_of_blocklengths[] = {1, 1};
    // MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
    // MPI_Datatype MPI_csr_matrix_t;    
    // MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
    //                     array_of_types, &MPI_csr_matrix_t);

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

    struct csr_matrix_t *A = malloc(sizeof(*A)); 
    if(my_rank == 0){
        A->n = 10;
        A->nz = 10;
        A->Ap = malloc(5*sizeof(int));
        A->Aj = malloc(5*sizeof(int));
        A->Ax = malloc(5*sizeof(double));
        for(int i = 0;i < 5;i++){
            A->Ap = 10;
            A->Aj = 10;
            A->Ax = 10;    
        }
    }
    printf("avant broadcast (%d) : A->n = %d, A->nz = %d\n", my_rank ,A->n, A->nz);

    MPI_Bcast(&A->n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A->nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A->Ap, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A->Aj, 5, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A->Ax, 5, MPI_INT, 0, MPI_COMM_WORLD);

    printf("apres broadcast : A->n = %d, A->nz = %d\n", A->n, A->nz);
    printf("n = %d\n", n);
    // Finalize the MPI environment.
    MPI_Finalize();
}
