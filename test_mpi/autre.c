#include <stdio.h>
#include <stdlib.h>

struct csr_matrix_t {
    int n;          // dimension
    int nz;         // number of non-zero entries
    int *Ap;        // row pointers
    int *Aj;        // column indices
    double *Ax;     // actual coefficient
};

int main(){
    struct csr_matrix_t *A = malloc(sizeof(*A)); 
    A->n = 9;
    A->nz = 10;
    struct csr_matrix_t *B = A;
    printf("B->n = %d, B->nz = %d\n", A->n, A->nz);
}