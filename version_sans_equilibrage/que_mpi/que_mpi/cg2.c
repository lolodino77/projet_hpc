/* 
 * Parallel implementation of the Conjugate Gradient Method.
 *
 * Authors : Laurent Dang-Vu et Thomas Genin
 *
 * v1.02 (2020-04-3)
 *
 * CHANGE LOG:
 *    v1.01 : fix a minor printing bug in load_mm (incorrect CSR matrix size)
 *    v1.02 : use https instead of http in "PRO-TIP"
 *  
 * USAGE: 
 * 	$ ./cg --matrix bcsstk13.mtx                # loading matrix from file
 *      $ ./cg --matrix bcsstk13.mtx > /dev/null    # ignoring solution
 *	$ ./cg < bcsstk13.mtx > /dev/null           # loading matrix from stdin
 *      $  zcat matrix.mtx.gz | ./cg                # loading gziped matrix from
 *      $ ./cg --matrix bcsstk13.mtx --seed 42      # changing right-hand side
 *      $ ./cg --no-check < bcsstk13.mtx            # no safety check
 *
 * PRO-TIP :
 *      # downloading and uncompressing the matrix on the fly
 *	$ curl --silent https://hpc.fil.cool/matrix/bcsstk13.mtx.gz | zcat | ./cg
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>		/* chronometrage */
#include <string.h>		/* pour memset */
#include <err.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <mpi.h>

#include "mmio.h"

#define THRESHOLD 1e-8		// maximum tolerance threshold

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

struct csr_matrix_t {
	int n;			// dimension
	int nz;			// number of non-zero entries
	int *Ap;		// row pointers
	int *Aj;		// column indices
	double *Ax;		// actual coefficient
};

/*************************** Utility functions ********************************/

/* Seconds (wall-clock time) since an arbitrary point in the past */
double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}//calcule tps qui secoule

/* Pseudo-random function to initialize b (rumors says it comes from the NSA) */
#define ROR(x, r) ((x >> r) | (x << (64 - r)))
#define ROL(x, r) ((x << r) | (x >> (64 - r)))
#define R(x, y, k) (x = ROR(x, 8), x += y, x ^= k, y = ROL(y, 3), y ^= x)
double PRF(int i, unsigned long long seed)
{
	unsigned long long y = i, x = 0xBaadCafe, b = 0xDeadBeef, a = seed;
	R(x, y, b);
	for (int i = 0; i < 31; i++) {
		R(a, b, i);
		R(x, y, b);
	}
	x += i;
	union { double d; unsigned long long l;	} res;
	res.l = ((x << 2) >> 2) | (((1 << 10) - 1ll) << 52);
	return 2 * (res.d - 1.5);
}

/*************************** Matrix IO ****************************************/

/* Load MatrixMarket sparse symetric matrix from the file descriptor f */
struct csr_matrix_t *load_mm(FILE * f, int *nnz2)//construct
{
	MM_typecode matcode;
	int n, m, nnz;

	/* -------- STEP 1 : load the matrix in COOrdinate format */
	double start = wtime();

	/* read the header, check format */
	if (mm_read_banner(f, &matcode) != 0)
		errx(1, "Could not process Matrix Market banner.\n");
	if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
		errx(1, "Matrix Market type: [%s] not supported (only sparse matrices are OK)", mm_typecode_to_str(matcode));
	if (!mm_is_symmetric(matcode) || !mm_is_real(matcode))
		errx(1, "Matrix type [%s] not supported (only real symmetric are OK)", mm_typecode_to_str(matcode));
	if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
		errx(1, "Cannot read matrix size");
	*nnz2 = nnz;
	fprintf(stderr, "[IO] Loading [%s] %d x %d with %d nz in triplet format\n", mm_typecode_to_str(matcode), n, n, nnz);
	fprintf(stderr, "     ---> for this, I will allocate %.1f MByte\n", 1e-6 * (40.0 * nnz + 8.0 * n));

	/* Allocate memory for the COOrdinate representation of the matrix (lower-triangle only) */
	int *Ti = malloc(nnz * sizeof(*Ti));
	int *Tj = malloc(nnz * sizeof(*Tj));
	double *Tx = malloc(nnz * sizeof(*Tx));
	if (Ti == NULL || Tj == NULL || Tx == NULL)
		err(1, "Cannot allocate (triplet) sparse matrix");

	/* Parse and load actual entries */
	for (int u = 0; u < nnz; u++) {
		int i, j;
		double x;
		if (3 != fscanf(f, "%d %d %lg\n", &i, &j, &x))
			errx(1, "parse error entry %d\n", u);
		Ti[u] = i - 1;	/* MatrixMarket is 1-based */
		Tj[u] = j - 1;
		/*
		 * Uncomment this to check input (but it slows reading)
		 * if (i < 1 || i > n || j < 1 || j > i)
		 *	errx(2, "invalid entry %d : %d %d\n", u, i, j); 
		 */
		Tx[u] = x;
	}

	double stop = wtime();
	fprintf(stderr, "     ---> loaded in %.1fs\n", stop - start);

	/* -------- STEP 2: Convert to CSR (compressed sparse row) representation ----- */
	start = wtime();

	/* allocate CSR matrix */
	struct csr_matrix_t *A = malloc(sizeof(*A));
	if (A == NULL)
		err(1, "malloc failed");
	int *w = malloc((n + 1) * sizeof(*w));
	int *Ap = malloc((n + 1) * sizeof(*Ap));
	int *Aj = malloc(2 * nnz * sizeof(*Ap));
	double *Ax = malloc(2 * nnz * sizeof(*Ax));
	if (w == NULL || Ap == NULL || Aj == NULL || Ax == NULL)
		err(1, "Cannot allocate (CSR) sparse matrix");

	/* the following is essentially a bucket sort */

	/* Count the number of entries in each row */
	for (int i = 0; i < n; i++)
		w[i] = 0;
	for (int u = 0; u < nnz; u++) {
		int i = Ti[u];
		int j = Tj[u];
		w[i]++;
		if (i != j)	/* the file contains only the lower triangular part */
			w[j]++;
	}

	/* Compute row pointers (prefix-sum) */
	int sum = 0;
	for (int i = 0; i < n; i++) {
		Ap[i] = sum;
		sum += w[i];
		w[i] = Ap[i];
	}
	Ap[n] = sum;

	/* Dispatch entries in the right rows */
	for (int u = 0; u < nnz; u++) {
		int i = Ti[u];
		int j = Tj[u];
		double x = Tx[u];
		Aj[w[i]] = j;
		Ax[w[i]] = x;
		w[i]++;
		if (i != j) {	/* off-diagonal entries are duplicated */
			Aj[w[j]] = i;
			Ax[w[j]] = x;
			w[j]++;
		}
	}

	/* release COOrdinate representation */
	free(w);
	free(Ti);
	free(Tj);
	free(Tx);
	stop = wtime();
	fprintf(stderr, "     ---> converted to CSR format in %.1fs\n", stop - start);
	fprintf(stderr, "     ---> CSR matrix size = %.1fMbyte\n", 1e-6 * (24. * nnz + 4. * n));

	A->n = n;
	A->nz = sum;
	A->Ap = Ap;
	A->Aj = Aj;
	A->Ax = Ax;
	return A;
}

/*************************** Matrix accessors *********************************/

/* Copy the diagonal of A into the vector d. */
void extract_diagonal(const struct csr_matrix_t *A, double *d){
	int n = A->n;
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = 0; i < n; i++) {
		d[i] = 0.0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++)
			if (i == Aj[u])
				d[i] += Ax[u];
	}
}

void extract_diagonal_part(const struct csr_matrix_t *A, double *d, int n_part, int i_ini)
{
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = i_ini; i < i_ini + n_part; i++) {
		d[i] = 0.0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++)
			if (i == Aj[u])
				d[i] += Ax[u];
	}
}

/* Matrix-vector product (with A in CSR format) : y = Ax */
void sp_gemv(const struct csr_matrix_t *A, const double *x, double *y)
{
	int n = A->n;
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = 0; i < n; i++) {
		y[i] = 0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++) {
			int j = Aj[u];
			double A_ij = Ax[u];
			y[i] += A_ij * x[j];
		}
	}
}

void sp_gemv_part(const struct csr_matrix_t *A, const double *x, double *y, int n_part, int i_ini)
{
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = i_ini; i < i_ini + n_part; i++) {
		y[i - i_ini] = 0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++) {
			int j = Aj[u]; //j = Aj[Ap[i]]
			double A_ij = Ax[u];
			y[i - i_ini] += A_ij * x[j];
		}
	}
	// printf("produit matriciel reussi\n");
}

/*************************** Vector operations ********************************/

/* dot product */
double dot(const int n, const double *x, const double *y)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;
}

/* euclidean norm (a.k.a 2-norm) */
double norm(const int n, const double *x)
{
	return sqrt(dot(n, x, x));
}

double dot_part( const double *x, const double *y, int i_ini, int n_part)
{
	double sum = 0.0;
	for (int i = i_ini; i < i_ini + n_part; i++)
		sum += x[i] * y[i];
	return sum;
}

/* euclidean norm (a.k.a 2-norm) */
double norm_part( const double *x, int i_ini, int n_part)
{
	return sqrt(dot_part(x, x, i_ini, n_part));
}

/******************************* main program *********************************/

/* options descriptor */
struct option longopts[6] = {
	{"seed", required_argument, NULL, 's'},
	{"rhs", required_argument, NULL, 'r'},
	{"matrix", required_argument, NULL, 'm'},
	{"solution", required_argument, NULL, 'o'},
	{"no-check", no_argument, NULL, 'c'},
	{NULL, 0, NULL, 0}
};

int main(int argc, char **argv)
{
	/* Initialisation de variables */
	int my_rank; //rank of the process
	int P; //number of process
	int i_block = 0; //numero du bloc courant en train d'etre calcul??
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
	MPI_Status status;
	int dest;
	int bTmp = 0; //numero du dernier bloc du vecteur x qui vient d'??tre calcul??
	int idTmp; //num??ro du processus dont le ma??tre vient de recevoir le travail
	int tagFin;
	double debut, fin;
	enum tagType {INDICE, TRAITEMENT, STOP, DOT_RZ, DOT_PQ, MATPROD};
	debut = my_gettimeofday();

	/* Parse command-line options */
	long long seed = 0;
	char *rhs_filename = NULL;
	char *matrix_filename = NULL;
	char *solution_filename = NULL;
	int safety_check = 1;
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 's':
			seed = atoll(optarg);
			break;
		case 'r':
			rhs_filename = optarg;
			break;
		case 'm':
			matrix_filename = optarg;
			break;
		case 'o':
			solution_filename = optarg;
			break;
		case 'c':
			safety_check = 0;
			break;
		default:
			errx(1, "Unknown option");
		}
	}

	// printf("hello i am process %s number %d\n", processor_name, my_rank);

/* Broadcast de la matrice A */
	int n = 0;
	int nnz = 0;
	int *nnz2 = malloc(sizeof(int));
	struct csr_matrix_t *A = malloc(sizeof(*A));

	/* Load the matrix */
	if(my_rank == 0){
		FILE *f_mat = stdin;
		if (matrix_filename) {
			f_mat = fopen(matrix_filename, "r");
			if (f_mat == NULL)
				err(1, "cannot matrix file %s", matrix_filename);
		}
		A = load_mm(f_mat, nnz2);
		n = A->n;
		nnz = *nnz2;
	}
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(my_rank != 0){
		A->n = n;
		A->Ap = malloc((n+1)*sizeof(int));
		A->Aj = malloc(2 * nnz*sizeof(int));
		A->Ax = malloc(2 * nnz*sizeof(double));
    }
    MPI_Bcast(&A->nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(A->Ap, n+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(A->Aj, 2*nnz, MPI_INT, 0, MPI_COMM_WORLD); //bug la
    MPI_Bcast(A->Ax, 2*nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
	// /* Allocate memory */
	int quotient = n / P;
	int reste = n % P;
	int n_part = quotient; //nombre d'elements par bloc d'un vecteur de taille n
								  //bloc = partie du vecteur calcul??e lors d'un calcul d'un processeur
    if(my_rank == P-1){
        n_part = quotient + reste;
    }    
	int i_ini = my_rank * quotient; //indice duquel on part pour calculer une partie du vecteur solution x
	// printf("n_part = %d/%d = %d\n", n, P, n_part);

	int recvcounts[P]; //taille du petit tableau de chaque processeur, dans l'ordre croissant de my_rank
    // printf("recvcounts, p = %d :\n", P);
    for(int i = 0;i < P-1;i++){ 
        recvcounts[i] = quotient;
        // printf("%d ", recvcounts[i]);
    }
    recvcounts[P-1] = quotient + reste;   
    // printf("%d\n", recvcounts[P-1]);

    int displs[P]; 
    displs[0] = 0;
    for(int i = 1;i < P;i++){
        displs[i] = i * quotient; 
    }
    // printf("displs : \n");
    // for(int i = 0; i < P; i ++){
    //     printf("%d ", displs[i]);
    // }
    // printf("\n");

	double *mem = malloc(7 * n * sizeof(double));
	if(mem == NULL)
		err(1, "cannot allocate dense vectors");
	double *x = mem;	/* solution vector */
	double *b = mem + n;	/* right-hand side */
	double *scratch = mem + 2 * n;	/* workspace for cg_solve() */

	// /* Prepare right-hand size */
	if (rhs_filename) {	/* load from file */
		FILE *f_b = fopen(rhs_filename, "r");
		if (f_b == NULL)
			err(1, "cannot open %s", rhs_filename);
		fprintf(stderr, "[IO] Loading b from %s\n", rhs_filename);
		for (int i = 0; i < n; i++) {
			if (1 != fscanf(f_b, "%lg\n", &b[i]))
				errx(1, "parse error entry %d\n", i);
		}
		fclose(f_b);
	} else {
		for (int i = 0; i < n; i++)
			b[i] = PRF(i, seed);
	}

	/* solve Ax == b with MPI, witn p processors*/
	double *r = scratch;	        // residue
	double *z = scratch + n;	// preconditioned-residue
	double *p = scratch + 2*n;	// search direction
	double *q = scratch + 3 * n;	// q == Ap
	double *d = scratch + 4 * n;	// diagonal entries of A (Jacobi preconditioning)
	double *q_part = malloc(n_part*sizeof(double)); /* une partie ou bloc du vecteur x */
	double rz_part;
	double pq_part;
	double start = wtime();
	double last_display = start;
	int iter = 0;
	double alpha = 0.0;
	double beta = 0.0;
	double rz = 0.0;
	double pq = 0.0;
	int nz = A->nz;

	// double *d_part;
	// extract_diagonal_part(A, d_part, n_part, i_ini);
	// MPI_Allgather(d_part, n_part, MPI_DOUBLE, d, n_part, MPI_DOUBLE, MPI_COMM_WORLD); /* q <-- A.p */
	
	fprintf(stderr, "[CG] Starting iterative solver\n");
	fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (12.0 * nz + 52.0 * n));
	fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n", 2. * nz, 12. * n);

	extract_diagonal(A, d);

	/* Initialisation des vecteurs */
	for (int i = 0; i < n; i++)
		x[i] = 0.0;
	for (int i = 0; i < n; i++)	// r <-- b - Ax == b
		r[i] = b[i];
	for (int i = 0; i < n; i++)	// z <-- M^(-1).r
		z[i] = r[i] / d[i];
	for (int i = 0; i < n; i++)	// p <-- z
		p[i] = z[i];
	printf("p[n] = %lf\n", p[n-1]);
	printf("z[n] = %lf\n", z[n-1]);
	printf("r[n] = %lf\n", r[n-1]);
	printf("x[n] = %lf\n", x[n-1]);

	/*Algorithme du gradient conjugu?? */
	rz_part = dot_part(r, z, i_ini, n_part);
	MPI_Allreduce(&rz_part, &rz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// rz = dot(r,z)	
	while (norm(n, r) > THRESHOLD){ 
		/* loop invariant : rz = dot(r, z) */
		double old_rz = rz;

	    sp_gemv_part(A, p, q_part, n_part, i_ini);
        MPI_Allgatherv(q_part, n_part, MPI_DOUBLE, q, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD); /* q <-- A.p */    

		pq_part = dot_part(p, q, i_ini, n_part);
		MPI_Allreduce(&pq_part, &pq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// rz = dot(r,z)	
		
		alpha = old_rz / pq;		
		for (int i = 0; i < n; i++)	// x <-- x + alpha*p
			x[i] += alpha * p[i];
		for (int i = 0; i < n; i++)	// r <-- r - alpha*q
			r[i] -= alpha * q[i]; //A*p
		for (int i = 0; i < n; i++)	// z <-- M^(-1).r
			z[i] = r[i] / d[i];
		
		rz_part = dot_part(r, z, i_ini, n_part);
		MPI_Allreduce(&rz_part, &rz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// rz = dot(r,z)	

		beta = rz / old_rz;
		for (int i =0; i < n; i++)	// p <-- z + beta*p
			p[i] = z[i] + beta * p[i];
		iter++;
		double t = wtime();
		if (t - last_display > 0.5) {
			double rate = iter / (t - start);	// iterations per s.
			double GFLOPs = 1e-9 * rate * (2 * nz + 12 * n);
			fprintf(stderr, "\r     ---> error : %2.2e, iter : %d (%.1f it/s, %.2f GFLOPs)", norm(n, r), iter, rate, GFLOPs);
			fflush(stdout);
			last_display = t;
		}
	}    
	fprintf(stderr, "\n     ---> Finished in %.1fs and %d iterations\n", wtime() - start, iter);



	/* Check result */
	if (safety_check) {
		double *y = scratch;
		sp_gemv(A, x, y);	// y = Ax
		for (int i = 0; i < n; i++){	// y = Ax - b
			y[i] -= b[i];
			// printf("y[%d] = %lf ", i, y[i]);
		}
		fprintf(stderr, "[check] max error = %2.2e\n", norm(n, y));
	}

	/* Dump the solution vector */

	if(my_rank == 0){
		FILE *f_x = stdout;
		// printf("solution filename = %s\n", solution_filename);
		if (solution_filename != NULL) {
			f_x = fopen(solution_filename, "w");
			if (f_x == NULL)
				err(1, "cannot open solution file %s", solution_filename);
			fprintf(stderr, "[IO] writing solution to %s\n", solution_filename);
		}
		// for (int i = 0; i < n; i++)
		// 	fprintf(f_x, "%a\n", x[i]);
	}

	free(mem);
	free(q_part);
	free(nnz2);
	free(A);
	// return EXIT_SUCCESS;

	MPI_Finalize();
}
