/* 
 * Sequential implementation of the Conjugate Gradient Method.
 *
 * Authors : Lilia Ziane Khodja & Charles Bouillaguet
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

// void maitre_esclave_root2(int* x){
// /* Premier tour des ordres envoyes aux esclaves */
// 	for(int i = 1;i < nbProc;i++){
// 		dest = i;
// 		MPI_Send(&i_block, 1, MPI_INT, dest, INDICE, MPI_COMM_WORLD);
// 		i_block += 1;			
// 	}
// }

void maitre_esclave_root_produit_scalaire(double* x, double* a, double* b, double* x_part, int tagMission, int nbProc, int n, int n_part, int nbOfBlock){
	//tagMission d??crit la mission en cours qui est calcul du produit scalaire rz
	//, du produit scalaire pq ou produit matriciel A*p = q 
/* Premier tour des ordres envoyes aux esclaves */
	int i_ini = 0; //indice duquel on part pour calculer une partie du vecteur solution x
	int i_block = 0; //numero du bloc courant en train d'etre calcul??
	MPI_Status status;
	int dest;
	int bTmp = 0; //numero du dernier bloc du vecteur x qui vient d'??tre calcul??
	int idTmp;
	enum tagType {INDICE, TRAITEMENT, STOP, DOT_RZ, DOT_PQ, MATPROD};

	for(int i = 1;i < nbProc;i++){
		dest = i;
		MPI_Send(&i_block, 1, MPI_INT, dest, tagMission, MPI_COMM_WORLD);
		i_block += 1;					
		MPI_Send(&a, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
		MPI_Send(&b, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);

		
		// if(tagMission == DOT_RZ){
		// 	MPI_Send(&r, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
		// 	MPI_Send(&z, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);
		// }
		// else if(tagMission == DOT_PQ){
		// 	MPI_Send(&p, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
		// 	MPI_Send(&q, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);
		// }
	}

	/* Envoi des ordres */
	while(i_block != nbOfBlock){
		//le maitre recoit le num??ro du dernier bloc calcul??
		MPI_Recv(&bTmp, 1, MPI_INT, MPI_ANY_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);	
		//le maitre recoit le dernier bloc calcul??
		MPI_Recv(x_part, sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		//Le ma??tre envoie la mission suivante ?? l'esclave
		dest = status.MPI_SOURCE;
		MPI_Send(&i_block, 1, MPI_INT, dest, tagMission, MPI_COMM_WORLD);
		i_block += 1;
		// if(tagMission == DOT_RZ){
		// 	MPI_Send(&r, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
		// 	MPI_Send(&z, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);
		// }
		// else if(tagMission == DOT_PQ){
		// 	MPI_Send(&p, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
		// 	MPI_Send(&q, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);
		// }
		
		/* Le maitre recopie le contenu de la partie du vecteur qu'il a re??u */
		*x += *x_part;
	}

	/* Reception des derniers travaux des esclaves */
	for(int i = 1;i < nbProc;i++){
		MPI_Recv(&bTmp, 1, MPI_INT, MPI_ANY_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);	
		MPI_Recv(x_part, sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		//On recup??re le num??ro du processus (idTmp) qui vient d'envoyer son travail au ma??tre
		idTmp = status.MPI_SOURCE;
		//On recopie le travail de l'esclave
		*x += *x_part;
		//On dit ?? l'esclave de ne plus travailler
		MPI_Send(&idTmp, 1, MPI_INT, idTmp, STOP, MPI_COMM_WORLD);			
	}
}

void maitre_esclave_root_produit_matriciel(double* x, double* a, double* x_part, int tagMission, int nbProc, int n, int n_part, int nbOfBlock){
	//tagMission d??crit la mission en cours qui est calcul du produit scalaire rz
	//, du produit scalaire pq ou produit matriciel A*p = q 
/* Premier tour des ordres envoyes aux esclaves */
	int i_ini = 0; //indice duquel on part pour calculer une partie du vecteur solution x
	int i_block = 0; //numero du bloc courant en train d'etre calcul??
	MPI_Status status;
	int dest;
	int bTmp = 0; //numero du dernier bloc du vecteur x qui vient d'??tre calcul??
	int idTmp;
	enum tagType {INDICE, TRAITEMENT, STOP, DOT_RZ, DOT_PQ, MATPROD};

	MPI_Send(a, n*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);
	for(int i = 1;i < nbProc;i++){
		dest = i;
		MPI_Send(&i_block, 1, MPI_INT, dest, tagMission, MPI_COMM_WORLD);
		i_block += 1;			
	}

	/* Envoi des ordres */
	while(i_block != nbOfBlock){
		//le maitre recoit le num??ro du dernier bloc calcul??
		MPI_Recv(&bTmp, 1, MPI_INT, MPI_ANY_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);	
		//le maitre recoit le dernier bloc calcul??
		MPI_Recv(x_part, n_part*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		dest = status.MPI_SOURCE;
		MPI_Send(&i_block, 1, MPI_INT, dest, tagMission, MPI_COMM_WORLD);
		i_block += 1;
		
		/* Le maitre recopie le contenu de la partie du vecteur qu'il a re??u */
		for(int i = 0;i<n_part;i++){
			x[bTmp*n_part + i] = x_part[i];
		}
	}

	/* Reception des derniers travaux des esclaves */
	for(int i = 1;i < nbProc;i++){
		MPI_Recv(&bTmp, 1, MPI_INT, MPI_ANY_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);	
		MPI_Recv(x_part, n_part*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		//On recup??re le num??ro du processus (idTmp) qui vient d'envoyer son travail au ma??tre
		idTmp = status.MPI_SOURCE;
		//On recopie le travail de l'esclave
		for(int j = 0;j<n_part;j++){
			(x + bTmp*n_part)[j] = x_part[j];
		}	
		//On dit ?? l'esclave de ne plus travailler
		MPI_Send(&idTmp, 1, MPI_INT, idTmp, STOP, MPI_COMM_WORLD);			
	}
}

/*************************** Matrix accessors *********************************/

/* Copy the diagonal of A into the vector d. */
void extract_diagonal(const struct csr_matrix_t *A, double *d, int n_part, int i_ini)
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
void sp_gemv(const struct csr_matrix_t *A, const double *x, double *y, int n)
{
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = 0; i < n; i++) {
		y[i] = 0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++) {
			int j = Aj[u]; //j = Aj[Ap[i]]
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
		y[i] = 0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++) {
			int j = Aj[u]; //j = Aj[Ap[i]]
			double A_ij = Ax[u];
			y[i] += A_ij * x[j];
		}
	}
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

/*********************** conjugate gradient algorithm *************************/

/* Solve Ax == b (the solution is written in x). Scratch must be preallocated of size 6n */
// void cg_solve(const struct csr_matrix_t *A, const double *b, double *x, const double epsilon, double *scratch, int n_part, int n)
// {
// 	int nz = A->nz;

// 	fprintf(stderr, "[CG] Starting iterative solver\n");
// 	fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (12.0 * nz + 52.0 * n_part));
// 	fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n", 2. * nz, 12. * n_part);

// 	//	double *mem = malloc(7 * n * sizeof(double));
// 	// double *x = mem;	/* solution vector */
// 	// double *b = mem + n;	/* right-hand side */
// 	//	double *scratch = mem + 2 * n;	/* workspace for cg_solve() */
// 	double *r = scratch;	        // residue
// 	double *z = scratch + n;	// preconditioned-residue
// 	double *p = scratch + 2*n;	// search direction
// 	double *q = scratch + 3 * n;	// q == Ap
// 	double *d = scratch + 4 * n;	// diagonal entries of A (Jacobi preconditioning)

// 	/* Isolate diagonal */
// 	extract_diagonal(A, d, n_part, i_ini);

// 	/* 
// 	 * This function follows closely the pseudo-code given in the (english)
// 	 * Wikipedia page "Conjugate gradient method". This is the version with
// 	 * preconditionning.
// 	 */

// 	/* We use x == 0 --- this avoids the first matrix-vector product. */
// 	for (int i = 0; i < n_part; i++)
// 		x[i] = 0.0;
// 	for (int i = 0; i < n; i++)	// r <-- b - Ax == b
// 		r[i] = b[i];
// 	for (int i = 0; i < n; i++)	// z <-- M^(-1).r
// 		z[i] = r[i] / d[i];
// 	for (int i = 0; i < n; i++)	// p <-- z
// 		p[i] = z[i];

// 	double rz = dot_part(r, z, i_ini, n_part);
// 	double start = wtime();
// 	double last_display = start;
// 	int iter = 0;
// 	while (norm_part(r,i_ini,n_part) > epsilon){ ///////PAS SUR SUR QUELLE CONDITION METTRE
// 		/* loop invariant : rz = dot(r, z) */
// 		double old_rz = rz;
// 		sp_gemv(A, p, q, n);	/* q <-- A.p */
// 		double alpha = old_rz / dot_part(p, q, i_ini, n_part);
// 		for (int i = 0; i < n_part; i++)	// x <-- x + alpha*p
// 			x[i] += alpha * p[i + i_ini];
// 		for (int i =0; i < n; i++)	// r <-- r - alpha*q
// 			r[i] -= alpha * q[i]; //A*p
// 		for (int i = 0; i < n; i++)	// z <-- M^(-1).r
// 			z[i] = r[i] / d[i];
// 		rz = dot_part(r, z, i_ini, n_part);	// restore invariant : rz = dot(r, z)
// 		double beta = rz / old_rz;
// 		for (int i =0; i < n; i++)	// p <-- z + beta*p
// 			p[i] = z[i] + beta * p[i];
// 		iter++;
// 		double t = wtime();
// 		if (t - last_display > 0.5) {
// 			/* verbosity */
// 			double rate = iter / (t - start);	// iterations per s.
// 			double GFLOPs = 1e-9 * rate * (2 * nz + 12 * n);
// 			fprintf(stderr, "\r     ---> error : %2.2e, iter : %d (%.1f it/s, %.2f GFLOPs)", norm(n, r), iter, rate, GFLOPs);
// 			fflush(stdout);
// 			last_display = t;
// 		}
// 	}
// 	fprintf(stderr, "\n     ---> Finished in %.1fs and %d iterations\n", wtime() - start, iter);
// }

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
	int nbProc; //number of process
	int i_ini = 0; //indice duquel on part pour calculer une partie du vecteur solution x
	int i_block = 0; //numero du bloc courant en train d'etre calcul??
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
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

	printf("hello i am process %s number %d\n", processor_name, my_rank);

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
    // printf("%d %d %lf %d %d %lf\n", A->Ap[0], A->Aj[0], A->Ax[0], A->Ap[50], A->Aj[50], A->Ax[50]);
  
	// /* Allocate memory */
	// n = taille du vecteur x
	//n(cfd1) = 70 656 = n
	int n_part = 92;//nombre d'elements par bloc du vecteur x
								  //bloc = partie du vecteur calcul??e lors d'un calcul d'un processeur
	int nbOfBlock = n/n_part;//nombre de blocs du vecteur x = nb de calculs a effectuer
	int reste = n % n_part;//indique si on rajoute un bloc si besoin pour calculer 
			//le reste du vecteur si la division a un reste
	double *x_part = malloc(n_part*sizeof(double)); /* une partie ou bloc du vecteur x */
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

	/* solve Ax == b with MPI, witn nbProc processors*/
	double *r = scratch;	        // residue
	double *z = scratch + n;	// preconditioned-residue
	double *p = scratch + 2*n;	// search direction
	double *q = scratch + 3 * n;	// q == Ap
	double *d = scratch + 4 * n;	// diagonal entries of A (Jacobi preconditioning)
	double *q_part = malloc(n_part*sizeof(double));
	double *rz_part = calloc(1,sizeof(double));
	double *pq_part = calloc(1,sizeof(double));

	if(my_rank == 0){		
		double start = wtime();
		double last_display = start;
		double alpha = 0.0;
		double beta = 0.0;
		double *rz = calloc(1,sizeof(double));
		double *pq = calloc(1,sizeof(double));

	// /* Initialisation des vecteurs */
		for (int i = 0; i < n_part; i++)
			x[i] = 0.0;
		for (int i = 0; i < n; i++)	// r <-- b - Ax == b
			r[i] = b[i];
		for (int i = 0; i < n; i++)	// z <-- M^(-1).r
			z[i] = r[i] / d[i];
		for (int i = 0; i < n; i++)	// p <-- z
			p[i] = z[i];

	// /*Algorithme du gradient conjugu?? */
		int nz = A->nz;
		int recvcount = n_part*nbProc;
		int iter = 0;	

		fprintf(stderr, "[CG] Starting iterative solver\n");
		fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (12.0 * nz + 52.0 * n));
		fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n", 2. * nz, 12. * n);

		maitre_esclave_root_produit_scalaire(rz, r, z, rz_part, DOT_RZ, nbProc, n, n_part, nbOfBlock); // rz = dot(r,z)
		printf("rz = %lf\n", *rz);
		while (norm(n, r) > THRESHOLD){
		/* loop invariant : rz = dot(r, z) */
			double old_rz = *rz;
			maitre_esclave_root_produit_matriciel(q, p, q_part, MATPROD, nbProc, n,  n_part, nbOfBlock); /* q <-- A.p */
			maitre_esclave_root_produit_scalaire(pq, p, q, pq_part, DOT_PQ, nbProc, n,  n_part, nbOfBlock); // pq <-- dot(p, q)
			alpha = old_rz / *pq;		
			for (int i = 0; i < n_part; i++)	// x <-- x + alpha*p
				x[i] += alpha * p[i + i_ini];
			for (int i =0; i < n; i++)	// r <-- r - alpha*q
				r[i] -= alpha * q[i]; //A*p
			for (int i = 0; i < n; i++)	// z <-- M^(-1).r
				z[i] = r[i] / d[i];
			maitre_esclave_root_produit_scalaire(rz, r, z, rz_part, DOT_RZ, nbProc, n,  n_part, nbOfBlock); // rz = dot(r,z)
			beta = *rz / old_rz;
			for (int i =0; i < n; i++)	// p <-- z + beta*p
				p[i] = z[i] + beta * p[i];
			iter++;
			double t = wtime();
			if (t - last_display > 0.5) {
				/* verbosity */
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
			//sp_gemv(const struct csr_matrix_t *A, const double *x, double *y, int n, int i_ini)
			//sp_gemv(A, p, q, n, i_ini);
			sp_gemv(A, x, y, n);	// y = Ax
			for (int i = 0; i < n; i++)	// y = Ax - b
				y[i] -= b[i];
			fprintf(stderr, "[check] max error = %2.2e\n", norm(n, y));
		}

		/* Dump the solution vector */
		FILE *f_x = stdout;
		if (solution_filename != NULL) {
			f_x = fopen(solution_filename, "w");
			if (f_x == NULL)
				err(1, "cannot open solution file %s", solution_filename);
			fprintf(stderr, "[IO] writing solution to %s\n", solution_filename);
		}
		for (int i = 0; i < n; i++)
			fprintf(f_x, "%a\n", x[i]);
		return EXIT_SUCCESS;
	 }
	else{// si le processus n'est pas le ma??tre mais un esclave
		MPI_Recv(&i_block, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if(status.MPI_TAG == DOT_RZ){
				/* Calcul d'une partie du produit scalaire (pour une partie des composantes) */
				MPI_Recv(r, n*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
				MPI_Recv(z, n*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		}
		else if(status.MPI_TAG == DOT_PQ){
				/* Calcul d'une partie du produit scalaire (pour une partie des composantes) */
				MPI_Recv(p, n*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
				MPI_Recv(q, n*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		}
		else if(status.MPI_TAG == MATPROD){
				MPI_Recv(p, n*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);		
		}
		while(1){
			MPI_Recv(&i_block, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			i_ini = i_block * n_part;

			if(status.MPI_TAG == DOT_RZ){
				*rz_part = dot_part(r, z, i_ini, n_part);
				MPI_Send(&rz_part, 1, MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
			}
			else if(status.MPI_TAG == DOT_PQ){
				*pq_part = dot_part(p, q, i_ini, n_part);
				MPI_Send(&pq_part, 1, MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
			}
			else if(status.MPI_TAG == MATPROD){
				sp_gemv_part(A, p, q_part, n_part, i_ini);
				MPI_Send(q_part, n_part*sizeof(double), MPI_DOUBLE, dest, TRAITEMENT, MPI_COMM_WORLD);	
			}
			else if(status.MPI_TAG == STOP){
				break;
			}

			/* Envoi le num??ro du bloc calcul?? */
			bTmp = i_block; //indice temporaire du dernier bloc trait??
			dest = 0;
			MPI_Send(&bTmp, 1, MPI_INT, dest, TRAITEMENT, MPI_COMM_WORLD);
		}
	}

	/* Affichage de la sortie */
	fin = my_gettimeofday();
	fprintf(stderr, "Temps total de calcul du processeur %d : %g sec\n", my_rank, fin - debut);
	//printf("my_rank = %d, fini\n",my_rank);

	MPI_Finalize();
}
