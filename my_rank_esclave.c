// scratch ddéjà défini
int i_ini = my_rank*n_part;
double *r = scratch;	        // residue
double *z = scratch + n;	// preconditioned-residue
double *p = scratch + 2*n;	// search direction
double *q = scratch + 3 * n;	// q == Ap
double *d = scratch + 4 * n;	// diagonal entries of A (Jacobi preconditioning)
double *q_part = malloc(n_part*sizeof(double)); /* une partie ou bloc du vecteur x */
double rz_part;
double pq_part;

extract_diagonal_part(A, d_part, n_part, i_ini);
MPI_Allgather(d_part, n_part, MPI_DOUBLE, d, n_part, MPI_DOUBLE, MPI_COMM_WORLD); /* q <-- A.p */

double start = wtime();
double last_display = start;
int iter = 0;
double alpha = 0.0;
double beta = 0.0;
double rz = 0.0;
double pq = 0.0;
/* Initialisation des vecteurs */
for (int i = 0; i < n_part; i++)
	x[i] = 0.0;
for (int i = 0; i < n; i++)	// r <-- b - Ax == b
	r[i] = b[i];
for (int i = 0; i < n; i++)	// z <-- M^(-1).r
	z[i] = r[i] / d[i];
for (int i = 0; i < n; i++)	// p <-- z
	p[i] = z[i];

/*Algorithme du gradient conjugué */
// int recvcount = n_part*nbProc;
double start = wtime();
double last_display = start;
int iter = 0;
MPI_Allreduce(MPI_IN_PLACE, &rz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// rz = dot(r,z)	
while (norm(n, r) > THRESHOLD){ ///////PAS SUR SUR QUELLE CONDITION METTRE
	/* loop invariant : rz = dot(r, z) */
	double old_rz = rz;
    MPI_Allgather(q_part, n_part, MPI_DOUBLE, q, n_part, MPI_DOUBLE, MPI_COMM_WORLD); /* q <-- A.p */
	MPI_Allreduce(&pq_part, &pq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// rz = dot(r,z)	
	alpha = old_rz / pq;		
	for (int i = 0; i < n_part; i++)	// x <-- x + alpha*p
		x[i] += alpha * p[i + i_ini];
	for (int i =0; i < n; i++)	// r <-- r - alpha*q
		r[i] -= alpha * q[i]; //A*p
	for (int i = 0; i < n; i++)	// z <-- M^(-1).r
		z[i] = r[i] / d[i];
	MPI_Allreduce(&rz_part, &rz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);// rz = dot(r,z)	
	beta = rz / old_rz;
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

free(scratch);
free(q_part);
