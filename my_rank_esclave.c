// scratch ddéjà défini
double *r = scratch;	        // residue
double *z = scratch + n;	// preconditioned-residue
double *p = scratch + 2*n;	// search direction
double *q = scratch + 3 * n;	// q == Ap
double *d = scratch + 4 * n;	// diagonal entries of A (Jacobi preconditioning)
double *q_part = malloc(n_part*sizeof(double)); /* une partie ou bloc du vecteur x */
double rz = 0.0;
extract_diagonal(A, d, n_part, i_ini);

if(my_rank != 0){
	MPI_Recv(&i_block, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	i_ini = i_block * n_part;

	if(status.MPI_TAG == DOT_RZ){
		/* Calcul d'une partie du produit scalaire (pour une partie des composantes) */
		MPI_Recv(r, n_part*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		MPI_Recv(z, n_part*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		rz = dot_part(r, z, i_ini, n_part);
		/* Envoi le résultat du calcul au maître et le numéro du bloc calculé */
		dest = 0;
		//Il y aura un reduce.
	}
	else if(status.MPI_TAG == DOT_PQ){
		/* Calcul d'une partie du produit scalaire (pour une partie des composantes) */
		MPI_Recv(p, n_part*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		MPI_Recv(q, n_part*sizeof(double), MPI_DOUBLE, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
		pq = dot_part(p, q, i_ini, n_part);
		//Il y aura un reduce de la racine.
	}
	else if(status.MPI_TAG == MATPROD){
		void sp_gemv_part(A, p, q_part, n_part, i_ini);
		//Il y aura un reduce de la racine.
	}

	/* Envoi le numéro du bloc calculé */
	bTmp = i_block; //indice temporaire du dernier bloc traité
	dest = 0;
	MPI_Send(&bTmp, 1, MPI_INT, dest, TRAITEMENT, MPI_COMM_WORLD);
}



	double rz = dot_part(r, z, i_ini, n_part);
	double start = wtime();
	double last_display = start;
	int iter = 0;
	while (norm_part(r,i_ini,n_part) > THRESHOLD){ ///////PAS SUR SUR QUELLE CONDITION METTRE
		/* loop invariant : rz = dot(r, z) */
		double old_rz = rz;
		sp_gemv(A, p, q, n);	/* q <-- A.p */
		double alpha = old_rz / dot_part(p, q, i_ini, n_part);
		for (int i = 0; i < n_part; i++)	// x <-- x + alpha*p
			x[i] += alpha * p[i + i_ini];
		for (int i =0; i < n; i++)	// r <-- r - alpha*q
			r[i] -= alpha * q[i]; //A*p
		for (int i = 0; i < n; i++)	// z <-- M^(-1).r
			z[i] = r[i] / d[i];
		rz = dot_part(r, z, i_ini, n_part);	// restore invariant : rz = dot(r, z)
		double beta = rz / old_rz;
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



if(my_rank == 0){
	double start = wtime();
	double last_display = start;
	int iter = 0;

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
	int recvcount = n_part*nbProc;
	MPI_Reduce(MPI_IN_PLACE, &rz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);// rz = dot(r,z)	
    MPI_Gather(q_part, n_part, MPI_INT, q, recvcount, MPI_INT, root, MPI_COMM_WORLD);



}