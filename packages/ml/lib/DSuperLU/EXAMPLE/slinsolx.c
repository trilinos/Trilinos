

/*
 * -- SuperLU routine (version 1.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include "ssp_defs.h"
#include "util.h"

main(int argc, char *argv[])
{
    char           fact[1], equed[1], trans[1], refact[1];
    SuperMatrix  A, L, U;
    SuperMatrix  B, X;
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
    float         *a;
    int            *asub, *xa;
    int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    factor_param_t iparam;
    int            info, lwork, nrhs, ldb, ldx, panel_size, relax;
    int            m, n, nnz, permc_spec;
    float         *rhsb, *rhsx, *xact;
    float         *R, *C;
    float         *ferr, *berr;
    float         u, rpg, rcond;
    int            i, j, rowequ, colequ, firstfact;
    mem_usage_t    mem_usage;
    void    parse_command_line();

    /* Defaults */
    lwork = 0;
    *fact      = 'N';
    *equed     = 'N';
    *trans     = 'N';
    *refact    = 'N';
    nrhs       = 1;
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    u          = 1.0;
    parse_command_line(argc, argv, &lwork, &panel_size, &relax, &u,
		       fact, trans, refact);
    firstfact = lsame_(fact, "F") || lsame_(refact, "Y");

    iparam.panel_size        = panel_size;
    iparam.relax             = relax;
    iparam.diag_pivot_thresh = u;
    iparam.drop_tol          = -1;
    
    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    ABORT("SLINSOLX: cannot allocate work[]");
	}
    }

    
    sreadhb(&m, &n, &nnz, &a, &asub, &xa);
    
    sCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, _S, GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    if ( !(rhsb = floatMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsx = floatMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
    sCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, DN, _S, GE);
    sCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, DN, _S, GE);
    xact = floatMalloc(n * nrhs);
    ldx = n;
    sGenXtrue(n, nrhs, xact, ldx);
    sFillRHS(trans, nrhs, xact, ldx, &A, &B);
    
    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     */    	
    permc_spec = 0;
    get_perm_c(permc_spec, &A, perm_c);

    if ( !(R = (float *) SUPERLU_MALLOC(A.nrow * sizeof(float))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (float *) SUPERLU_MALLOC(A.ncol * sizeof(float))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (float *) SUPERLU_MALLOC(nrhs * sizeof(float))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (float *) SUPERLU_MALLOC(nrhs * sizeof(float))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");

    
    /* Solve the system and compute the condition number
       and error bounds using dgssvx.      */
    
    sgssvx(fact, trans, refact, &A, &iparam, perm_c, perm_r, etree,
	   equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond,
	   ferr, berr, &mem_usage, &info);

    printf("sgssvx(): info %d\n", info);

    if ( info == 0 || info == n+1 ) {

	printf("Recip. pivot growth = %e\n", rpg);
	printf("Recip. condition number = %e\n", rcond);
	printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
	for (i = 0; i < nrhs; ++i) {
	    printf("%8d%16e%16e\n", i+1, ferr[i], berr[i]);
	}
	       
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       mem_usage.expansions);
	     
	fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork >= 0 ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    }
}

/*  
 * Parse command line options to get relaxed snode size, panel size, etc.
 */
void
parse_command_line(int argc, char *argv[], int *lwork, int *w, int *relax, 
                   float *u, char *fact, char *trans, char *refact)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "hl:w:r:u:f:t:p:e:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-l <int> - length of work[*] array\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-r <int> - max size of relaxed supernodes\n");
	    printf("\t-u <int> - pivoting threshold\n");
	    printf("\t-f <'N', 'F' or 'E'> - fact (see linsolvex.c)\n");
	    printf("\t-t <'N' or 'T'> - solve transposed system or not\n");
	    printf("\t-p <'N' or 'Y'> - refactor or not\n");
	    exit(1);
	    break;
	  case 'l': *lwork = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'r': *relax = atoi(optarg); 
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
	  case 'f': *fact = *optarg;
	            break;
	  case 't': *trans = *optarg;
	            break;
	  case 'p': *refact = *optarg;
	            break;
  	}
    }
}
