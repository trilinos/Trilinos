

/*
 * -- SuperLU routine (version 1.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/*
 * File name:		cdrive.c
 * Purpose:             MAIN test program
 */
#include <string.h>
#include "csp_defs.h"
#include "util.h"

#define NTESTS    5      /* Number of test types */
#define NTYPES    11     /* Number of matrix types */
#define NTRAN     2    
#define THRESH    20.0
#define FMT1      "%10s:n=%d, test(%d)=%12.5g\n"
#define	FMT2      "%10s:fact=%c, trans=%c, refact=%c, equed=%c, n=%d, imat=%d, test(%d)=%12.5g\n"
#define FMT3      "%10s:info=%d, izero=%d, n=%d, nrhs=%d, imat=%d, nfail=%d\n"


main(int argc, char *argv[])
{
/* 
 * Purpose
 * =======
 *
 * CDRIVE is the main test program for the COMPLEX linear 
 * equation driver routines CGSSV and CGSSVX.
 * 
 * The program is invoked by a shell script file -- ctest.csh.
 * The output from the tests are written into a file -- ctest.out.
 *
 * =====================================================================
 */
    complex         *a, *a_save;
    int            *asub, *asub_save;
    int            *xa, *xa_save;
    SuperMatrix  A, B, X, L, U;
    SuperMatrix  ASAV, AC;
    NCformat       *Astore, *Ustore;
    SCformat       *Lstore;
    factor_param_t iparam;
    mem_usage_t    mem_usage;
    int            *perm_r; /* row permutation from partial pivoting */
    int            *perm_c, *pc_save; /* column permutation */
    int            *etree;
    complex  zero = {0.0, 0.0};
    float         thresh, drop_tol;
    float         *R, *C;
    float         *ferr, *berr;
    float         *rwork;
    complex	   *wwork;
    void           *work;
    int            info, lwork, nrhs, panel_size, relax;
    int            m, n, nnz;
    complex         *xact;
    complex         *rhsb, *solx, *bsav;
    int            ldb, ldx;
    float         rpg, rcond;
    int            i, j, k1, rowequ, colequ;
    float         rowcnd, colcnd, amax;
    int            maxsuper, rowblk, colblk;
    int            prefact, nofact, equil, iequed, notrans, norefact;
    int            nt, nrun, nfail, nerrs, imat, fimat, nimat;
    int            nfact, ifact, nrefact, irefact, itran;
    int            kl, ku, mode, lda;
    int            zerot, izero, ioff;
    float         anorm, cndnum;
    complex         *Afull;
    float         result[NTESTS];
    double         *utime;
    flops_t        *ops;
    extern SuperLUStat_t SuperLUStat;
    void    parse_command_line();
    static char    matrix_type[8];
    static char    fact[1], trans[1], equed[1], refact[1],
                   path[3], sym[1], dist[1];

    /* Fixed set of parameters */
    int            iseed[]  = {1988, 1989, 1990, 1991};
    static char    equeds[]  = {'N', 'R', 'C', 'B'};
    static char    facts[]   = {'F', 'N', 'E'};
    static char    refacts[] = {'Y', 'N'};
    static char    transs[]  = {'N', 'T'};

    /* Some function prototypes */
    extern int sp_cget02(char *, int, int, int, SuperMatrix *, complex *,
                         int, complex *, int, float *resid);
    extern int sp_cget04(int, int, complex *, int, 
                         complex *, int, float rcond, float *resid);
    extern int sp_cget07(char *, int, int, SuperMatrix *, complex *, int,
                         complex *, int, complex *, int, 
                         float *, float *, float *);
    extern int clatb4_(char *, int *, int *, int *, char *, int *, int *, 
	               float *, int *, float *, char *);
    extern int clatms_(int *, int *, char *, int *, char *, float *d,
                       int *, float *, float *, int *, int *,
                       char *, complex *, int *, complex *, int *);
    extern int sp_cconvert(int, int, complex *, int, int, int,
	                   complex *a, int *, int *, int *);


    /* Executable statements */

    strcpy(path, "CGE");
    nrun  = 0;
    nfail = 0;
    nerrs = 0;


    /* Defaults */
    lwork      = 0;
    n          = 1;
    nrhs       = 1;
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    strcpy(matrix_type, "LA");
    parse_command_line(argc, argv, matrix_type, &n,
		       &panel_size, &relax, &nrhs, &maxsuper,
		       &rowblk, &colblk, &lwork);
    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    fprintf(stderr, "expert: cannot allocate %d bytes\n", lwork);
	    exit (-1);
	}
    }

    iparam.panel_size = panel_size;
    iparam.relax      = relax;
    iparam.diag_pivot_thresh = 1.0;
    iparam.drop_tol   = 0.0;
    
    if ( strcmp(matrix_type, "LA") == 0 ) {
	/* Test LAPACK matrix suite. */
	m = n;
	lda = MAX(n, 1);
	nnz = n * n;        /* upper bound */
	fimat = 1;
	nimat = NTYPES;
	Afull = complexCalloc(lda * n);
	callocateA(n, nnz, &a, &asub, &xa);
    } else {
	/* Read a sparse matrix */
	fimat = nimat = 0;
	creadhb(&m, &n, &nnz, &a, &asub, &xa);
    }

    callocateA(n, nnz, &a_save, &asub_save, &xa_save);
    rhsb = complexMalloc(m * nrhs);
    bsav = complexMalloc(m * nrhs);
    solx = complexMalloc(n * nrhs);
    ldb  = m;
    ldx  = n;
    cCreate_Dense_Matrix(&B, m, nrhs, rhsb, ldb, DN, _C, GE);
    cCreate_Dense_Matrix(&X, n, nrhs, solx, ldx, DN, _C, GE);
    xact = complexMalloc(n * nrhs);
    etree   = intMalloc(n);
    perm_r  = intMalloc(n);
    perm_c  = intMalloc(n);
    pc_save = intMalloc(n);
    R       = (float *) SUPERLU_MALLOC(m*sizeof(float));
    C       = (float *) SUPERLU_MALLOC(n*sizeof(float));
    ferr    = (float *) SUPERLU_MALLOC(nrhs*sizeof(float));
    berr    = (float *) SUPERLU_MALLOC(nrhs*sizeof(float));
    j = MAX(m,n) * MAX(4,nrhs);    
    rwork   = (float *) SUPERLU_MALLOC(j*sizeof(float));
    for (i = 0; i < j; ++i) rwork[i] = 0.;
    if ( !R ) ABORT("SUPERLU_MALLOC fails for R");
    if ( !C ) ABORT("SUPERLU_MALLOC fails for C");
    if ( !ferr ) ABORT("SUPERLU_MALLOC fails for ferr");
    if ( !berr ) ABORT("SUPERLU_MALLOC fails for berr");
    if ( !rwork ) ABORT("SUPERLU_MALLOC fails for rwork");
    wwork   = complexCalloc( MAX(m,n) * MAX(4,nrhs) );

    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i] = i;

    for (imat = fimat; imat <= nimat; ++imat) { /* All matrix types */
	
	if ( imat ) {

	    /* Skip types 5, 6, or 7 if the matrix size is too small. */
	    zerot = (imat >= 5 && imat <= 7);
	    if ( zerot && n < imat-4 )
		continue;
	    
	    /* Set up parameters with CLATB4 and generate a test matrix
	       with CLATMS.  */
	    clatb4_(path, &imat, &n, &n, sym, &kl, &ku, &anorm, &mode,
		    &cndnum, dist);

	    clatms_(&n, &n, dist, iseed, sym, &rwork[0], &mode, &cndnum,
		    &anorm, &kl, &ku, "No packing", Afull, &lda,
		    &wwork[0], &info);

	    if ( info ) {
		printf(FMT3, "CLATMS", info, izero, n, nrhs, imat, nfail);
		continue;
	    }

	    /* For types 5-7, zero one or more columns of the matrix
	       to test that INFO is returned correctly.   */
	    if ( zerot ) {
		if ( imat == 5 ) izero = 1;
		else if ( imat == 6 ) izero = n;
		else izero = n / 2 + 1;
		ioff = (izero - 1) * lda;
		if ( imat < 7 ) {
		    for (i = 0; i < n; ++i) Afull[ioff + i] = zero;
		} else {
		    for (j = 0; j < n - izero + 1; ++j)
			for (i = 0; i < n; ++i)
			    Afull[ioff + i + j*lda] = zero;
		}
	    } else {
		izero = 0;
	    }

	    /* Convert to sparse representation. */
	    sp_cconvert(n, n, Afull, lda, kl, ku, a, asub, xa, &nnz);

	} else {
	    izero = 0;
	    zerot = 0;
	}
	
	cCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, _C, GE);

	/* Save a copy of matrix A in ASAV */
	cCreate_CompCol_Matrix(&ASAV, m, n, nnz, a_save, asub_save, xa_save,
			      NC, _C, GE);
	cCopy_CompCol_Matrix(&A, &ASAV);
	
	/* Form exact solution. */
	cGenXtrue(n, nrhs, xact, ldx);
	
	for (iequed = 0; iequed < 4; ++iequed) {
	    *equed = equeds[iequed];
	    if (iequed == 0) nfact = 3;
	    else nfact = 1;

	    for (ifact = 0; ifact < nfact; ++ifact) {
		*fact = facts[ifact];
		if (ifact == 0) nrefact = 1;
		else nrefact = 2;

		for (irefact = 0; irefact < nrefact; ++irefact) {
		    *refact = refacts[irefact];
		    norefact  = lsame_(refact, "N");
		    prefact   = lsame_(fact, "F") || ( ! norefact );
		    nofact    = lsame_(fact, "N");
		    equil     = lsame_(fact, "E");

		    /* Restore the matrix A. */
		    cCopy_CompCol_Matrix(&ASAV, &A);
			
		    if ( zerot ) {
			if ( prefact ) continue;
		    } else if ( ! nofact ) {
			if ( equil || iequed ) {
			    /* Compute row and column scale factors to
			       equilibrate matrix A.    */
			    cgsequ(&A, R, C, &rowcnd, &colcnd,
				      &amax, &info);

			    /* Force equilibration. */
			    if ( !info && n > 0 ) {
				if ( lsame_(equed, "R") ) {
				    rowcnd = 0.;
				    colcnd = 1.;
				} else if ( lsame_(equed, "C") ) {
				    rowcnd = 1.;
				    colcnd = 0.;
				} else if ( lsame_(equed, "B") ) {
				    rowcnd = 0.;
				    colcnd = 0.;
				}
			    }
			
			    /* Equilibrate the matrix. */
			    claqgs(&A, R, C, rowcnd, colcnd, amax, equed);
			}
		    }
		    
		    if ( prefact ) {	/* First time factor */
			
			StatInit(panel_size, relax);
			
			/* Preorder the matrix, obtain the column etree. */
			sp_preorder("N", &A, perm_c, etree, &AC);

			/* Factor the matrix AC. */
			cgstrf("N", &AC, iparam.diag_pivot_thresh,
			       iparam.drop_tol, iparam.relax,
			       iparam.panel_size, etree,
			       work, lwork, perm_r, perm_c, &L, &U, &info);

			if ( info ) { 
                            printf("** First factor: info %d, equed %c\n",
				   info, *equed);
                            if ( lwork == -1 ) {
                                printf("** Estimated memory: %d bytes\n",
                                        info - n);
                                exit(0);
                            }
                        }
	
                        Destroy_CompCol_Permuted(&AC);
			StatFree();
			
		    } /* if .. first time factor */
		    
		    for (itran = 0; itran < NTRAN; ++itran) {
			*trans = transs[itran];

			/* Restore the matrix A. */
			cCopy_CompCol_Matrix(&ASAV, &A);
			
 			/* Set the right hand side. */
			cFillRHS(trans, nrhs, xact, ldx, &A, &B);
			cCopy_Dense_Matrix(m, nrhs, rhsb, ldb, bsav, ldb);

			/*----------------
			 * Test cgssv
			 *----------------*/
			if ( nofact && norefact && itran == 0) {
                            /* Not yet factored, and untransposed */
	
			    cCopy_Dense_Matrix(m, nrhs, rhsb, ldb, solx, ldx);
			    cgssv(&A, perm_c, perm_r, &L, &U, &X, &info);
			    
			    if ( info && info != izero ) {
                                printf(FMT3, "cgssv",
				       info, izero, n, nrhs, imat, nfail);
			    } else {
                                /* Reconstruct matrix from factors and
	                           compute residual. */
                                sp_cget01(m, n, &A, &L, &U, perm_r,
					       &result[0]);
				nt = 1;
				if ( izero == 0 ) {
				    /* Compute residual of the computed
				       solution. */
				    cCopy_Dense_Matrix(m, nrhs, rhsb, ldb,
						       wwork, ldb);
				    sp_cget02(trans, m, n, nrhs, &A, solx,
                                              ldx, wwork,ldb, &result[1]);
				    nt = 2;
				}
				
				/* Print information about the tests that
				   did not pass the threshold.      */
				for (i = 0; i < nt; ++i) {
				    if ( result[i] >= THRESH ) {
					printf(FMT1, "cgssv", n, i,
					       result[i]);
					++nfail;
				    }
				}
				nrun += nt;
			    } /* else .. info == 0 */

			    /* Restore perm_c. */
			    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i];

		            if (lwork == 0) {
			        Destroy_SuperNode_Matrix(&L);
			        Destroy_CompCol_Matrix(&U);
			    }
			} /* if .. end of testing cgssv */
    
			/*----------------
			 * Test cgssvx
			 *----------------*/
    
			/* Equilibrate the matrix if fact = 'F' and
			   equed = 'R', 'C', or 'B'.   */
			if ( iequed > 0 && n > 0 ) {
			    claqgs(&A, R, C, rowcnd, colcnd, amax, equed);
			}
			
			/* Solve the system and compute the condition number
			   and error bounds using cgssvx.      */
			cgssvx(fact, trans, refact, &A, &iparam, perm_c,
			       perm_r, etree, equed, R, C, &L, &U, work,
			       lwork, &B, &X, &rpg, &rcond, ferr, berr,
			       &mem_usage, &info);

			if ( info && info != izero ) {
			    printf(FMT3, "cgssvx",
				   info, izero, n, nrhs, imat, nfail);
                            if ( lwork == -1 ) {
                                printf("** Estimated memory: %d bytes\n",
                                        mem_usage.total_needed);
                                exit(0);
                            }
			} else {
			    if ( !prefact ) {
			    	/* Reconstruct matrix from factors and
	 			   compute residual. */
                                sp_cget01(m, n, &A, &L, &U, perm_r,
					       &result[0]);
				k1 = 0;
			    } else {
			   	k1 = 1;
			    }

			    if ( !info ) {
				/* Compute residual of the computed solution.*/
				cCopy_Dense_Matrix(m, nrhs, bsav, ldb,
						  wwork, ldb);
				sp_cget02(trans, m, n, nrhs, &ASAV, solx, ldx,
					  wwork, ldb, &result[1]);

				/* Check solution from generated exact
				   solution. */
				sp_cget04(n, nrhs, solx, ldx, xact, ldx, rcond,
					  &result[2]);

				/* Check the error bounds from iterative
				   refinement. */
				sp_cget07(trans, n, nrhs, &ASAV, bsav, ldb,
					  solx, ldx, xact, ldx, ferr, berr,
					  &result[3]);

				/* Print information about the tests that did
				   not pass the threshold.    */
				for (i = k1; i < NTESTS; ++i) {
				    if ( result[i] >= THRESH ) {
					printf(FMT2, "cgssvx",
					       *fact, *trans, *refact, *equed,
					       n, imat, i, result[i]);
					++nfail;
				    }
				}
				nrun += NTESTS;
			    } /* if .. info == 0 */
			} /* else .. end of testing cgssvx */

		    } /* for itran ... */

		    if ( lwork == 0 ) {
			Destroy_SuperNode_Matrix(&L);
			Destroy_CompCol_Matrix(&U);
		    }

		} /* for irefact ... */
	    } /* for ifact ... */
	} /* for iequed ... */
#if 0    
    if ( !info ) {
	PrintPerf(&L, &U, &mem_usage, rpg, rcond, ferr, berr, equed);
    }
#endif    

    } /* for imat ... */

    /* Print a summary of the results. */
    PrintSumm("CGE", nfail, nrun, nerrs);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (bsav);
    SUPERLU_FREE (solx);    
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (pc_save);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    SUPERLU_FREE (rwork);
    SUPERLU_FREE (wwork);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    Destroy_CompCol_Matrix(&A);
    Destroy_CompCol_Matrix(&ASAV);
    if ( lwork > 0 ) {
	SUPERLU_FREE (work);
	Destroy_SuperMatrix_Store(&L);
	Destroy_SuperMatrix_Store(&U);
    }

    return 0;
}

/*  
 * Parse command line options to get relaxed snode size, panel size, etc.
 */
void
parse_command_line(int argc, char *argv[], char *matrix_type,
		   int *n, int *w, int *relax, int *nrhs,
		   int *maxsuper, int *rowblk, int *colblk, int *lwork)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "ht:n:w:r:s:m:b:c:l:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-r <int> - granularity of relaxed supernodes\n");
	    exit(1);
	    break;
	  case 't': strcpy(matrix_type, optarg);
	            break;
	  case 'n': *n = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'r': *relax = atoi(optarg); 
	            break;
	  case 's': *nrhs = atoi(optarg); 
	            break;
	  case 'm': *maxsuper = atoi(optarg); 
	            break;
	  case 'b': *rowblk = atoi(optarg); 
	            break;
	  case 'c': *colblk = atoi(optarg); 
	            break;
	  case 'l': *lwork = atoi(optarg); 
	            break;
  	}
    }
}
