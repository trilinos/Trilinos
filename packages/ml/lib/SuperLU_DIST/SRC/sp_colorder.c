#include "superlu_ddefs.h"


void check_perm(char *, int_t, int_t *);

void
sp_colorder(superlu_options_t *options,  SuperMatrix *A, int_t *perm_c, 
	    int_t *etree, SuperMatrix *AC)
{
/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 *
 * Purpose
 * =======
 *
 * sp_colorder() permutes the columns of the original matrix. It performs
 * the following steps:
 *
 *    1. Apply column permutation perm_c[] to A's column pointers to form AC;
 *
 *    2. If options->Fact = DOFACT, then
 *       (1) Compute column elimination tree etree[] of AC'AC;
 *       (2) Post order etree[] to get a postordered elimination tree etree[],
 *           and a postorder permutation post[];
 *       (3) Apply post[] permutation to columns of AC;
 *       (4) Overwrite perm_c[] with the product perm_c * post.
 *
 * Arguments
 * =========
 *
 * options (input) superlu_options_t*
 *         Specifies whether or not the elimination tree will be re-used.
 *         If options->Fact == DOFACT, this means first time factor A, 
 *         etree is computed and output.
 *         Otherwise, re-factor A, etree is input, unchanged on exit.
 *
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A can be:
 *         Stype = NC or NCP; Dtype = _D; Mtype = GE. In the future,
 *         more general A can be handled.
 *
 * perm_c  (input/output) int*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *         If options->Fact == DOFACT, perm_c is both input and output.
 *         On output, it is changed according to a postorder of etree.
 *         Otherwise, perm_c is input.
 *         
 * etree   (input/output) int*
 *         Elimination tree of Pc*(A'+A)*Pc', dimension A->ncol.
 *         If options->Fact == DOFACT, etree is an output argument,
 *         otherwise it is an input argument.
 *         Note: etree is a vector of parent pointers for a forest whose
 *         vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
 *
 * AC      (output) SuperMatrix*
 *         The resulting matrix after applied the column permutation
 *         perm_c[] to matrix A. The type of AC can be:
 *         Stype = NCP; Dtype = D; Mtype = GE.
 *
 */

    NCformat  *Astore;
    NCPformat *ACstore;
    int_t       *iwork, *post;
    register  int_t n, i;
#if ( DEBUGlevel>=1 )
    int iam;
    MPI_Comm_rank( MPI_COMM_WORLD, &iam );
    CHECK_MALLOC(iam, "Enter sp_colorder()");
#endif

    n = A->ncol;
    
    /* Apply column permutation perm_c to A's column pointers so to
       obtain NCP format in AC = A*Pc.  */
    AC->Stype       = NCP;
    AC->Dtype       = A->Dtype;
    AC->Mtype       = A->Mtype;
    AC->nrow        = A->nrow;
    AC->ncol        = A->ncol;
    Astore          = A->Store;
    ACstore = AC->Store = (void *) SUPERLU_MALLOC( sizeof(NCPformat) );
    if ( !ACstore ) ABORT("SUPERLU_MALLOC fails for ACstore");
    ACstore->nnz    = Astore->nnz;
    ACstore->nzval  = Astore->nzval;
    ACstore->rowind = Astore->rowind;
    ACstore->colbeg = (int_t*) SUPERLU_MALLOC(n*sizeof(int_t));
    if ( !(ACstore->colbeg) ) ABORT("SUPERLU_MALLOC fails for ACstore->colbeg");
    ACstore->colend = (int_t*) SUPERLU_MALLOC(n*sizeof(int_t));
    if ( !(ACstore->colend) ) ABORT("SUPERLU_MALLOC fails for ACstore->colend");

#if ( DEBUGlevel>=3 )
    if ( !iam ) {
	print_int_vec("pre_order:", n, perm_c);
	check_perm("Initial perm_c", n, perm_c);
    }
#endif      

    for (i = 0; i < n; i++) {
	ACstore->colbeg[perm_c[i]] = Astore->colptr[i]; 
	ACstore->colend[perm_c[i]] = Astore->colptr[i+1];
    }
	
    if ( options->Fact == DOFACT ) { 
	/* Factor A "from scratch" -- we also compute the etree, and
	 * make perm_c consistent with the postorder of the etree.
	 */

	iwork = (int_t*) SUPERLU_MALLOC((n+1)*sizeof(int_t)); 
	if ( !iwork ) ABORT("SUPERLU_MALLOC fails for iwork[]");

	if ( A->nrow != A->ncol ) { /* Rectangular matrix */
	    /* Compute the column etree of A*Pc'. */
	    sp_coletree(ACstore->colbeg, ACstore->colend, ACstore->rowind,
			A->nrow, A->ncol, etree);
	} else {
	    /* Compute the etree of Pc*(A'+A)*Pc'. */
	    int_t *b_colptr, *b_rowind, bnz, j;
	    int_t *c_colbeg, *c_colend;

	    /* Form B = A + A'. */
	    a_plus_at(n, Astore->nnz, Astore->colptr, Astore->rowind,
		      &bnz, &b_colptr, &b_rowind);

	    /* Form C = Pc*B*Pc'. */
	    c_colbeg = (int_t*) SUPERLU_MALLOC(n*sizeof(int_t));
	    c_colend = (int_t*) SUPERLU_MALLOC(n*sizeof(int_t));
	    if (!(c_colbeg) || !(c_colend) )
		ABORT("SUPERLU_MALLOC fails for c_colbeg/c_colend");
	    for (i = 0; i < n; i++) {
		c_colbeg[perm_c[i]] = b_colptr[i]; 
		c_colend[perm_c[i]] = b_colptr[i+1];
	    }
	    for (j = 0; j < n; ++j) {
		for (i = c_colbeg[j]; i < c_colend[j]; ++i) {
		    b_rowind[i] = perm_c[b_rowind[i]];
		}
	    }

	    /* Compute etree of C. */
	    sp_symetree(c_colbeg, c_colend, b_rowind, n, etree);

	    SUPERLU_FREE(b_colptr);
	    if ( bnz ) SUPERLU_FREE(b_rowind);
	    SUPERLU_FREE(c_colbeg);
	    SUPERLU_FREE(c_colend);
	}
#if ( DEBUGlevel>=3 )
	if ( !iam ) print_int_vec("etree:", n, etree);
#endif	
	
	/* Post order etree */
	post = (int_t *) TreePostorder(n, etree);
	/* for (i = 0; i < n+1; ++i) inv_post[post[i]] = i;
	   iwork = post; */

	/* Renumber etree in postorder */
	for (i = 0; i < n; ++i) iwork[post[i]] = post[etree[i]];
	for (i = 0; i < n; ++i) etree[i] = iwork[i];

#if ( DEBUGlevel>=3 )
	if ( !iam ) print_int_vec("postorder etree:", n, etree);
#endif
	
	/* Postmultiply A*Pc by post[] */
	for (i = 0; i < n; ++i) iwork[post[i]] = ACstore->colbeg[i];
	for (i = 0; i < n; ++i) ACstore->colbeg[i] = iwork[i];
	for (i = 0; i < n; ++i) iwork[post[i]] = ACstore->colend[i];
	for (i = 0; i < n; ++i) ACstore->colend[i] = iwork[i];

	for (i = 0; i < n; ++i)
	    iwork[i] = post[perm_c[i]];  /* product of perm_c and post */
	for (i = 0; i < n; ++i) perm_c[i] = iwork[i];

#if ( DEBUGlevel>=3 )
	if ( !iam ) {
	    print_int_vec("Pc*post:", n, perm_c);
	    check_perm("final perm_c", n, perm_c);
	}
#endif

	SUPERLU_FREE (post);
	SUPERLU_FREE (iwork);

    } /* if options->Fact == DOFACT ... */


#if ( DEBUGlevel>=1 )
    /* Memory allocated but not freed:
       ACstore, ACstore->colbeg, ACstore->colend  */
    CHECK_MALLOC(iam, "Exit sp_colorder()");
#endif

} /* SP_COLORDER */

void
check_perm(char *what, int_t n, int_t *perm)
{
    register int_t i;
    int_t          *marker;
    marker = (int_t *) intCalloc(n);

    for (i = 0; i < n; ++i) {
	if ( marker[perm[i]] == 1 || perm[i] >= n ) {
	    printf("%s: Not a valid PERM[%d] = %d\n", what, i, perm[i]);
	    ABORT("check_perm");
	} else {
	    marker[perm[i]] = 1;
	}
    }

    SUPERLU_FREE(marker);
}
