/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include "superlu_zdefs.h"


/* Variables external to this file */
extern LU_stack_t stack;



/*
 * mem_usage consists of the following fields:
 *    - for_lu (float)
 *      The amount of space used in bytes for the L\U data structures.
 *    - total (float)
 *      The amount of space needed in bytes to perform factorization.
 *    - expansions (int)
 *      Number of memory expansions during the LU factorization.
 */
int_t zQuerySpace_dist(int_t n, LUstruct_t *LUstruct, gridinfo_t *grid,
		       mem_usage_t *mem_usage)
{
    register int_t dword, gb, iword, k, maxsup, nb, nsupers;
    int_t *index, *xsup;
    int iam, mycol, myrow;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;

    iam = grid->iam;
    myrow = MYROW( iam, grid );
    mycol = MYCOL( iam, grid );
    iword = sizeof(int_t);
    dword = sizeof(doublecomplex);
    maxsup = sp_ienv_dist(3);
    nsupers = Glu_persist->supno[n-1] + 1;
    xsup = Glu_persist->xsup;
    mem_usage->for_lu = 0;

    /* For L factor */
    nb = CEILING( nsupers, grid->npcol ); /* Number of local column blocks */
    for (k = 0; k < nb; ++k) {
	gb = k * grid->npcol + mycol; /* Global block number. */
	if ( gb < nsupers ) {
	    index = Llu->Lrowind_bc_ptr[k];
	    if ( index ) {
		mem_usage->for_lu += (float)
		    ((BC_HEADER + index[0]*LB_DESCRIPTOR + index[1]) * iword);
		mem_usage->for_lu += (float)(index[1]*SuperSize( gb )*dword);
	    }
	}
    }

    /* For U factor */
    nb = CEILING( nsupers, grid->nprow ); /* Number of local row blocks */
    for (k = 0; k < nb; ++k) {
	gb = k * grid->nprow + myrow; /* Global block number. */
	if ( gb < nsupers ) {
	    index = Llu->Ufstnz_br_ptr[k];
	    if ( index ) {
		mem_usage->for_lu += (float)(index[2] * iword);
		mem_usage->for_lu += (float)(index[1] * dword);
	    }
	}
    }

    /* Working storage to support factorization */
    mem_usage->total = mem_usage->for_lu;
    mem_usage->total +=
	(float)(( Llu->bufmax[0] + Llu->bufmax[2] ) * iword +
		( Llu->bufmax[1] + Llu->bufmax[3] + maxsup ) * dword );
    /**** another buffer to use mpi_irecv in pdgstrf_irecv.c ****/
    mem_usage->total +=
	(float)( Llu->bufmax[0] * iword +  Llu->bufmax[1] * dword );
    mem_usage->total += (float)( maxsup * maxsup + maxsup) * iword;
    k = CEILING( nsupers, grid->nprow );
    mem_usage->total += (float)(2 * k * iword);

    return 0;
} /* zQuerySpace_dist */


/*
 * Allocate storage for original matrix A
 */
void
zallocateA_dist(int_t n, int_t nnz, doublecomplex **a, int_t **asub, int_t **xa)
{
    *a    = (doublecomplex *) doublecomplexMalloc_dist(nnz);
    *asub = (int_t *) intMalloc_dist(nnz);
    *xa   = (int_t *) intMalloc_dist(n+1);
}


doublecomplex *doublecomplexMalloc_dist(int_t n)
{
    doublecomplex *buf;
    buf = (doublecomplex *) SUPERLU_MALLOC(n * sizeof(doublecomplex)); 
    return (buf);
}

doublecomplex *doublecomplexCalloc_dist(int_t n)
{
    doublecomplex *buf;
    register int_t i;
    doublecomplex zero = {0.0, 0.0};
    buf = (doublecomplex *) SUPERLU_MALLOC(n * sizeof(doublecomplex));
    if ( !buf ) return (buf);
    for (i = 0; i < n; ++i) buf[i] = zero;
    return (buf);
}

