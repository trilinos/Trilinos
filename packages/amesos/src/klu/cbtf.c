#include "klu_kernel.h"

/* the btf fortran interface to mc13 and mc21 */
int btf_ (int *Ai, int *Ap, int *nz, int *n, int *A1p,
    int *A1i, int *Perm, int *Lenc, int *Cp, int *Bp,
    int *Cperm, int *Rperm, int *Work, int *nblocks) ;

int cbtf
(
    int n,
    int Ap [ ],
    int Ai [ ],
    int P [ ],
    int Q [ ],
    int R [ ]
)
{
    int nz, *A1p, *A1i, *Perm, *Cperm, *Rperm, *Lenc, *Work, *Cp, *Bp,
	nblocks, k ;

    nblocks = 0 ;
    nz = Ap [n] ;
    if (nz <= 0)
    {
	/* singular */
	return (nblocks) ;
    }

    /* allocate workspace */
    A1p = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    A1i = (int *) ALLOCATE ((nz+1) * sizeof (int)) ;
    Lenc  = (int *) ALLOCATE (n * sizeof (int)) ;
    Perm  = (int *) ALLOCATE (n * sizeof (int)) ;
    Cp    = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Work  = (int *) ALLOCATE ((4*n) * sizeof (int)) ;

    /* allocate 1-based outputs */
    Cperm = (int *) ALLOCATE (n * sizeof (int)) ;
    Rperm = (int *) ALLOCATE (n * sizeof (int)) ;
    Bp    = (int *) ALLOCATE (n * sizeof (int)) ;

    if ((A1p == (int *) NULL) || (A1i == (int *) NULL) ||
	(Lenc == (int *) NULL) || (Perm == (int *) NULL) ||
	(Cp == (int *) NULL) || (Work == (int *) NULL) ||
	(Cperm == (int *) NULL) || (Rperm == (int *) NULL) ||
	(Bp == (int *) NULL))
    {
	/* out of memory */
	nblocks = -1 ;
    }

    if (nblocks != -1)
    {
	(void) btf_ (Ai, Ap, &nz, &n, A1p, A1i, Perm, Lenc, Cp, Bp,
	    Cperm, Rperm, Work, &nblocks) ;
    }

    /* free workspace */
    FREE (A1p, int) ;
    FREE (A1i, int) ;
    FREE (Lenc, int) ;
    FREE (Work, int) ;
    FREE (Cp, int) ;
    FREE (Perm, int) ;

    /* convert 1-based outputs to 0-based outputs (if successful) */
    if (nblocks > 0)
    {
	for (k = 0 ; k < n ; k++)
	{
	    P [k] = Rperm [k] - 1 ;
	    Q [k] = Cperm [k] - 1 ;
	}
	for (k = 0 ; k < nblocks ; k++)
	{
	    R [k] = Bp [k] - 1 ;
	}
	R [nblocks] = n ;
    }

    /* free the 1-based outputs */
    FREE (Cperm, int) ;
    FREE (Rperm, int) ;
    FREE (Bp, int) ;

    return (nblocks) ;
}
