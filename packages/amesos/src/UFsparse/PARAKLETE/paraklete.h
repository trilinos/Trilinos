/* ========================================================================== */
/* === Paraklete/paraklete.h ================================================ */
/* ========================================================================== */

/*
 * PARAKLETE version 0.2: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/*
#define HACK 8
*/

#include "cholmod.h"
#include "cholmod_internal.h"

#ifndef NMPI
#include "mpi.h"
#define MPI(statement) statement
#else
#define MPI(statement)
#endif

#define TAG0 0

/* To enable debugging, uncomment this line */
/* #undef NDEBUG  */

/* ========================================================================== */
/* === paraklete_common ===================================================== */
/* ========================================================================== */

/* Workspace (for each node) and global parameters */

typedef struct paraklete_common_struct
{
    double tol_diag ;
    double tol_offdiag ;
    double growth ;

    int status ;

    int nproc ;
    int myid ;
    int dump ;		/* debug print level */

    FILE *file ;	/* per-process output file, for debugging only */

    cholmod_common cm ;	/* workspace for each node */
}
paraklete_common ;


/* ========================================================================== */
/* === paraklete_symbolic =================================================== */
/* ========================================================================== */

/* The LUsymbolic object consists of 3 malloc'd blocks: the LUsymbolic struct
 * itself, Mem_n of size 3*n, and Mem_c of size 7*ncompenents+2 */

typedef struct paraklete_symbolic_struct
{
    int n ;		/* dimension of the matrix to factorize */

    int ncomponents ;	/* number of components in the separator tree. */

    int *Mem_n ;	/* size 3*n (malloc'd block for Cperm, Cinv, Cparent */

    int *Mem_c ;	/* size 7*ncomponents+2 (malloc'd block for Cstart,
			 * Child, Childp, Clnz, Cn, Cnz, and Sched) */

    int *Cperm ;	/* size n.  Symmetric nested-dissection fill-reducing
			 * ordering of the original matrix.  Cperm [k] = i if
			 * row/column i of the original matrix becomes row/
			 * column k of the permuted matrix.  Does not take into
			 * account numerical pivoting. */

    int *Cinv ;		/* size n.  The inverse of Cperm.  Cinv [i] = k if
			 * Cperm [k] = i. */

    int *Cparent ;	/* size n, but only entries 0 to ncomponents-1 used.
			 * Cparent [c] is the parent of node c in the separator
			 * tree, or EMPTY if c is a root. */

    int *Cstart ;	/* size ncomponents+1.  Cstart[c] is the first candidate
			 * pivot row/column in node c, in the matrix permuted
			 * with Cperm.  Cstart [c+1]-1 is the last candidate
			 * pivot row/column in node c. */

    int *Child ;	/* size ncomponents. */
    int *Childp ;	/* size ncomponents+1.
			 * Child [Childp [c] ... Childp [c+1]-1] are the
			 * children of node c.
			 */

    int *Cn ;		/* size ncomponents.  Cn [c] is the total dimension of
		         * node c, including the pivots (Cstart [c] to
			 * Cstart [c+1]-1), and the sum of the number of pivots
			 * of all ancestors of node c in the separator tree.  If
			 * no pivots are lost during numerical factorization,
			 * Cn [c] is the dimension of the LU factors of node c.
			 */

    int *Clnz ;		/* size ncomponents.  Clnz [c] is an estimate of the
			 * number of nonzeros in the LU factors of node c,
			 * assuming no numerical pivoting and assuming A+A' is
			 * factorized. */

    int *Cnz ;		/* size ncomponents.  Cnz [c] is the number of entries
			 * in the original matrix A that start in node c. */

    int *Sched ;	/* size ncomponents.  Sched [c] is the id of the process
			 * that factorizes node c */
}
paraklete_symbolic ;


/* ========================================================================== */
/* === paraklete_node ======================================================= */
/* ========================================================================== */

/* LU factorization and Schur complement for one node of the separator tree */

/* Return pointers to the integer and real parts of column j of a matrix A
 * stored with merged row indices and numerical values.
 * Assumes sizeof(double) = 2*sizeof(int).
 */

#define GET_COLUMN(Ap,Anz,Aix,j,Ai,Ax,len) \
{ \
    int getp = Ap [j] ; \
    len = Anz [j] ; \
    Ai = (int *) (Aix + getp) ; \
    Ax = (double *) (Aix + getp + ((len + 1) / 2)) ; \
}

typedef struct paraklete_node_struct
{

    /* ---------------------------------------------------------------------- */
    /* header information */
    /* ---------------------------------------------------------------------- */

    /* The header is sent in a single message, from child to parent */

#define PK_HEADER 8		/* size of header */

    int header [PK_HEADER] ;	/* contents defined below: */

#define PK_STATUS header[0]	/* status, one of the following states: */
#define PK_OK 0			/* OK */
#define PK_SINGULAR (-1)	/* singular (only applies to a root node of
				 * the separator tree) */
#define PK_OUT_OF_MEMORY (-2)	/* out of memory */
#define PK_UNKNOWN (-3)		/* not computed or not known */

#define PK_SN header[1]		/* dimension of the Schur complement */

#define PK_SNZ header[2]	/* number of nonzeros in S */

#define PK_NFOUND header[3]	/* number of pivots found for this node. */

#define PK_NLOST header[4]	/* # pivots lost, equal to npiv - nfound */

#define PK_NN header[5]		/* dimension of the LU factors of this node.
				 * Equal to Cn [c] c if no pivots fail. 
				 * Failed pivot candidates of the children are
				 * added to Cn [c] to obtain nn.  */

#define PK_NPIV header[6]	/* number of candidate pivots for this node.
				 * Equal to Cstart [c+1] - Cstart [c] if no
				 * pivots fail. * Failed pivots from the
				 * children are added to obtain npiv. */

#define PK_SSIZE header[7]	/* size of the sx array */

    /* ---------------------------------------------------------------------- */
    /* local copy of the matrix to factorize for just this node */
    /* ---------------------------------------------------------------------- */

    cholmod_sparse *A ;	/* dimension Cn [c] */
    cholmod_sparse *C ;	/* dimension nn, = A + Schur complements of children */

    /* ---------------------------------------------------------------------- */
    /* LU factors */
    /* ---------------------------------------------------------------------- */

    int lnz ;		/* number of nonzeros in L */

    int unz ;		/* number of nonzeros in U */

    int nlost_in ;	/* sum of # of pivots lost in children */

    int *llen ;		/* size npiv.  llen [j] is the number of entries in
			 * column j of L. */

    int *ulen ;		/* size npiv.  ulen [j] is the number of entries in
			 * column j of U */

    int *lp ;		/* size nn. lp [j] is the start of column j of L in the
			 * ix array. */

    int *up ;		/* size nn.  up [j] is the start of column j of U in the
			 * ix array */

    int *Plocal ;	/* size npiv. Plocal [k] = i if local row i is the kth
			 * row of the LU factors for this node.  Both k and i
			 * are in the range 0 to npiv-1 */

    int *Qlocal ;	/* size npiv.
			 * Qlocal [k] = j if local col j is the kth col of the
			 * LU factors for this node. */

    int *Pglobal ;	/* size npiv. Pglobal [k] = i if local row i is the kth
			 * row of the LU factors of A(Cperm,Cperm).
			 * k is in the range 0 to npiv-1,
			 * i is in the range 0 to (LU->n)-1 */

    int *Qglobal ;	/* size npiv. Qglobal [k] = j if local col j is the kth
			 * col of the LU factors of A(Cperm,Cperm).
			 * k is in the range 0 to npiv-1,
			 * j is in the range 0 to (LU->n)-1 */

    int *Pinv ;		/* size npiv.  The inverse of Plocal.
			 * Pinv [i] = k if local row i is the kth row of the
			 * LU factors for this node. */

    int *Qinv ;		/* size npiv.  The inverse of Qlocal.
			 * Qinv [j] = k if local col j is the kth col of the
			 * LU factors for this node. */

    int lusize ;	/* size of ix array */

    double *ix ;	/* indices and numerical values of the LU factors */

    int *Lost ;		/* size nchild */
    int *Lostp ;	/* size nchild+1 */
    int nchild ;

    /* To traverse column j of L:
     *
     *	    GET_COLUMN (LUnode->lp, LUnode->llen, LUnode->ix, Li, Lx, len) ;
     *	    for (p = 0 ; p < len ; p++)
     *	    {
     *		i = Li [p] ;
     *		lij = Lx [p] ;
     *	    }
     *
     * To traverse column j of U:
     *
     *	    GET_COLUMN (LUnode->up, LUnode->ulen, LUnode->ix, Ui, Ux, len) ;
     *	    for (p = 0 ; p < len ; p++)
     *	    {
     *		i = Ui [p] ;
     *		uij = Ux [p] ;
     *	    }
     */

    int nrealloc ;
    int noffdiag ;
    double umin ;
    double umax ;

    /* ---------------------------------------------------------------------- */
    /* Schur complement */
    /* ---------------------------------------------------------------------- */

    /* The Schur complement has dimension sn = nn-nfound */

    double *sx ;	/* indices and numerical values of Schur complmement */
    int *sp ;		/* size sn.  sp [j] is the start of column j in sx */
    int *slen ;		/* size sn.  sp [j] is the start of column j in sx */

    /* ---------------------------------------------------------------------- */
    /* solution to Ly=Pb and Ux=y */
    /* ---------------------------------------------------------------------- */

    int nk ;		/* Cstart [c+1] - Cstart [c], nominal pivots at node */
    double *B ;		/* size nk = Cstart [c+1] - Cstart [c] */
    double *X ;		/* size nn */

    /* ---------------------------------------------------------------------- */
    /* MPI message handling */
    /* ---------------------------------------------------------------------- */

    MPI (MPI_Request req ; )
    MPI (MPI_Request *Req ; )	/* size nchild */

    /* workspace at this node */
    int *W2 ;			/* size 2*nlost_in */

} paraklete_node ;


/* ========================================================================== */
/* === paraklete_numeric ==================================================== */
/* ========================================================================== */

/* All nodes share this data structure.  LUnode [0..ncomponents] is an array
 * of the LU factorization of each node, and status information about that node.
 * If a process owns node c, it uses LUnode [c] to hold permanent parts of the
 * global LU factors.  If it does not own node c, it uses LUnode [c] to hold
 * information about a child or parent in the separator tree (such as the
 * Schur complement from a child node owned by another process.
 *
 * Some information is shared globally: n, ncomponents, and the global
 * permutation (P, Q, Pinv, Qinv).   This structure also contains workspace
 * held by process 0 only (W, Ep2, and E).
 */

typedef struct paraklete_numeric_struct
{
    int n ;
    int ncomponents ;

    int *P ;	    /* size n. row permutation */
    int *Q ;	    /* size n. column permutation.  Both P and Q include
		     * LUsymbolic->Cperm, so that L*U=P*A*Q is the final
		     * factorization. */

    int *Pinv ;
    int *Qinv ;

    paraklete_node **LUnode ;	    /* size ncomponents.  LUnode [c] is a
				     * pointer to the LU factorization of
				     * node c. */

    /* workspace for process 0 only */
    double *W ;		/* size n */
    int *Ep2 ;		/* size n+1 */
    cholmod_sparse *E ;	/* n-by-n, E = A (Cperm, Cperm) */

} paraklete_numeric ;

/* ========================================================================== */

paraklete_symbolic *paraklete_analyze
(
    cholmod_sparse *A,
    paraklete_common *Common
) ;

paraklete_numeric *paraklete_factorize
(
    /* inputs, not modified */
    cholmod_sparse *A,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

int paraklete_kernel
(
    cholmod_sparse *A,
    paraklete_node *LUnode,
    paraklete_common *Common
) ;

int paraklete_factorize_node
(
    int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

int paraklete_solve
(
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    double *B,
    paraklete_common *Common
) ;

int paraklete_lsolve_node
(
    int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

int paraklete_usolve_node
(
    int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

void paraklete_free_symbolic
(
    paraklete_symbolic **LUsymbolicHandle,
    paraklete_common *Common
) ;

void paraklete_free_numeric
(
    paraklete_numeric **LUHandle,
    paraklete_common *Common
) ;

/* ========================================================================== */
/* === Paraklete version ==================================================== */
/* ========================================================================== */

#define PARAKLETE_DATE "May 23, 2006"
#define PARAKLETE_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define PARAKLETE_MAIN_VERSION 0
#define PARAKLETE_SUB_VERSION 2
#define PARAKLETE_VERSION \
    PARAKLETE_VERSION_CODE(PARAKLETE_MAIN_VERSION,PARAKLETE_SUB_VERSION)

/* ========================================================================== */
/* === debugging definitions ================================================ */
/* ========================================================================== */

#undef DEBUG

extern int my_tries ;

#ifndef NDEBUG

#define PR0(params) { (void) fprintf params ; fflush (Common->file) ; }
#define PR1(params) { if (Common->dump >= 1) (void) fprintf params ; fflush (Common->file) ; }
#define PR2(params) { if (Common->dump >= 2) (void) fprintf params ; fflush (Common->file) ; }
#define PR3(params) { if (Common->dump >= 3) (void) fprintf params ; fflush (Common->file) ; }
#define DEBUG(statement) statement

#else

#define PR0(params)
#define PR1(params)
#define PR2(params)
#define PR3(params)
#define DEBUG(statement)

#endif
