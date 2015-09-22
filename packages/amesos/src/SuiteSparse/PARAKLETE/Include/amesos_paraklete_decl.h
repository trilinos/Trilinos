/* ========================================================================== */
/* === Paraklete/paraklete.h ================================================ */
/* ========================================================================== */

/* Include file for Paraklete library API.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* use the UF_long version of CHOLMOD and KLU, where Int is UF_long ("long"
 * on all platforms except Windows, where it is int64 but where Paraklete has
 * not been ported) */

/* Include file for user programs that call paraklete_* routines */

#ifndef AMESOS_PARAKLETE_DECL_H
#define AMESOS_PARAKLETE_DECL_H

#ifndef DLONG
#define DLONG
#endif

#include "amesos_cholmod.h"
#include "amesos_cholmod_internal.h"
#include "amesos_btf_decl.h"
#include "amesos_klu_decl.h"

#ifndef NMPI
#include "mpi.h"
#define MPI(statement) statement
#else
#define MPI(statement)
#endif

/* all integers are UF_long, or "long" on all but Windows */
#define MPI_Int MPI_LONG

/* use the long version of KLU and BTF */
#define KLU_analyze amesos_klu_l_analyze
#define KLU_factor amesos_klu_l_factor
#define KLU_solve amesos_klu_l_solve
#define KLU_symbolic klu_l_symbolic
#define KLU_numeric klu_l_numeric
#define KLU_common klu_l_common
#define KLU_free_symbolic amesos_klu_l_free_symbolic
#define KLU_free_numeric amesos_klu_l_free_numeric
#define KLU_defaults amesos_klu_l_defaults

#define BTF_order amesos_btf_l_order

#define TAG0 0

/* turn off debugging */
#ifndef NDEBUG
#define NDEBUG
#endif

/* To enable debugging, uncomment this line */
/*
#undef NDEBUG
*/

/* ========================================================================== */
/* === paraklete_common ===================================================== */
/* ========================================================================== */

/* Workspace (for each node) and global parameters */

typedef struct paraklete_common_struct
{
    double tol_diag ;       /* diagonal pivot tolerance (0 to 1)
                             * default 0.01 */
    double tol_offdiag ;    /* off-diagonal pivot tolerance (0 to 1),
                             * default 1.0 */
    double growth ;         /* memory reallocation factor, default 2.0 */

    Int status ;            /* PK_OK, PK_SINGULAR, etc, ... (see below) */

/* the first two error codes match with CHOLMOD error codes: */
#define PK_OK CHOLMOD_OK                        /* status is OK */
#define PK_OUT_OF_MEMORY CHOLMOD_OUT_OF_MEMORY  /* out of memory */
#define PK_TOO_LARGE CHOLMOD_TOO_LARGE          /* Int overflow occurred */
#define PK_SINGULAR (-5)        /* singular */
#define PK_UNKNOWN (-6)         /* not computed or not known */

    Int nproc ;             /* number of processors */
    Int myid ;              /* who am I? (in range 0 to nproc-1) */

    Int nleaves ;           /* number of leaves in separator tree.  Should
                             * be >= nprocs for parallelism.  Larger value is
                             * OK (can affect ordering quality, up or down).
                             * Default: nprocs */

    Int dump ;		    /* debug print level, for debugging only */
    FILE *file ;	    /* per-process output file, for debugging only */

    cholmod_common cm ;	    /* CHOLMOD workspace for each node */

    KLU_common km ;         /* KLU common */
}
paraklete_common ;


/* ========================================================================== */
/* === paraklete_symbolic =================================================== */
/* ========================================================================== */

/* The LUsymbolic object consists of 3 malloc'd blocks: the LUsymbolic struct
 * itself, Mem_n of size 3*n, and Mem_c of size 7*ncomponents+2 */

typedef struct paraklete_symbolic_struct
{
    Int n ;		/* dimension of the matrix to factorize */

    Int ncomponents ;	/* final number of components in the separator tree. */
    Int ncomp0 ;        /* initial (from cholmod_nested_dissection) */
    Int ncomp1 ;        /* after collapsing (from cholmod_collapse_septree) */

    /* ---------------------------------------------------------------------- */
    /* size-3n memory space, in Mem_n: */
    /* ---------------------------------------------------------------------- */

    Int *Mem_n ;	/* size 3*n (malloc'd block for Cperm, RpermInv,
			 * Cparent */

    Int *Cperm ;	/* size n.  Symmetric nested-dissection fill-reducing
			 * ordering of the original matrix.  Cperm [k] = i if
			 * row/column i of the original matrix becomes row/
			 * column k of the permuted matrix.  Does not take into
			 * account numerical pivoting (as returned by
			 * paraklete_analyze).
			 * -----------------------------------------------------
			 * paraklete_reanalyze combines the fill-reducing
			 * ordering with column permutations arising from column
			 * pivoting in paraklete_factorize, and stores this
			 * as Cperm (a copy of LU->Q, from the numeric object).
			 */

    Int *RpermInv ;	/* size n.  The inverse of Rperm, or Cperm if Rperm is
			 * NULL.  RpermInv [i] = k if Rperm [k] = i. */

    Int *Cparent ;	/* size n, but only entries 0 to ncomponents-1 used.
			 * Cparent [c] is the parent of node c in the separator
			 * tree, or EMPTY if c is a root. */

    /* ---------------------------------------------------------------------- */
    /* size-(7*ncomponents+2) memory space, in Mem_c: */
    /* ---------------------------------------------------------------------- */

    Int *Mem_c ;	/* size 7*ncomponents+2 (malloc'd block for Cstart,
			 * Child, Childp, Clnz, Cn, Cnz, and Sched) */

    Int *Cstart ;	/* size ncomponents+1.  Cstart[c] is the first candidate
			 * pivot row/column in node c, in the matrix permuted
			 * with Cperm.  Cstart [c+1]-1 is the last candidate
			 * pivot row/column in node c. */

    Int *Child ;	/* size ncomponents. */

    Int *Childp ;	/* size ncomponents+1.
			 * Child [Childp [c] ... Childp [c+1]-1] are the
			 * children of node c.
			 */

    Int *Cn ;		/* size ncomponents.  Cn [c] is the total dimension of
		         * node c, including the pivots (Cstart [c] to
			 * Cstart [c+1]-1), and the sum of the number of pivots
			 * of all ancestors of node c in the separator tree.  If
			 * no pivots are lost during numerical factorization,
			 * Cn [c] is the dimension of the LU factors of node c.
			 */

    Int *Clnz ;		/* size ncomponents.  Clnz [c] is an estimate of the
			 * number of nonzeros in the LU factors of node c,
			 * assuming no numerical pivoting and assuming A+A' is
			 * factorized. */

    Int *Cnz ;		/* size ncomponents.  Cnz [c] is the number of entries
			 * in the original matrix A that start in node c. */

    Int *Sched ;	/* size ncomponents.  Sched [c] is the id of the process
			 * that factorizes node c */

    /* ---------------------------------------------------------------------- */
    /* size-n memory space, allocated by itself */
    /* ---------------------------------------------------------------------- */

    Int *Rperm ;	/* size n if not NULL.  NULL for output of
			 * paraklete_analyze (in which case Rperm and the
			 * initial Cperm are identical).
			 * -----------------------------------------------------
			 * For paraklete_reanalyze, this is a copy of LU->P from
			 * the LU numeric object from paraklete_factorize. */
}
paraklete_symbolic ;


/* ========================================================================== */
/* === paraklete_node ======================================================= */
/* ========================================================================== */

/* LU factorization and Schur complement for one node of the separator tree */

/* Return pointers to the integer and real parts of column j of a matrix A
 * stored with merged row indices and numerical values.
 * Assumes sizeof(double) = 2*sizeof(Int).
 */

#define GET_COLUMN(Ap,Anz,Aix,j,Ai,Ax,len) \
{ \
    Int getp = Ap [j] ; \
    len = Anz [j] ; \
    Ai = (Int *) (Aix + getp) ; \
    Ax = (double *) (Aix + getp + ((len + 1) / 2)) ; \
}

typedef struct paraklete_node_struct
{

    /* ---------------------------------------------------------------------- */
    /* header information */
    /* ---------------------------------------------------------------------- */

    /* The header is sent in a single message, from child to parent */

#define PK_HEADER 8		/* size of header */

    Int header [PK_HEADER] ;	/* contents defined below: */

#define PK_STATUS header[0]	/* status (see Common->status above) */

#define PK_SN header[1]		/* dimension of the Schur complement */

#define PK_SNZ header[2]	/* number of nonzeros in S */

#define PK_NFOUND header[3]	/* number of pivots found for this node. */

#define PK_NLOST header[4]	/* # pivots lost, equal to npiv - nfound */

#define PK_NN header[5]		/* dimension of the LU factors of this node. */
				/* Equal to Cn [c] c if no pivots fail. 
				 * Failed pivot candidates of the children are
				 * added to Cn [c] to obtain nn.  */

#define PK_NPIV header[6]	/* number of candidate pivots for this node. */
				/* Equal to Cstart [c+1] - Cstart [c] if no
				 * pivots fail. Failed pivots from the
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

    Int lnz ;		/* number of nonzeros in L */

    Int unz ;		/* number of nonzeros in U */

    Int nlost_in ;	/* sum of # of pivots lost in children */

    Int *llen ;		/* size npiv.  llen [j] is the number of entries in
			 * column j of L. */

    Int *ulen ;		/* size npiv.  ulen [j] is the number of entries in
			 * column j of U */

    Int *lp ;		/* size nn. lp [j] is the start of column j of L in the
			 * ix array. */

    Int *up ;		/* size nn.  up [j] is the start of column j of U in the
			 * ix array */

    Int *Plocal ;	/* size npiv. Plocal [k] = i if local row i is the kth
			 * row of the LU factors for this node.  Both k and i
			 * are in the range 0 to npiv-1 */

    Int *Qlocal ;	/* size npiv.
			 * Qlocal [k] = j if local col j is the kth col of the
			 * LU factors for this node. */

    Int *Pglobal ;	/* size npiv. Pglobal [k] = i if local row i is the kth
			 * row of the LU factors of A(Cperm,Cperm).
			 * k is in the range 0 to npiv-1,
			 * i is in the range 0 to (LU->n)-1 */

    Int *Qglobal ;	/* size npiv. Qglobal [k] = j if local col j is the kth
			 * col of the LU factors of A(Cperm,Cperm).
			 * k is in the range 0 to npiv-1,
			 * j is in the range 0 to (LU->n)-1 */

    Int *Pinv ;		/* size npiv.  The inverse of Plocal.
			 * Pinv [i] = k if local row i is the kth row of the
			 * LU factors for this node. */

    Int *Qinv ;		/* size npiv.  The inverse of Qlocal.
			 * Qinv [j] = k if local col j is the kth col of the
			 * LU factors for this node. */

    Int lusize ;	/* size of ix array */

    double *ix ;	/* indices and numerical values of the LU factors */

    Int *Lost ;		/* size nchild */
    Int *Lostp ;	/* size nchild+1 */
    Int nchild ;

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

    Int nrealloc ;
    Int noffdiag ;
    double umin ;
    double umax ;

    /* ---------------------------------------------------------------------- */
    /* Schur complement */
    /* ---------------------------------------------------------------------- */

    /* The Schur complement has dimension sn = nn-nfound */

    double *sx ;	/* indices and numerical values of Schur complmement */
    Int *sp ;		/* size sn.  sp [j] is the start of column j in sx */
    Int *slen ;		/* size sn.  sp [j] is the start of column j in sx */

    /* ---------------------------------------------------------------------- */
    /* solution to Ly=Pb and Ux=y */
    /* ---------------------------------------------------------------------- */

    Int nk ;		/* Cstart [c+1] - Cstart [c], nominal pivots at node */
    double *B ;		/* size nk = Cstart [c+1] - Cstart [c] */
    double *X ;		/* size nn */

    /* ---------------------------------------------------------------------- */
    /* MPI message handling */
    /* ---------------------------------------------------------------------- */

    MPI (MPI_Request req ; )
    MPI (MPI_Request *Req ; )	/* size nchild */

    /* workspace at this node */
    Int *W2 ;			/* size 2*nlost_in */

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

#define PARAKLETE_MAGIC (-42)

typedef struct paraklete_numeric_struct
{
    Int magic ;     /* always equal to PARAKLETE_MAGIC, and less than zero */
    Int n ;
    Int ncomponents ;

    Int *P ;	    /* size n. row permutation */
    Int *Q ;	    /* size n. column permutation.  Both P and Q include
		     * LUsymbolic->Cperm, so that L*U=P*A*Q is the final
		     * factorization. */

    Int *Pinv ;
    Int *Qinv ;

    paraklete_node **LUnode ;	    /* size ncomponents.  LUnode [c] is a
				     * pointer to the LU factorization of
				     * node c. */

    /* workspace for process 0 only */
    double *W ;		/* size n */
    Int *Ep2 ;		/* size n+1 */
    cholmod_sparse *E ;	/* n-by-n, E = A (Cperm, Cperm) */

} paraklete_numeric ;

/* ========================================================================== */
/* === paraklete_btf_symbolic =============================================== */
/* ========================================================================== */

typedef struct paraklete_btf_symbolic_struct
{
    Int n ;		/* dimension of A */
    Int nblocks ;	/* number of diagonal blocks in BTF form */
    Int cnz ;		/* # of entries in diagonal blocks of A(p,q) */
    Int fnz ;		/* # of entries in off-diagonal blocks of A(p,q) */
    Int *Mem_n ;	/* contains Qbtf, Pbinv, and Rbtf */
    /* Int *Pbtf ; */	/* BTF row permutation, size n */
    Int *Qbtf ;		/* BTF column permutation, size n */
    Int *Pbinv ;	/* inverse of Pbtf, size n */
    Int *Rbtf ;		/* BTF block boundaries, size n+1 (0..nblocks used) */

    /* symbolic analysis of each diagaonal block (NULL for singletons) */
    void **LUsymbolic ;	/* array of pointers, size n */
			/* only first nblocks entries used */
}
paraklete_btf_symbolic ;

/* ========================================================================== */
/* === paraklete_btf_numeric ================================================ */
/* ========================================================================== */

typedef struct paraklete_btf_numeric_struct
{
    Int nblocks ;	    /* number of diagonal blocks in BTF form */
    double *Singleton ;	    /* singleton values, size nblocks */

    cholmod_sparse *F ;	    /* off-diagonal entries, size n-by-n */

    /* symbolic analysis of each diagaonal block (NULL for singletons) */
    void **LUnumeric ;      /* array of pointers, size nblocks */
}
paraklete_btf_numeric ;

/* ========================================================================== */

paraklete_btf_symbolic *amesos_paraklete_btf_analyze
(
    cholmod_sparse *A,
    paraklete_common *Common
) ;

paraklete_btf_numeric *amesos_paraklete_btf_factorize
(
    cholmod_sparse *A,
    paraklete_btf_symbolic *LU_btf_symbolic,
    paraklete_common *Common
) ;

Int amesos_paraklete_btf_solve
(
    paraklete_btf_numeric *LU_btf_numeric,
    paraklete_btf_symbolic *LU_btf_symbolic,
    double *B,
    paraklete_common *Common
) ;

void amesos_paraklete_btf_free_symbolic
(
    paraklete_btf_symbolic **LU_btf_symbolic_handle,
    paraklete_common *Common
) ;

void amesos_paraklete_btf_free_numeric
(
    paraklete_btf_numeric **LU_btf_numeric_handle,
    paraklete_common *Common
) ;

paraklete_btf_symbolic *amesos_paraklete_btf_alloc_symbolic
(
    Int n,
    Int nblocks,
    paraklete_common *Common
) ;

/* ========================================================================== */

paraklete_symbolic *amesos_paraklete_analyze
(
    cholmod_sparse *A,
    paraklete_common *Common
) ;

paraklete_numeric *amesos_paraklete_factorize
(
    cholmod_sparse *A,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

Int amesos_paraklete_kernel
(
    cholmod_sparse *A,
    paraklete_node *LUnode,
    paraklete_common *Common
) ;

Int amesos_paraklete_factorize_node
(
    Int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

Int amesos_paraklete_solve
(
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    double *B,
    paraklete_common *Common
) ;

Int amesos_paraklete_lsolve_node
(
    Int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

Int amesos_paraklete_usolve_node
(
    Int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

paraklete_symbolic *amesos_paraklete_alloc_symbolic
(
    Int n,
    Int ncomponents,
    Int do_Rperm,
    paraklete_common *Common
) ;

void amesos_paraklete_free_symbolic
(
    paraklete_symbolic **LUsymbolicHandle,
    paraklete_common *Common
) ;

void amesos_paraklete_free_numeric
(
    paraklete_numeric **LUHandle,
    paraklete_common *Common
) ;

paraklete_symbolic *amesos_paraklete_reanalyze
(
    cholmod_sparse *A,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
) ;

void amesos_paraklete_start (Int nproc, Int myid, paraklete_common *Common) ;
void amesos_paraklete_finish (paraklete_common *Common) ;

#define PARAKLETE_ERROR(status,message) \
    amesos_paraklete_error (status, __FILE__, __LINE__, message, Common) ;

void amesos_paraklete_error (Int status, char *filename, Int line, char *message,
    paraklete_common *Common) ;

/* ========================================================================== */
/* === Paraklete version ==================================================== */
/* ========================================================================== */

#define PARAKLETE_DATE "Nov 27, 2007"
#define PARAKLETE_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define PARAKLETE_MAIN_VERSION 0
#define PARAKLETE_SUB_VERSION 3
#define PARAKLETE_VERSION \
    PARAKLETE_VERSION_CODE(PARAKLETE_MAIN_VERSION,PARAKLETE_SUB_VERSION)

/* ========================================================================== */
/* === debugging definitions ================================================ */
/* ========================================================================== */

#undef DEBUG
#undef ASSERT
#undef PR0
#undef PR1
#undef PR2
#undef PR3

extern Int my_tries ;

#ifndef NDEBUG

#include <assert.h>
#define PR0(params) { (void) fprintf params ; fflush (Common->file) ; }
#define PR1(params) { if (Common->dump >= 1) (void) fprintf params ; fflush (Common->file) ; }
#define PR2(params) { if (Common->dump >= 2) (void) fprintf params ; fflush (Common->file) ; }
#define PR3(params) { if (Common->dump >= 3) (void) fprintf params ; fflush (Common->file) ; }
#define DEBUG(statement) statement
#define ASSERT(expression) (assert (expression))

#else

#define PR0(params)
#define PR1(params)
#define PR2(params)
#define PR3(params)
#define DEBUG(statement)
#define ASSERT(expression)

#endif

#endif /* AMESOS_PARAKLETE_DECL_H */
