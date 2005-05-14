/* ========================================================================== */
/* === Include/cholmod_cholesky.h =========================================== */
/* ========================================================================== */

/* CHOLMOD Cholesky module.
 *
 * Sparse Cholesky routines: analysis, factorization, and solve.
 *
 * The primary routines are all that a user requires to order, analyze, and
 * factorize a sparse symmetric positive definite matrix A (or A*A'), and
 * to solve Ax=b (or A*A'x=b).  The primary routines rely on the secondary
 * routines, the CHOLMOD Core module, and the AMD and COLAMD packages.  They
 * make optional use of the CHOLMOD Supernodal and Partition modules, the
 * METIS package, and the CCOLAMD package.
 *
 * Primary routines:
 * -----------------
 *
 * cholmod_analyze		order and analyze (simplicial or supernodal)
 * cholmod_factorize		simplicial or supernodal Cholesky factorization
 * cholmod_solve		solve a linear system (simplicial or supernodal)
 *
 * Secondary routines:
 * ------------------
 *
 * cholmod_analyze_p		analyze, with user-provided permutation or f set
 * cholmod_factorize_p		factorize, with user-provided permutation or f
 * cholmod_etree		find the elimination tree
 * cholmod_rowcolcounts		compute the row/column counts of L
 * cholmod_amd			order using AMD
 * cholmod_colamd		order using COLAMD
 * cholmod_rowfac		incremental simplicial factorization
 * cholmod_resymbol		recompute the symbolic pattern of L
 * cholmod_resymbol_noperm	recompute the symbolic pattern of L, no L->Perm
 * cholmod_row			compute the pattern of kth row of L
 * 
 * Requires the Core module, and two packages: AMD and COLAMD.
 * Optionally uses the Supernodal and Partition modules.
 * Not required by any CHOLMOD module.
 */

#ifndef CHOLMOD_CHOLESKY_H
#define CHOLMOD_CHOLESKY_H

#include "cholmod_config.h"
#include "cholmod_core.h"

#ifndef NPARTITION
#include "cholmod_partition.h"
#endif

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

/* -------------------------------------------------------------------------- */
/* cholmod_analyze */
/* -------------------------------------------------------------------------- */

cholmod_factor *cholmod_analyze
(
    cholmod_sparse *A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_analyze_p */
/* -------------------------------------------------------------------------- */

cholmod_factor *cholmod_analyze_p
(
    cholmod_sparse *A,
    void *UserPerm,
    void *fset,
    size_t fsize,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_factorize */
/* -------------------------------------------------------------------------- */

int cholmod_factorize
(
    cholmod_sparse *A,
    cholmod_factor *L,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_factorize_p */
/* -------------------------------------------------------------------------- */

int cholmod_factorize_p
(
    cholmod_sparse *A,
    cholmod_scalar beta,
    void *fset,
    size_t fsize,
    cholmod_factor *L,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_solve */
/* -------------------------------------------------------------------------- */

#define CHOLMOD_A    0		/* solve Ax=b */
#define CHOLMOD_LDLt 1		/* solve LDL'x=b */
#define CHOLMOD_LD   2		/* solve LDx=b */
#define CHOLMOD_DLt  3		/* solve DL'x=b */
#define CHOLMOD_L    4		/* solve Lx=b */
#define CHOLMOD_Lt   5		/* solve L'x=b */
#define CHOLMOD_D    6		/* solve Dx=b */
#define CHOLMOD_P    7		/* permute x=Px */
#define CHOLMOD_Pt   8		/* permute x=P'x */

/* D is identity for LL' factorizations. */

int cholmod_solve
(
    int sys,
    cholmod_factor *L,
    cholmod_dense *X,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_etree */
/* -------------------------------------------------------------------------- */

int cholmod_etree
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in packed row form */

    /* outputs, not defined on input */
    void *Parent,	/* size nrow.  Parent [j] = p if p is the parent of j */

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_rowcolcounts */
/* -------------------------------------------------------------------------- */

int cholmod_rowcolcounts
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in packed column form */
    void *fset,		/* if present and A->stype=0, then LDL'=F*F' is
			 * analyzed, where F = A (:, fset [0..fsize-1]) */
    size_t fsize,		/* number of columns in the set f */
    void *Parent,	/* size nrow.  Parent [i] = p if p is the parent of i */
    void *Post,	/* size nrow.  Post [k] = i if i is the kth node in
			 * the postordered etree. */

    /* outputs, not defined on input */
    void *RowCount,	/* size nrow. RowCount [i] = # entries in the ith row of
			 * L, including the diagonal. */
    void *ColCount,	/* size nrow. ColCount [i] = # entries in the ith
			 * column of L, including the diagonal. */
    void *First,	/* size nrow.  First [i] = k is the least postordering
			 * of any descendant of i. */
    void *Level,	/* size nrow.  Level [i] is the length of the path from
			 * i to the root, with Level [root] = 0. */

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_amd */
/* -------------------------------------------------------------------------- */

int cholmod_amd
(
    cholmod_sparse *A,	    /* user's input matrix, n-by-n */
    void *Perm,		    /* size n, output permutation */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_colamd */
/* -------------------------------------------------------------------------- */

int cholmod_colamd
(
    cholmod_sparse *A,	    /* user's input matrix, nrow-by-ncol */
    void *fset,
    size_t fsize,
    void *Perm,		    /* size nrow, output permutation */
    cholmod_common *Common
) ;


/* -------------------------------------------------------------------------- */
/* cholmod_rowfac */
/* -------------------------------------------------------------------------- */

int cholmod_rowfac
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in packed column form */
    cholmod_sparse *F,	/* used for LDL'=AA' case only, stored in row form */
    cholmod_scalar beta,/* factorize beta*I+A or beta*I+AA' */
    size_t kstart,		/* first row to factorize */
    size_t kend,		/* last row to factorize is kend-1 */

    /* input/output: */
    cholmod_factor *L,

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_row_subtree */
/* -------------------------------------------------------------------------- */

/* Find the nonzero pattern of x for the system Lx=b where L = (0:k-1,0:k-1)
 * and b = kth column of A or A*A' (rows 0 to k-1 only) */

int cholmod_row_subtree
(
    /* inputs */
    cholmod_sparse *A,
    cholmod_sparse *F,	    /* F = A' */
    size_t krow,
    void *Parent_p,

    /* output */
    cholmod_sparse *R,	    /* pattern of L (k,:) */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_resymbol */
/* -------------------------------------------------------------------------- */

/* cholmod_resymbol is the same as cholmod_resymbol_noperm, except that it
 * first permutes A according to L->Perm.  A can be upper/lower/unsymmetric,
 * in contrast to cholmod_resymbol_noperm (which can be lower or unsym). */

int cholmod_resymbol
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in column form */
    void *fset,		/* if present and A->stype=0, then LDL'=F*F' is
			 * analyzed, where F = A (:, fset [0..fsize-1]).
			 * Otherwise, F=A. */
    size_t fsize,		/* number of columns in the set f */

    /* The factors L and D.  Some terms are inputs, some both input & output */
    cholmod_factor *L,	/* nrow-by-nrow */
    int pack,		/* if TRUE, convert L from unpacked to packed.  Note
			 * that a packed L remains packed and a dynamic L
			 * remains dynamic. */

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_resymbol_noperm */
/* -------------------------------------------------------------------------- */

int cholmod_resymbol_noperm
(
    cholmod_sparse *A,
    void *fset,
    size_t fsize,
    cholmod_factor *L,
    int pack,
    cholmod_common *Common
) ;

#endif
