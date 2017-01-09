/* ========================================================================== */
/* === Include/cholmod_partition.h ========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_partition.h.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod_partition.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* CHOLMOD Partition module.
 *
 * Graph partitioning and graph-partition-based orderings.  Includes an
 * interface to CCOLAMD and CSYMAMD, constrained minimum degree ordering
 * methods which order a matrix following constraints determined via nested
 * dissection.
 *
 * cholmod_nested_dissection	CHOLMOD nested dissection ordering
 * cholmod_metis		METIS nested dissection ordering (METIS_NodeND)
 * cholmod_ccolamd		interface to CCOLAMD ordering
 * cholmod_csymamd		interface to CSYMAMD ordering
 * cholmod_camd			interface to CAMD ordering
 * cholmod_bisect		graph partitioner (currently based on METIS)
 * cholmod_metis_bisector	direct interface to METIS_NodeComputeSeparator
 *
 * Requires the Core and Cholesky modules, and three packages: METIS, CAMD,
 * and CCOLAMD.  Optionally used by the Cholesky module.
 *
 * Note that METIS does not have a version that uses UF_long integers.  If you
 * try to use cholmod_nested_dissection, cholmod_metis, cholmod_bisect, or
 * cholmod_metis_bisector on a matrix that is too large, an error code will be
 * returned.  METIS does have an "idxtype", which could be redefined as UF_long,
 * if you wish to edit METIS or use compile-time flags to redefine idxtype.
 */

#ifndef AMESOS_CHOLMOD_PARTITION_H
#define AMESOS_CHOLMOD_PARTITION_H

#include "amesos_cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_nested_dissection */
/* -------------------------------------------------------------------------- */

/* Order A, AA', or A(:,f)*A(:,f)' using CHOLMOD's nested dissection method
 * (METIS's node bisector applied recursively to compute the separator tree
 * and constraint sets, followed by CCOLAMD using the constraints).  Usually
 * finds better orderings than METIS_NodeND, but takes longer.
 */

UF_long amesos_cholmod_nested_dissection	/* returns # of components */
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    int *Perm,		/* size A->nrow, output permutation */
    int *CParent,	/* size A->nrow.  On output, CParent [c] is the parent
			 * of component c, or EMPTY if c is a root, and where
			 * c is in the range 0 to # of components minus 1 */
    int *Cmember,	/* size A->nrow.  Cmember [j] = c if node j of A is
			 * in component c */
    /* --------------- */
    cholmod_common *Common
) ;

UF_long amesos_cholmod_l_nested_dissection (cholmod_sparse *, UF_long *, size_t,
    UF_long *, UF_long *, UF_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_metis */
/* -------------------------------------------------------------------------- */

/* Order A, AA', or A(:,f)*A(:,f)' using METIS_NodeND. */

int amesos_cholmod_metis
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int postorder,	/* if TRUE, follow with etree or coletree postorder */
    /* ---- output --- */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int amesos_cholmod_l_metis (cholmod_sparse *, UF_long *, size_t, int, UF_long *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_ccolamd */
/* -------------------------------------------------------------------------- */

/* Order AA' or A(:,f)*A(:,f)' using CCOLAMD. */

int amesos_cholmod_ccolamd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int *Cmember,	/* size A->nrow.  Cmember [i] = c if row i is in the
			 * constraint set c.  c must be >= 0.  The # of
			 * constraint sets is max (Cmember) + 1.  If Cmember is
			 * NULL, then it is interpretted as Cmember [i] = 0 for
			 * all i */
    /* ---- output --- */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int amesos_cholmod_l_ccolamd (cholmod_sparse *, UF_long *, size_t, UF_long *,
    UF_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_csymamd */
/* -------------------------------------------------------------------------- */

/* Order A using CSYMAMD. */

int amesos_cholmod_csymamd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    /* ---- output --- */
    int *Cmember,	/* size nrow.  see cholmod_ccolamd above */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int amesos_cholmod_l_csymamd (cholmod_sparse *, UF_long *, UF_long *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_camd */
/* -------------------------------------------------------------------------- */

/* Order A using CAMD. */

int amesos_cholmod_camd
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    int *Cmember,	/* size nrow.  see cholmod_ccolamd above */
    int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
) ;

int amesos_cholmod_l_camd (cholmod_sparse *, UF_long *, size_t, UF_long *, UF_long *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_bisect */
/* -------------------------------------------------------------------------- */

/* Finds a node bisector of A, A*A', A(:,f)*A(:,f)'. */

UF_long amesos_cholmod_bisect	/* returns # of nodes in separator */
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to bisect */
    int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int compress,	/* if TRUE, compress the graph first */
    /* ---- output --- */
    int *Partition,	/* size A->nrow.  Node i is in the left graph if
			 * Partition [i] = 0, the right graph if 1, and in the
			 * separator if 2. */
    /* --------------- */
    cholmod_common *Common
) ;

UF_long amesos_cholmod_l_bisect (cholmod_sparse *, UF_long *, size_t, int, UF_long *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_metis_bisector */
/* -------------------------------------------------------------------------- */

/* Find a set of nodes that bisects the graph of A or AA' (direct interface
 * to METIS_NodeComputeSeparator). */

UF_long amesos_cholmod_metis_bisector	/* returns separator size */
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to bisect */
    int *Anw,		/* size A->nrow, node weights */
    int *Aew,		/* size nz, edge weights */
    /* ---- output --- */
    int *Partition,	/* size A->nrow.  see cholmod_bisect above. */
    /* --------------- */
    cholmod_common *Common
) ;

UF_long amesos_cholmod_l_metis_bisector (cholmod_sparse *, UF_long *, UF_long *,
    UF_long *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_collapse_septree */
/* -------------------------------------------------------------------------- */

/* Collapse nodes in a separator tree. */

UF_long amesos_cholmod_collapse_septree
(
    /* ---- input ---- */
    size_t n,		/* # of nodes in the graph */
    size_t ncomponents,	/* # of nodes in the separator tree (must be <= n) */
    double nd_oksep,    /* collapse if #sep >= nd_oksep * #nodes in subtree */
    size_t nd_small,    /* collapse if #nodes in subtree < nd_small */
    /* ---- in/out --- */
    int *CParent,	/* size ncomponents; from cholmod_nested_dissection */
    int *Cmember,	/* size n; from cholmod_nested_dissection */
    /* --------------- */
    cholmod_common *Common
) ;

UF_long amesos_cholmod_l_collapse_septree (size_t, size_t, double, size_t, UF_long *,
    UF_long *, cholmod_common *) ;

#endif
