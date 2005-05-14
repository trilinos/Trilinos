/* ========================================================================== */
/* === Include/cholmod_partition.h ========================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Partition version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD Partition module.
 *
 * Graph partitioning and graph-partition-based orderings.  Includes an
 * interface to CCOLAMD and CSYMAMD, constrained minimum degree ordering
 * methods which order a matrix following constraints determined via nested
 * dissection.
 *
 * Primary routines:
 * -----------------
 * cholmod_nested_dissection	CHOLMOD nested dissection ordering
 * cholmod_metis		METIS nested dissection ordering (METIS_NodeND)
 *
 * Secondary routines:
 * -------------------
 * cholmod_ccolamd		interface to CCOLAMD ordering
 * cholmod_csymamd		interface to CSYMAMD ordering
 * cholmod_bisect		graph partitioner (currently based on METIS)
 * cholmod_metis_bisector	direct interface to METIS_NodeComputeSeparator
 *
 * Requires the Core module, and two packages: METIS and CCOLAMD.
 * Optionally used by the Cholesky module.
 */

#ifndef CHOLMOD_PARTITION_H
#define CHOLMOD_PARTITION_H

#include "cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_nested_dissection */
/* -------------------------------------------------------------------------- */

long cholmod_nested_dissection
(
    cholmod_sparse *A,	/* order A if symmetric, A*A' if unsymmetric */
    void *fset,
    size_t fsize,

    /* outputs, contents undefined on input */
    void *Perm,	/* size n.  Perm [k] = j if node j is kth node in the
			 * permuted matrix */
    void *CParent,	/* size n.  On output, CParent [c] is the parent
			 * of componebt c, or EMPTY if c is a root. */
    void *Cmember,	/* size n.  Cmember [j] = c if node j of A is
			 * in component c */

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_metis */
/* -------------------------------------------------------------------------- */

int cholmod_metis
(
    cholmod_sparse *A,
    void *fset,
    size_t fsize,

    /* outputs, contents undefined on input */
    void *Perm,		/* size n.  Perm [k] = j if node j is kth node in the
			 * permuted matrix */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_ccolamd */
/* -------------------------------------------------------------------------- */

int cholmod_ccolamd
(
    /* inputs, not modified */
    cholmod_sparse *A,	    /* user's input matrix, nrow-by-ncol */
    void *fset,
    size_t fsize,

    void *Cmember,	    /* size nrow.  Cmember [i] = c if row i is in the
			     * constraint set c.  c must be >= 0.  The # of
			     * constraint sets is max (Cmember) + 1.  If
			     * Cmember is NULL, then it is interpretted as
			     * Cmember [i] = 0 for all i */

    /* output, not defined on input */
    void *Perm,	    /* size nrow, output permutation */

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_csymamd */
/* -------------------------------------------------------------------------- */

int cholmod_csymamd
(
    cholmod_sparse *A,	    /* user's input matrix, nrow-by-ncol */
    void *Cmember,	    /* size nrow.  see cholmod_ccolamd above */
    void *Perm,	    /* size nrow, output permutation */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_bisect */
/* -------------------------------------------------------------------------- */

long cholmod_bisect
(
    /* input only, not modified on output: */
    cholmod_sparse *A,		/* bisect A if symmetric, A*A' if unsymmetric */
    void *fset,
    size_t fsize,
    int compress,

    /* output, contents undefined on input */
    void *Partition,		/* size A->nrow */

    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_metis_bisector */
/* -------------------------------------------------------------------------- */

long cholmod_metis_bisector
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    void *Anw,		    /* size n, node weights */
    void *Aew,		    /* size nz, edge weights */

    /* output, undefined on input */
    void *Partition,	    /* size n */

    cholmod_common *Common
) ;

#endif
