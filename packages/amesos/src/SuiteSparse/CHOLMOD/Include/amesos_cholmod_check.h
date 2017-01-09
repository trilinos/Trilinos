/* ========================================================================== */
/* === Include/cholmod_check.h ============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_check.h.  Copyright (C) 2005-2006, Timothy A. Davis
 * CHOLMOD/Include/cholmod_check.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* CHOLMOD Check module.
 *
 * Routines that check and print the 5 basic data types in CHOLMOD, and 3 kinds
 * of integer vectors (subset, perm, and parent), and read in matrices from a
 * file:
 *
 * cholmod_check_common	    check/print the Common object
 * cholmod_print_common
 *
 * cholmod_check_sparse	    check/print a sparse matrix in column-oriented form
 * cholmod_print_sparse
 *
 * cholmod_check_dense	    check/print a dense matrix
 * cholmod_print_dense
 *
 * cholmod_check_factor	    check/print a Cholesky factorization
 * cholmod_print_factor
 *
 * cholmod_check_triplet    check/print a sparse matrix in triplet form
 * cholmod_print_triplet
 *
 * cholmod_check_subset	    check/print a subset (integer vector in given range)
 * cholmod_print_subset
 *
 * cholmod_check_perm	    check/print a permutation (an integer vector)
 * cholmod_print_perm
 *
 * cholmod_check_parent	    check/print an elimination tree (an integer vector)
 * cholmod_print_parent
 *
 * cholmod_print_common and cholmod_check_common are the only two routines that
 * you may call after calling cholmod_finish.
 *
 * Requires the Core module.  Not required by any CHOLMOD module, except when
 * debugging is enabled (in which case all modules require the Check module).
 *
 */

#ifndef AMESOS_CHOLMOD_CHECK_H
#define AMESOS_CHOLMOD_CHECK_H

#include "amesos_cholmod_core.h"
#include <stdio.h>

/* -------------------------------------------------------------------------- */
/* cholmod_check_common:  check the Common object */
/* -------------------------------------------------------------------------- */

int cholmod_check_common
(
    cholmod_common *Common
) ;

int cholmod_l_check_common (cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_common:  print the Common object */
/* -------------------------------------------------------------------------- */

int cholmod_print_common
(
    /* ---- input ---- */
    char *name,		/* printed name of Common object */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_common (char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_sparse:  check a sparse matrix */
/* -------------------------------------------------------------------------- */

int cholmod_check_sparse
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_sparse (cholmod_sparse *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_sparse */
/* -------------------------------------------------------------------------- */

int cholmod_print_sparse
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* sparse matrix to print */
    char *name,		/* printed name of sparse matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_sparse (cholmod_sparse *, char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_dense:  check a dense matrix */
/* -------------------------------------------------------------------------- */

int cholmod_check_dense
(
    /* ---- input ---- */
    cholmod_dense *X,	/* dense matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_dense (cholmod_dense *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_dense:  print a dense matrix */
/* -------------------------------------------------------------------------- */

int cholmod_print_dense
(
    /* ---- input ---- */
    cholmod_dense *X,	/* dense matrix to print */
    char *name,		/* printed name of dense matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_dense (cholmod_dense *, char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_factor:  check a factor */
/* -------------------------------------------------------------------------- */

int cholmod_check_factor
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_factor (cholmod_factor *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_factor:  print a factor */
/* -------------------------------------------------------------------------- */

int cholmod_print_factor
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factor to print */
    char *name,		/* printed name of factor */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_factor (cholmod_factor *, char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_triplet:  check a sparse matrix in triplet form */
/* -------------------------------------------------------------------------- */

int cholmod_check_triplet
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_triplet (cholmod_triplet *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_triplet:  print a triplet matrix */
/* -------------------------------------------------------------------------- */

int cholmod_print_triplet
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to print */
    char *name,		/* printed name of triplet matrix */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_triplet (cholmod_triplet *, char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_subset:  check a subset */
/* -------------------------------------------------------------------------- */

int cholmod_check_subset
(
    /* ---- input ---- */
    int *Set,		/* Set [0:len-1] is a subset of 0:n-1.  Duplicates OK */
    UF_long len,	/* size of Set (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_subset (UF_long *, UF_long, size_t, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_subset:  print a subset */
/* -------------------------------------------------------------------------- */

int cholmod_print_subset
(
    /* ---- input ---- */
    int *Set,		/* Set [0:len-1] is a subset of 0:n-1.  Duplicates OK */
    UF_long len,	/* size of Set (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    char *name,		/* printed name of Set */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_subset (UF_long *, UF_long, size_t, char *,
    cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_perm:  check a permutation */
/* -------------------------------------------------------------------------- */

int cholmod_check_perm
(
    /* ---- input ---- */
    int *Perm,		/* Perm [0:len-1] is a permutation of subset of 0:n-1 */
    size_t len,		/* size of Perm (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_perm (UF_long *, size_t, size_t, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_perm:  print a permutation vector */
/* -------------------------------------------------------------------------- */

int cholmod_print_perm
(
    /* ---- input ---- */
    int *Perm,		/* Perm [0:len-1] is a permutation of subset of 0:n-1 */
    size_t len,		/* size of Perm (an integer array) */
    size_t n,		/* 0:n-1 is valid range */
    char *name,		/* printed name of Perm */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_perm (UF_long *, size_t, size_t, char *, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_parent:  check an elimination tree */
/* -------------------------------------------------------------------------- */

int cholmod_check_parent
(
    /* ---- input ---- */
    int *Parent,	/* Parent [0:n-1] is an elimination tree */
    size_t n,		/* size of Parent */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_check_parent (UF_long *, size_t, cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_parent */
/* -------------------------------------------------------------------------- */

int cholmod_print_parent
(
    /* ---- input ---- */
    int *Parent,	/* Parent [0:n-1] is an elimination tree */
    size_t n,		/* size of Parent */
    char *name,		/* printed name of Parent */
    /* --------------- */
    cholmod_common *Common
) ;

int cholmod_l_print_parent (UF_long *, size_t, char *, cholmod_common *) ;

#endif
