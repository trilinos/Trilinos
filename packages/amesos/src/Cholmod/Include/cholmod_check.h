/* ========================================================================== */
/* === Include/cholmod_check.h ============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Check version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD Check module.
 *
 * Routines that check and print the 5 basic data types in CHOLMOD, and 3 kinds
 * of integer vectors (subset, perm, and parent):
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
 * cholmod_check_parent	    check/print an elimination tree (an int. vector)
 * cholmod_print_parent
 * 
 *
 * cholmod_print_common and cholmod_check_common are the only two routines that
 * you may call after calling cholmod_finish.
 *
 * Requires the Core module.  Not required by any CHOLMOD module, except when
 * debugging is enabled (in which case all modules require the Check module).
 */

#ifndef CHOLMOD_CHECK_H
#define CHOLMOD_CHECK_H

#include "cholmod_core.h"

/* -------------------------------------------------------------------------- */
/* cholmod_check_common */
/* -------------------------------------------------------------------------- */

int cholmod_check_common
(
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_common */
/* -------------------------------------------------------------------------- */

int cholmod_print_common
(
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_sparse */
/* -------------------------------------------------------------------------- */

int cholmod_check_sparse
(
    cholmod_sparse *A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_sparse */
/* -------------------------------------------------------------------------- */

int cholmod_print_sparse
(
    cholmod_sparse *A,
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_dense */
/* -------------------------------------------------------------------------- */

int cholmod_check_dense
(
    cholmod_dense *X,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_dense */
/* -------------------------------------------------------------------------- */

int cholmod_print_dense
(
    cholmod_dense *X,
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_factor */
/* -------------------------------------------------------------------------- */

int cholmod_check_factor
(
    cholmod_factor *L,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_factor */
/* -------------------------------------------------------------------------- */

int cholmod_print_factor
(
    cholmod_factor *L,
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_triplet */
/* -------------------------------------------------------------------------- */

int cholmod_check_triplet
(
    cholmod_triplet *T,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_triplet */
/* -------------------------------------------------------------------------- */

int cholmod_print_triplet
(
    cholmod_triplet *T,
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_subset */
/* -------------------------------------------------------------------------- */

int cholmod_check_subset
(
    void *Set,
    size_t len,
    size_t n,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_subset */
/* -------------------------------------------------------------------------- */

int cholmod_print_subset
(
    void *Set,
    size_t len,
    size_t n,
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_perm */
/* -------------------------------------------------------------------------- */

int cholmod_check_perm
(
    void *Perm,
    size_t len,
    size_t n,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_perm */
/* -------------------------------------------------------------------------- */

int cholmod_print_perm
(
    void *Perm,
    size_t len,
    size_t n,
    char *name,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_parent */
/* -------------------------------------------------------------------------- */

int cholmod_check_parent
(
    void *Parent,
    size_t n,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_parent */
/* -------------------------------------------------------------------------- */

int cholmod_print_parent
(
    void *Parent,
    size_t n,
    char *name,
    cholmod_common *Common
) ;

#endif
