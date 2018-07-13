/* ========================================================================== */
/* === UFconfig.h =========================================================== */
/* ========================================================================== */

/* Configuration file for SuiteSparse: a Suite of Sparse matrix packages
 * (AMD, TRILINOS_COLAMD, CCOLAMD, CAMD, CHOLMOD, UMFPACK, CXSparse, and others).
 *
 * UFconfig.h provides the definition of the long integer.  On most systems,
 * a C program can be compiled in LP64 mode, in which long's and pointers are
 * both 64-bits, and int's are 32-bits.  Windows 64, however, uses the LLP64
 * model, in which int's and long's are 32-bits, and long long's and pointers
 * are 64-bits.
 *
 * SuiteSparse packages that include long integer versions are
 * intended for the LP64 mode.  However, as a workaround for Windows 64
 * (and perhaps other systems), the long integer can be redefined.
 *
 * If _WIN64 is defined, then the __int64 type is used instead of long.
 *
 * The long integer can also be defined at compile time.  For example, this
 * could be added to UFconfig.mk:
 *
 * CFLAGS = -O -D'UF_long=long long' -D'UF_long_max=9223372036854775801' \
 *   -D'UF_long_id="%lld"'
 *
 * This file defines UF_long as either long (on all but _WIN64) or
 * __int64 on Windows 64.  The intent is that a UF_long is always a 64-bit
 * integer in a 64-bit code.  ptrdiff_t might be a better choice than long;
 * it is always the same size as a pointer.
 *
 * This file also defines the SUITESPARSE_VERSION and related definitions.
 *
 * Copyright (c) 2007, University of Florida.  No licensing restrictions
 * apply to this file or to the UFconfig directory.  Author: Timothy A. Davis.
 */

#ifndef _UFCONFIG_H
#define _UFCONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <limits.h>

/* ========================================================================== */
/* === UF_long ============================================================== */
/* ========================================================================== */

#ifndef UF_long

#ifdef _WIN64

#define UF_long __int64
#define UF_long_max _I64_MAX
#define UF_long_id "%I64d"

#else

#define UF_long long
#define UF_long_max LONG_MAX
#define UF_long_id "%ld"

#endif
#endif

/* ========================================================================== */
/* === SuiteSparse version ================================================== */
/* ========================================================================== */

/* SuiteSparse is not a package itself, but a collection of packages, some of
 * which must be used together (UMFPACK requires AMD, CHOLMOD requires AMD,
 * TRILINOS_COLAMD, CAMD, and CCOLAMD, etc).  A version number is provided here for the
 * collection itself.  The versions of packages within each version of
 * SuiteSparse are meant to work together.  Combining one packge from one
 * version of SuiteSparse, with another package from another version of
 * SuiteSparse, may or may not work.
 *
 * SuiteSparse Version 3.0.0 contains the following packages:
 *
 *  AMD		    version 2.2.0
 *  CAMD	    version 2.2.0
 *  TRILINOS_COLAMD	    version 2.7.0
 *  CCOLAMD	    version 2.7.0
 *  CHOLMOD	    version 1.5.0
 *  CSparse	    version 2.2.0
 *  CXSparse	    version 2.2.0
 *  KLU		    version 1.0.0
 *  TRILINOS_BTF		    version 1.0.0
 *  LDL		    version 2.0.0
 *  UFconfig	    version number is the same as SuiteSparse
 *  UMFPACK	    version 5.1.0
 *  RBio	    version 1.1.0
 *  UFcollection    version 1.1.0
 *
 * Other package dependencies:
 *  BLAS	    required by CHOLMOD and UMFPACK
 *  LAPACK	    required by CHOLMOD
 *  METIS 4.0.1	    required by CHOLMOD (optional) and KLU (optional)
 */

#define TRILINOS_SUITESPARSE_DATE "May 31, 2007"
#define TRILINOS_SUITESPARSE_VER_CODE(main,sub) ((main) * 1000 + (sub))
#define TRILINOS_SUITESPARSE_MAIN_VERSION 3
#define TRILINOS_SUITESPARSE_SUB_VERSION 0
#define TRILINOS_SUITESPARSE_SUBSUB_VERSION 0
#define TRILINOS_SUITESPARSE_VERSION \
    TRILINOS_SUITESPARSE_VER_CODE(TRILINOS_SUITESPARSE_MAIN_VERSION,TRILINOS_SUITESPARSE_SUB_VERSION)

#ifdef __cplusplus
}
#endif
#endif
