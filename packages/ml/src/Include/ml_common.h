/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Include file that should be included in every header file.           */
/* ******************************************************************** */
/* Author        : Jonathan Hu                                          */
/* Date          : July, 2003                                           */
/* ******************************************************************** */

#ifndef __MLCOMMON__
#define __MLCOMMON__

#ifdef HAVE_CONFIG_H 
#include "ml_config.h"

#ifdef HAVE_ML_AZTEC
#define AZTEC
#endif

#ifdef HAVE_ML_SUPERLU
#define SUPERLU
#endif

#ifdef HAVE_MPI
#define ML_MPI
#endif

#if defined(HAVE_ML_EXTERNAL_MPI_FUNCTIONS) && defined(HAVE_MPI)
#define ML_USING_MPI_FUNCTIONS
#endif

#ifdef HAVE_BLAS
#define USE_VENDOR_BLAS
#endif

#ifdef HAVE_LAPACK
#define USE_VENDOR_LAPACK
#endif

#ifdef HAVE_ML_ENRICH
#define ML_ENRICH
#endif

#ifdef HAVE_ML_NEW_T_PE
#define ML_NEW_T_PE
#endif

#ifdef HAVE_ML_COMPLEX_MAXWELL
#define GREG
#endif

#ifdef HAVE_ML_TIMING
#define ML_TIMING
#endif


#endif /*ifdef HAVE_CONFIG_H*/

#endif
