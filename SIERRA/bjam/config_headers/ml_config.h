/* src/cmake/ml_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */

/* Define the Fortran name mangling to be used for the BLAS */
#define F77_BLAS_MANGLE(name,NAME) name ## _

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

/* define if bool is a built-in type */
/* #undef HAVE_BOOL */

/* Define if the C complier supports __PRETTY_FUNCTION__ */
/* #undef HAVE_CFUNC */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define if want to build libcheck */
/* #undef HAVE_LIBCHECK */

/* Define if want to build with ml enabled */
#define HAVE_ML_AMESOS

/* Define if want to build with ml enabled */
#define HAVE_ML_AZTECOO

/* Define if want to build with ml_benchmarking enabled */
/* #undef HAVE_ML_BENCHMARKING */

/* Define if want to build with ml_complex_maxwell enabled */
/* #undef HAVE_ML_COMPLEX_MAXWELL */

/* Define if want to build with ml_enrich enabled */
/* #undef HAVE_ML_ENRICH */

/* Define if want to build with ml enabled */
#define HAVE_ML_EPETRA

/* Define if want to build with ml enabled */
#define HAVE_ML_EPETRAEXT

/* Define if want to build with ml_flops enabled */
/* #undef HAVE_ML_FLOPS */

/* Define if want to build with ml enabled */
#define HAVE_ML_GALERI

/* Define if want to build with ml enabled */
#define HAVE_ML_IFPACK

/* Define if want to build with ml enabled */
/* #undef HAVE_ML_ISORROPIA */

/* Define if want to build ml-matlab */
/* #undef HAVE_ML_MATLAB */

/* Define if want to build with ml_memory_checking enabled */
/* #undef HAVE_ML_MEMORY_CHECK */

/* Define if want to build with ml_metis enabled */
/* #undef HAVE_ML_METIS */

/* Define if want to build with ml_newtpe enabled */
/* #undef HAVE_ML_NEW_T_PE */

/* Define if want to build with ml_parasails enabled */
/* #undef HAVE_ML_PARASAILS */

/* Define if want to build with ml_parmetis enabled */
/* #undef HAVE_ML_PARMETIS */

/* Define if want to build with ml_parpack enabled */
/* #undef HAVE_ML_PARPACK */

/* Define if want to build with ml_superlu enabled */
/* #undef HAVE_ML_SUPERLU */

/* Define if want to build with ml_superlu2 enabled */
/* #undef HAVE_ML_SUPERLU2_0 */

/* Define if want to build with ml_superlu4 enabled */
/* #undef HAVE_ML_SUPERLU4_0 */

/* Define if want to build with ml_superlu_dist enabled */
/* #undef HAVE_ML_SUPERLU_DIST */

/* Define if want to build with ml enabled */
#define HAVE_ML_TEUCHOS

/* Define if want to build with ml_timing enabled */
/* #undef HAVE_ML_TIMING */

/* Define if want to build with ml_zoltan enabled */
#define HAVE_ML_ZOLTAN

/* Define if want to build with ml_zoltan3 enabled */
/* #undef HAVE_ML_ZOLTAN_THREE */

/* define if we want to use MPI */
#define HAVE_MPI

/* define if the compiler accepts the new for scoping rules */
#define HAVE_NEW_FOR_SCOPING

/* Define if want to build with petsc enabled */
/* #undef HAVE_PETSC */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Support for 64-bit integers */
/* #undef ML_BIG_INT */

/* mallinfo with ML */
/* #undef ML_MALLINFO */

/* Support for multiple right hand sides. */
/* #undef ML_MULTIPLE_RHS_BLOCK_FACTOR */

/* Define to 1 if you have the ANSI C header files. */
/* #undef STDC_HEADERS */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Turn on Teko smoothers, this is a circular dependency use at your own risk */
/* #undef HAVE_ML_TekoSmoothers */
