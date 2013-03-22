/* src/include/Zoltan_config.h.in.  Generated from configure.ac by autoheader.  */

/* KDD Copied F77 macros from packages/epetra/cmake/Epetra_config.h.in. */
/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */


/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* ZOLTAN_ID_TYPE is unsigned int */
/* #undef UNSIGNED_INT_GLOBAL_IDS */

/* ZOLTAN_ID_TYPE is unsigned long */
/* #undef UNSIGNED_LONG_GLOBAL_IDS */

/* ZOLTAN_ID_TYPE is unsigned long long */
/* #undef UNSIGNED_LONG_LONG_GLOBAL_IDS */

/* define if we want to use MPI */
#define HAVE_MPI

/* Define if want to build with nemesis_exodus enabled */
/* #undef HAVE_NEMESIS_EXODUS */

/* Define if want to build with parmetis enabled */
/* #undef HAVE_METIS */

/* Define if want to build with parmetis enabled */
#define HAVE_PARMETIS

/* Define if want to build with patoh enabled */
/* #undef HAVE_PATOH */

/* Define if want to build with scotch enabled */
/* #undef HAVE_SCOTCH */

/* Define if want to build with OVIS enabled */
/* #undef HAVE_OVIS */

/* Define if DON'T want support for MPI TPL */
#ifndef HAVE_MPI
#define NO_MPI_TPL
#endif

/* Define if want to build with zlib enabled */
/* #undef HAVE_GZIP */

/* Use to have only filename when debugging memory */
#define SHORT_FILE

/* HUND support */
/* #undef HAVE_ZOLTAN_HUND */

/* Revert to Old Hash function support */
/* #undef HAVE_ZOLTAN_KNUTH_HASH */

#ifdef HAVE_ZOLTAN_HUND
#define CEDRIC_2D_PARTITIONS
#endif
