
#ifndef _FEI_config_h_
#define _FEI_config_h_


/* Enable use of FETI library */
/* #undef HAVE_FEI_FETI */

/* Define to build with fei-amesos support enabled */
#define HAVE_FEI_AMESOS

/* Define to build with fei-aztecoo support enabled */
#define HAVE_FEI_AZTECOO

/* Define to build with fei-aztecoo support enabled */
#define HAVE_FEI_BELOS

/* Define to build with fei-ml support enabled */
#define HAVE_FEI_ML

/* Define to build with fei-ifpack support enabled */
#define HAVE_FEI_IFPACK

/* Define to build with fei-epetra support enabled */
#define HAVE_FEI_EPETRA

/* Define to build with fei-teuchos support enabled */
/* #undef HAVE_FEI_TEUCHOS */

/* Define if want to build examples */
/* #undef HAVE_EXAMPLES */

/* Define if want to build tests */
/* #undef HAVE_TESTS */

/* define if we want to use MPI */
#define HAVE_MPI

/* define if we want to use Boost  */
/* #undef HAVE_FEI_BOOST */

#define FEI_DEPRECATED
#define FEI_DEPRECATED_MSG(MSG)


#endif

