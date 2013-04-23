/* src/EpetraExt_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if want to build epetraext-block */
#define HAVE_BLOCK

/* Define if want to build epetraext-btf */
/* #undef HAVE_BTF */

/* Define if want to build epetraext-coloring */
#define HAVE_COLORING

/* Define if want to build epetraext-hdf5 */
/* #undef HAVE_EPETRAEXT_HDF5 */

/* Define if want to build with epetraext enabled */
#define HAVE_EPETRAEXT_TRIUTILS

/* Define if want to build epetraext-experimental */
/* #undef HAVE_EXPERIMENTAL */

/* Define if want to build with fatal_messages enabled */
/* #undef HAVE_FATAL_MESSAGES */

/* Define if want to build with format_io enabled */
/* #undef HAVE_FORMAT_IO */

/* Define if want to build epetraext-graphreorderings */
/* #undef HAVE_GRAPH_REORDERINGS */

/* Define if want to build epetraext-inout */
#define HAVE_INOUT

/* Define to 1 if you have the `eng' library (-leng). */
/* #undef HAVE_LIBENG */

/* Define to 1 if you have the `metis' library (-lmetis). */
/* #undef HAVE_LIBMETIS */

/* Define to 1 if you have the `mx' library (-lmx). */
/* #undef HAVE_LIBMX */

/* Define to 1 if you have the `parmetis' library (-lparmetis). */
#define HAVE_LIBPARMETIS

/* Define if want to build epetraext-matlab */
/* #undef HAVE_MATLAB */

/* Define if want to build epetraext-modelevaluator */
#define HAVE_MODEL_EVALUATOR

/* Define if want to build with petsc enabled */
/* #undef HAVE_PETSC */

/* Define if want to build with hypre enabled */
/* #undef HAVE_HYPRE */

/* Define if want to build epetraext-restrict */
#define HAVE_RESTRICT

/* Define if want to build with threads enabled */
/* #undef HAVE_THREADS */

/* Define if want to build epetraext-transform */
#define HAVE_TRANSFORM

/* Define if want to build with warning_messages enabled */
/* #undef HAVE_WARNING_MESSAGES */

/* Define if PyTrilinos is enabled */
/* #undef HAVE_PYTRILINOS */

/* define if STL map key is required to be const */
/* #undef MUST_CONST_STL_MAP_KEY */

/* define if we want to use MPI */
#define HAVE_MPI

#ifndef EPETRAEXT_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define EPETRAEXT_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define EPETRAEXT_DEPRECATED
#  endif
#endif

