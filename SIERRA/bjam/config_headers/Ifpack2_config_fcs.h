/* Ifpack2_config.h.in is preprocessed by CMake to create Ifpack2_config.h. */

/*
 * Build options
 */

/* Define if we want to build Ifpack2 with debug support */
/* #undef HAVE_IFPACK2_DEBUG */

/* Define if explicit (template) instantiation is enabled. */
/* #undef HAVE_IFPACK2_EXPLICIT_INSTANTIATION */

/* Define if experimental Ifpack2 code is enabled. */
/* #undef HAVE_IFPACK2_EXPERIMENTAL */

/* Define if Ifpack2::SupportGraph is enabled. */
/* #undef HAVE_IFPACK2_SUPPORTGRAPH */

/*
 * Package dependencies
 */

/* Define if we have Amesos2 */
#define HAVE_IFPACK2_AMESOS2

/* Define if we have Belos */
/* #undef HAVE_IFPACK2_BELOS */

/* Define if we have Galeri */
/* #undef HAVE_IFPACK2_GALERI */

/* Define if we have KokkosClassic */
/* #undef HAVE_IFPACK2_KOKKOSCLASSIC */

/* Define if we have ThyraTpetraAdapters */
#define HAVE_IFPACK2_THYRATPETRAADAPTERS

/* Define if we have Xpetra */
#define HAVE_IFPACK2_XPETRA

/* Define if we have Zoltan2 */
#define HAVE_IFPACK2_ZOLTAN2

/*
 * TPL (third-party library) dependencies
 */

/* Define if building with the Cholmod library */
/* #undef HAVE_IFPACK2_CHOLMOD */

/* Define if building with the Lemon library */
/* #undef HAVE_IFPACK2_LEMON */

/* Define if building with MPI (the Message-Passing Interface) */
/* #undef HAVE_IFPACK2_MPI */

/* Define if we have the QD extended-precision TPL */
/* #undef HAVE_IFPACK2_QD */

