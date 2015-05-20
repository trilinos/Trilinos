/* Define if want to build teuchos-complex */
#define HAVE_SACADO_COMPLEX

/* Define if want to build with uninit */
/* #undef HAVE_SACADO_UNINIT */

/* Define if ADOL-C is enabled */
/* #undef HAVE_ADOLC */

/* Define if ADIC is enabled */
/* #undef HAVE_ADIC */

/* define if we want to use MPI */
#define HAVE_MPI

/* Define if want to build sacado-examples */
/* #undef HAVE_SACADO_EXAMPLES */

/* Define if want to build sacado-tests */
/* #undef HAVE_SACADO_TESTS */

/* Define if want to build with teuchos enabled */
#define HAVE_SACADO_TEUCHOS

/* Define if want to build with kokkos-core enabled */
#define HAVE_SACADO_KOKKOSCORE

/* Define if want to build with TeuchosKokkosComm subpackage enabled */
#define HAVE_SACADO_TEUCHOSKOKKOSCOMM
#ifdef HAVE_SACADO_TEUCHOSKOKKOSCOMM
// For backwards compatibility
#  define HAVE_SACADO_KOKKOSMPICOMM
#endif // HAVE_SACADO_TEUCHOSKOKKOSCOMM

/* Define if want to enable Kokkos view specializations for Sacado */
#define HAVE_SACADO_VIEW_SPEC

/* define if the compiler is confused by std::sin, etc., within namespace
   Sacado::Rad */
/* #undef RAD_NO_USING_STDCC */

/* Define to enable extra debugging checks */
/* #undef SACADO_DEBUG */

/* Define if compiler supports c99 tr1 cmath functions */
#define HAS_C99_TR1_CMATH

/* Define to enable C++11 support*/
#define HAVE_SACADO_CXX11
