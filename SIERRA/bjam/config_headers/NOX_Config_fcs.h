
#define HAVE_MPI

/* #undef HAVE_NOX_DEBUG */

/* #undef WITH_PRERELEASE */


#define HAVE_LOCA_ANASAZI

#define HAVE_NOX_AMESOS

#define HAVE_NOX_BELOS

#define HAVE_NOX_EPETRAEXT

#define HAVE_NOX_ML_EPETRA

#define HAVE_NOX_STRATIMIKOS

/* #undef HAVE_NOX_TEKO */

/* #undef HAVE_NOX_THYRA */

/* template qualifier required for calling template methods from non-template
   code */
/* #undef INVALID_TEMPLATE_QUALIFIER */

/* define if STL map key is required to be const */
/* #undef MUST_CONST_STL_MAP_KEY */

/* If optional MF TPL library is available */
/* #undef HAVE_LOCA_MF */

#ifndef NOX_FUNC_TIME_MONITOR
#  define NOX_TEUCHOS_TIME_MONITOR
#  define NOX_FUNC_TIME_MONITOR(FUNCNAME) \
     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, NOX)
#  define NOX_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF) \
     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)
#endif

