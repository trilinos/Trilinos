/* ========================================================================== */
/* === ccolamd_global.c ===================================================== */
/* ========================================================================== */

/* ----------------------------------------------------------------------------
 * CCOLAMD Copyright (C), Univ. of Florida.  Authors: Timothy A. Davis,
 * Sivasankaran Rajamanickam, and Stefan Larimore
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* Global variables for CCOLAMD */

/* For TRILINOSSS_LIB_DLL_EXPORT: dllexport'd from trilinosss.dll, dllimport'd by
 * consumers on Windows shared builds. */
#include "trilinos_UFconfig.h"

#ifndef NPRINT
#ifdef MATLAB_MEX_FILE
#include "mex.h"
TRILINOSSS_LIB_DLL_EXPORT int (*trilinos_ccolamd_printf) (const char *, ...) = mexPrintf ;
#else
#include <stdio.h>
TRILINOSSS_LIB_DLL_EXPORT int (*trilinos_ccolamd_printf) (const char *, ...) = printf ;
#endif
#else
TRILINOSSS_LIB_DLL_EXPORT int (*trilinos_ccolamd_printf) (const char *, ...) = ((void *) 0) ;
#endif

