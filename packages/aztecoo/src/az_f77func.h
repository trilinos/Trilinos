#ifndef _AZ_F77FUNC_H_
#define _AZ_F77FUNC_H_

#include "az_aztec_defs.h"
#include <string.h>

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, 1
#else
#define CHAR_MACRO(char_var) &char_var
#endif

/* Define fcd (Fortran az_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define az_fcd fcd

#elif defined(INTEL_CXML)

#define PREFIX __stdcall
#define az_fcd char *, unsigned int

#elif defined(INTEL_MKL)

#define PREFIX
#define az_fcd char *

#endif 
/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif


#define F77_FUNC(lcase,UCASE) UCASE

#ifdef F77_FUNC_
#undef F77_FUNC_
#endif

#define F77_FUNC_(lcase,UCASE) UCASE

#else /* Define az_fcd for all other machines */

#define az_fcd char * 

#ifndef HAVE_CONFIG_H

#define PREFIX
#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) lcase
#define F77_FUNC_(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) lcase ## _
#define F77_FUNC_(lcase,UCASE) lcase ## __
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */
#endif /* HAVE_CONFIG_H */

#endif

#endif /* _AZ_F77FUNC_H_ */
