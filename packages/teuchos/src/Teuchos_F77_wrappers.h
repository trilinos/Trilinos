#ifndef _TEUCHOS_F77_WRAPPERS_H_
#define _TEUCHOS_F77_WRAPPERS_H_

/*! \file Teuchos_F77_wrappers.h  
    \brief Macros for portably calling Fortran77 from C/C++
*/

#include "Teuchos_ConfigDefs.hpp"

/* Define fcd (Fortran Teuchos_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)

#  if defined(CRAY_T3X)

#    include <fortran.h>
#    define F77_CALL_PREFIX
#    define FORTRAN_CHAR_1_ARG(ARG_NAME) fcd ARG_NAME
#    define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const fcd ARG_NAME

#  elif defined(INTEL_CXML)

#    define F77_CALL_PREFIX __stdcall 
#    define FORTRAN_CHAR_1_ARG(ARG_NAME) char* ARG_NAME, unsigned int
#    define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const char* ARG_NAME, unsigned int

#  elif defined(INTEL_MKL)

#    define F77_CALL_PREFIX
#    define FORTRAN_CHAR_1_ARG(ARG_NAME) char* ARG_NAME
#    define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const char* ARG_NAME

#  endif 

#endif

/* External macros */

#define FORTRAN_NAME_UL(UNAME,LNAME) F77_FUNC(UNAME,LNAME)

#define FORTRAN_FUNC_DECL_UL(TYPE,UFUNC_NAME,LFUNC_NAME) TYPE F77_CALL_PREFIX FORTRAN_NAME_UL(UFUNC_NAME,LFUNC_NAME)

#define FORTRAN_FUNC_CALL_UL(UFUNC_NAME,LFUNC_NAME) FORTRAN_NAME_UL(UFUNC_NAME,LFUNC_NAME)

#define FORTRAN_COMMMON_BLOCK_NAME_UL(UNAME,LNAME)  FORTRAN_NAME_UL(UFUNC_NAME,LFUNC_NAME)

#endif // _TEUCHOS_F77_WRAPPERS_H_
