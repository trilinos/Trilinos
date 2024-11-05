// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_F77_WRAPPERS_H_
#define _TEUCHOS_F77_WRAPPERS_H_

/*! \file Teuchos_F77_wrappers.h  
    \brief Macros for portably calling Fortran77 from C/C++
*/

#include "Teuchos_ConfigDefs.hpp"

/* Define fcd (Fortran Teuchos_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X)

#  include <fortran.h>
#  define F77_CALL_PREFIX
#  define FORTRAN_CHAR_1_ARG(ARG_NAME) fcd* ARG_NAME
#  define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const fcd& ARG_NAME
#  define FORTRAN_CHAR_1_ARG_CALL(ARG_NAME) ARG_NAME

#elif defined(INTEL_CXML)

#  define F77_CALL_PREFIX __stdcall 
#  define FORTRAN_CHAR_1_ARG(ARG_NAME) char* ARG_NAME, unsigned int
#  define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const char& ARG_NAME, unsigned int
#  define FORTRAN_CHAR_1_ARG_CALL(ARG_NAME) ARG_NAME, 1

#elif defined(INTEL_MKL)

#  define F77_CALL_PREFIX
#  define FORTRAN_CHAR_1_ARG(ARG_NAME) char* ARG_NAME
#  define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const char& ARG_NAME
#  define FORTRAN_CHAR_1_ARG_CALL(ARG_NAME) ARG_NAME, 1

#else

#  define F77_CALL_PREFIX
#  define FORTRAN_CHAR_1_ARG(ARG_NAME) char* ARG_NAME
#  define FORTRAN_CONST_CHAR_1_ARG(ARG_NAME) const char& ARG_NAME
#  define FORTRAN_CHAR_1_ARG_CALL(ARG_NAME) ARG_NAME

#endif

/* RAB: 20030924: ToDo: Find a way to see if const is supported or not by C or
   just decide that this will only be for C++ code and be done with it. */

/* External macros */

#define FORTRAN_NAME_UL(UNAME,LNAME) F77_FUNC(LNAME,UNAME)

#define FORTRAN_FUNC_DECL_UL(TYPE,UFUNC_NAME,LFUNC_NAME) TYPE F77_CALL_PREFIX FORTRAN_NAME_UL(UFUNC_NAME,LFUNC_NAME)

#define FORTRAN_FUNC_CALL_UL(UFUNC_NAME,LFUNC_NAME) FORTRAN_NAME_UL(UFUNC_NAME,LFUNC_NAME)

#define FORTRAN_FUNC_PTR_DECL_UL(TYPE,UFUNC_NAME,LFUNC_NAME) TYPE (F77_CALL_PREFIX *FORTRAN_NAME_UL(UFUNC_NAME,LFUNC_NAME))

#define FORTRAN_COMMMON_BLOCK_NAME_UL(UNAME,LNAME)  FORTRAN_NAME_UL(UNAME,LNAME)\

#define FORTRAN_NAME_UL_(UNAME,LNAME) F77_FUNC_(LNAME,UNAME)

#define FORTRAN_FUNC_DECL_UL_(TYPE,UFUNC_NAME,LFUNC_NAME) TYPE F77_CALL_PREFIX FORTRAN_NAME_UL_(UFUNC_NAME,LFUNC_NAME)

#define FORTRAN_FUNC_CALL_UL_(UFUNC_NAME,LFUNC_NAME) FORTRAN_NAME_UL_(UFUNC_NAME,LFUNC_NAME)

#define FORTRAN_FUNC_PTR_DECL_UL_(TYPE,UFUNC_NAME,LFUNC_NAME) TYPE (F77_CALL_PREFIX *FORTRAN_NAME_UL_(UFUNC_NAME,LFUNC_NAME))

#define FORTRAN_COMMMON_BLOCK_NAME_UL_(UNAME,LNAME)  FORTRAN_NAME_UL_(UNAME,LNAME)\

#ifdef __cplusplus

// These are the platform dependent C++ equivalents of fortran types
// RAB: 2003/11/20: ToDo: Move this into Teuchos namespace at some point
namespace FortranTypes {

typedef int							f_int;					// INTEGER
typedef float						f_real;					// REAL
typedef double						f_dbl_prec;				// DOUBLE PRECISION
typedef int							f_logical;				// LOGICAL
typedef char						f_char;					// CHARACTER*1
typedef unsigned int				f_char_len;				// length argument for a CHARACTER*(*)
//typedef std::complex<f_real>		f_complex;				// COMPLEX
//typedef std::complex<f_dbl_prec>	f_complex_16;			// COMPLEX*16

enum {	F_TRUE = true, F_FALSE = false }; // Let compiler figure this out!

#endif /* __cplusplus */

} // namespace FortranTypes

#endif // _TEUCHOS_F77_WRAPPERS_H_
