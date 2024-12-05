/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef _AZ_F77FUNC_H_
#define _AZ_F77FUNC_H_

#if defined(AztecOO_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The AztecOO package is deprecated"
#endif
#endif

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
#define PREFIX

#ifdef TRILINOS_NO_CONFIG_H

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) lcase
#define F77_FUNC_(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) lcase ## _
#if defined(SGI32) || defined(SGI64)
#define F77_FUNC_(lcase,UCASE) lcase ## _
#else
#define F77_FUNC_(lcase,UCASE) lcase ## __
#endif
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */
#endif /* TRILINOS_NO_CONFIG_H */

#endif

#endif /* _AZ_F77FUNC_H_ */
