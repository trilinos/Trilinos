// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef TEUCHOS_EXPAND_SCALAR_TYPE_MACROS_HPP
#define TEUCHOS_EXPAND_SCALAR_TYPE_MACROS_HPP


/*! \file Teuchos_UnitTestHelpers.hpp

\brief Macros for expanding code
*/


#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_TEUCHOS_FLOAT
#  define TEUCHOS_MACRO_EXPAND_FLOAT(INSTANT_MACRO)\
    INSTANT_MACRO(float)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT(CLASSNAME)\
    template class CLASSNAME<float>;
#else
#  define TEUCHOS_MACRO_EXPAND_FLOAT(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT(CLASSNAME)
#endif


#define TEUCHOS_MACRO_EXPAND_DOUBLE(INSTANT_MACRO)\
  INSTANT_MACRO(double)
#define TEUCHOS_CLASS_TEMPLATE_INSTANT_DOUBLE(CLASSNAME)\
  template class CLASSNAME<double>;

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_FLOAT)
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_FLOAT(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<float>)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME)\
     template class CLASSNAME<std::complex<float> >;
#else
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_FLOAT(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME)
#endif


#ifdef HAVE_TEUCHOS_COMPLEX
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_DOUBLE(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<double>)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)\
     template class CLASSNAME<std::complex<double> >;
#else
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_DOUBLE(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)
#endif


/** \brief Instantiate a macro template for the set of supported real scalar
 * types.
 */
#define TEUCHOS_MACRO_EXPAND_REAL_SCALAR_TYPES(MACRONAME) \
  TEUCHOS_MACRO_EXPAND_FLOAT(MACRONAME) \
  TEUCHOS_MACRO_EXPAND_DOUBLE(MACRONAME)


/** \brief Instantiate a macro template for the set of supported real and
 * complex scalar types.
 */
#define TEUCHOS_MACRO_EXPAND_SCALAR_TYPES(MACRONAME)\
  TEUCHOS_MACRO_EXPAND_FLOAT(MACRONAME) \
  TEUCHOS_MACRO_EXPAND_DOUBLE(MACRONAME) \
  TEUCHOS_MACRO_EXPAND_COMPLEX_FLOAT(MACRONAME) \
  TEUCHOS_MACRO_EXPAND_COMPLEX_DOUBLE(MACRONAME)


#endif  // TEUCHOS_EXPAND_SCALAR_TYPE_MACROS_HPP
