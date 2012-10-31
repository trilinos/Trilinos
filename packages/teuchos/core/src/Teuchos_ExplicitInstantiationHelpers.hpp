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

#ifndef TEUCHOS_EXPLICIT_INSTANTIATION_HELPERS_HPP
#define TEUCHOS_EXPLICIT_INSTANTIATION_HELPERS_HPP


/*! \file Teuchos_ExplicitInstantiationHelpers.hpp

\brief Macros for helping to templated classes create explicit instantiations.

*/


#include "Teuchos_ConfigDefs.hpp"


//
// 2007/07/10: rabartl: NOTE: Semicolons must only be used at the lowest level
// of final code to ensure that there will not be any empty semicolon lines
// that might issue a compiler warning or error. In general, I like to define
// macros that need a semicolon when you use them because my emacs mode will
// then do the indentation correctly.  However, this is not a big deal since
// these macros only get used in a final *.cpp file and at that point they are
// only used once in the entire mostly empty file.
//

#ifdef HAVE_TEUCHOS_FLOAT
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT(INSTANT_MACRO)\
    INSTANT_MACRO(float)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT(CLASSNAME)\
    template class CLASSNAME<float>;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT(CLASSNAME)
#endif


#define TEUCHOS_MACRO_TEMPLATE_INSTANT_DOUBLE(INSTANT_MACRO)\
  INSTANT_MACRO(double)
#define TEUCHOS_CLASS_TEMPLATE_INSTANT_DOUBLE(CLASSNAME)\
  template class CLASSNAME<double>;

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_FLOAT)
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<float>)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME)\
     template class CLASSNAME<std::complex<float> >;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME)
#endif


#ifdef HAVE_TEUCHOS_COMPLEX
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<double>)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)\
     template class CLASSNAME<std::complex<double> >;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)
#endif


/** \brief Instantiate a macro template for the set of supported real scalar
 * types.
 */
#define TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(MACRONAME) \
  TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT(MACRONAME) \
  TEUCHOS_MACRO_TEMPLATE_INSTANT_DOUBLE(MACRONAME)


/** \brief Instantiate a macro template for the set of supported real and
 * complex scalar types.
 */
#define TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(MACRONAME)\
  TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT(MACRONAME) \
  TEUCHOS_MACRO_TEMPLATE_INSTANT_DOUBLE(MACRONAME) \
  TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT(MACRONAME) \
  TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE(MACRONAME)


/** \brief Instantiate a class template for the set of supported real scalar
 * types.
 */
#define TEUCHOS_CLASS_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(CLASSNAME)\
  TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT(CLASSNAME) \
  TEUCHOS_CLASS_TEMPLATE_INSTANT_DOUBLE(CLASSNAME)


/** \brief Instantiate a class template for the set of supported real and
 * complex scalar types.
 */
#define TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(CLASSNAME)\
  TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT(CLASSNAME) \
  TEUCHOS_CLASS_TEMPLATE_INSTANT_DOUBLE(CLASSNAME) \
  TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME) \
  TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)


#endif  // TEUCHOS_EXPLICIT_INSTANTIATION_HELPERS_HPP
