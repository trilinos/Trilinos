// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#ifdef HAVE_TEUCHOS_INST_FLOAT
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


#ifdef HAVE_TEUCHOS_INST_COMPLEX_FLOAT
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<float>)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME)\
     template class CLASSNAME<std::complex<float> >;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_FLOAT(CLASSNAME)
#endif


#ifdef HAVE_TEUCHOS_INST_COMPLEX_DOUBLE
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<double>)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)\
     template class CLASSNAME<std::complex<double> >;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_COMPLEX_DOUBLE(CLASSNAME)
#endif

#ifdef HAVE_TEUCHOS_LONG_DOUBLE
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_LONG_DOUBLE(INSTANT_MACRO)\
     INSTANT_MACRO(long double)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(CLASSNAME)\
     template class CLASSNAME<long double>;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_LONG_DOUBLE(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(CLASSNAME)
#endif


#ifdef HAVE_TEUCHOSCORE_QUADMATH
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT128(INSTANT_MACRO)\
      INSTANT_MACRO(__float128)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT128(CLASSNAME)\
     template class CLASSNAME<__float128>;
#else
#  define TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT128(INSTANT_MACRO)
#  define TEUCHOS_CLASS_TEMPLATE_INSTANT_FLOAT128(CLASSNAME)
#endif // HAVE_TEUCHOSCORE_QUADMATH


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
