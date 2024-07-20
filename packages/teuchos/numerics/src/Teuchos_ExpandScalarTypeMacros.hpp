// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_EXPAND_SCALAR_TYPE_MACROS_HPP
#define TEUCHOS_EXPAND_SCALAR_TYPE_MACROS_HPP


/*! \file Teuchos_ExpandScaleTypeMacros.hpp

\brief Macros for expanding code
*/


#include "Teuchos_ExplicitInstantiationHelpers.hpp"

#ifdef HAVE_TEUCHOS_INST_FLOAT
#  define TEUCHOS_MACRO_EXPAND_FLOAT(INSTANT_MACRO)\
    INSTANT_MACRO(float)
#else
#  define TEUCHOS_MACRO_EXPAND_FLOAT(INSTANT_MACRO)
#endif


#define TEUCHOS_MACRO_EXPAND_DOUBLE(INSTANT_MACRO)\
  INSTANT_MACRO(double)

#ifdef HAVE_TEUCHOS_INST_COMPLEX_FLOAT
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_FLOAT(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<float>)
#else
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_FLOAT(INSTANT_MACRO)
#endif


#ifdef HAVE_TEUCHOS_INST_COMPLEX_DOUBLE
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_DOUBLE(INSTANT_MACRO)\
     INSTANT_MACRO(std::complex<double>)
#else
#  define TEUCHOS_MACRO_EXPAND_COMPLEX_DOUBLE(INSTANT_MACRO)
#endif

#ifdef HAVE_TEUCHOS_LONG_DOUBLE
#  define TEUCHOS_MACRO_EXPAND_LONG_DOUBLE(INSTANT_MACRO)\
     INSTANT_MACRO(long double)
#else
#  define TEUCHOS_MACRO_EXPAND_LONG_DOUBLE(INSTANT_MACRO)
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
  TEUCHOS_MACRO_EXPAND_COMPLEX_DOUBLE(MACRONAME) \
  TEUCHOS_MACRO_EXPAND_LONG_DOUBLE(MACRONAME)


#endif  // TEUCHOS_EXPAND_SCALAR_TYPE_MACROS_HPP
