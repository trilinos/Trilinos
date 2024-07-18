// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_UNIT_TEST_HELPERS_HPP
#define THYRA_UNIT_TEST_HELPERS_HPP


/*! \file Thyra_UnitTestHelpers.hpp

\brief Thyra-specific macros for helping to create concrete unit tests.
*/


#include "Thyra_ConfigDefs.hpp"
#include "Teuchos_UnitTestHelpers.hpp"


#ifdef HAVE_TEUCHOS_INST_FLOAT
#  define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME)\
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, float)
#else
#  define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME)
#endif

#define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_DOUBLE(TEST_GROUP, TEST_NAME)\
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, double)

#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT)
#  define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME)\
     typedef std::complex<float> ComplexFloat; \
     TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, ComplexFloat)
#else
#  define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME)
#endif

#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
#  define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME)\
     typedef std::complex<double> ComplexDouble; \
     TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, ComplexDouble)
#else
#  define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME)
#endif


/** \brief Instantiate a whole group of tests for supported real Scalar
 * types.
 */
#define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(TEST_GROUP, TEST_NAME)\
  THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME) \
  THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_DOUBLE(TEST_GROUP, TEST_NAME)


/** \brief Instantiate a whole group of tests for supported Scalar types. */
#define THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(TEST_GROUP, TEST_NAME)\
  THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME) \
  THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_DOUBLE(TEST_GROUP, TEST_NAME) \
  THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME) \
  THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME)


#endif  // THYRA_UNIT_TEST_HELPERS_HPP
