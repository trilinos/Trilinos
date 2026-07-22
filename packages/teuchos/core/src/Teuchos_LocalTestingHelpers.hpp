// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_LOCAL_TESTING_HELPERS_HPP
#define TEUCHOS_LOCAL_TESTING_HELPERS_HPP


/*! \file Teuchos_LocalTestingHelpers.hpp
 *
 * \brief Utilities to make writing tests easier.
 *
 * <b>WARNING!</b> These macros are not namespaced, so you must
 * only include it in *.cpp files for unit testing only!
 *
 */


#include "Teuchos_TestingHelpers.hpp"


/** \brief Echo the given statement before it is executed.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define ECHO( statement ) \
  TEUCHOS_ECHO( statement, out )


/** \brief Assert the given statement is true.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_ASSERT( v1 ) \
  TEUCHOS_TEST_ASSERT( v1, out, success )


/** \brief Assert the equality of v1 and constant v2.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_EQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success )


/** \brief Assert the equality of v1 and v2.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_EQUALITY( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY( v1, v2, out, success )


/** \brief Assert the inequality of v1 and constant v2.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_INEQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_INEQUALITY_CONST( v1, v2, out, success )


/** \brief Assert the inequality of v1 and v2.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_INEQUALITY( v1, v2 ) \
  TEUCHOS_TEST_INEQUALITY( v1, v2, out, success )


/** \brief Assert the relative floating-point equality of
 * rel_error(v1,v2) <= tol.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_FLOATING_EQUALITY( v1, v2, tol ) \
  TEUCHOS_TEST_FLOATING_EQUALITY( v1, v2, tol, out, success )


/** \brief Assert that two iterators are equal.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_ITER_EQUALITY( iter1, iter2 ) \
  TEUCHOS_TEST_ITER_EQUALITY( iter1, iter2, out, success )


/** \brief Assert that two iterators are NOT equal.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_ITER_INEQUALITY( iter1, iter2 ) \
  TEUCHOS_TEST_ITER_INEQUALITY( iter1, iter2, out, success )


/** \brief Assert that a[i] == val.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_ARRAY_ELE_EQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, false, out, local_success )


/** \brief Assert that a[i] != val.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_ARRAY_ELE_INEQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_INEQUALITY( a, i, val, false, out, local_success )


/** \brief Assert that v1 comp v2 (where comp = '==', '>=", "!=", etc).
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_COMPARE( v1, comp, v2 ) \
  TEUCHOS_TEST_COMPARE( v1, comp, v2, out, success )


/** \brief Assert that v1 comp v2 (where comp = '==', '>=", "!=", etc) where
 * the second object v2 is printed as value.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_COMPARE_CONST( v1, comp, v2 ) \
  TEUCHOS_TEST_COMPARE_CONST( v1, comp, v2, out, success )


/** \brief Assert that a1.size()==a2.size() and a[i]==b[i], i=0....
 *
 * Works for all object types that support a1[i], a1.size(), a2[j], and
 * a2.size() and types a1 and a2 can be different types!
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_COMPARE_ARRAYS( a1, a2 ) \
  { \
    const bool l_result = compareArrays(a1,#a1,a2,#a2,out); \
    if (!l_result) success = false; \
  }


/** \brief Assert that a1.size()==a2.size() and rel_error(a[i],b[i]) <= tol, i=0....
 *
 * Works for all object types that support a1[i], a1.size(), a2[j], and
 * a2.size() and types a1 and a2 can be different types!
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_COMPARE_FLOATING_ARRAYS( a1, a2, tol ) \
  { \
    const bool result = compareFloatingArrays(a1,#a1,a2,#a2,tol,out); \
    if (!result) success = false; \
  }

/** \brief Assert that a1.size()==a2.size() and |a[i]-b[i]| <= tol, i=0....
 *
 * Works for all object types that support a1[i], a1.size(), a2[j], and
 * a2.size(), but the element types of a1 and a2 must be the same and Teuchos::ScalarTraits
 * must have a specialization for this element type.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_ABSOLUTE_COMPARE_FLOATING_ARRAYS( a1, a2, tol ) \
  { \
    const bool result = compareFloatingArraysAbsolute(a1,#a1,a2,#a2,tol,out); \
    if (!result) success = false; \
  }

/** \brief Assert that the statement 'code' throws the exception 'ExceptType'
 * (otherwise the test fails).
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_THROW( code, ExceptType  ) \
  TEUCHOS_TEST_THROW( code, ExceptType, out, success  )


/** \brief Asserr that the statement 'code' does *not* thrown any excpetions.
 *
 * \ingroup Teuchos_UnitTestAssertMacros_grp
 */
#define TEST_NOTHROW( code  ) \
  TEUCHOS_TEST_NOTHROW( code, out, success  )


#endif  // TEUCHOS_LOCAL_TESTING_HELPERS_HPP
