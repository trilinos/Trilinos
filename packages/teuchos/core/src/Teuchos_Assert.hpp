// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_ASSERT_HPP
#define TEUCHOS_ASSERT_HPP


#include "Teuchos_TestForException.hpp"


/** \brief This macro is throws when an assert fails.
 *
 * \note <tt>The std::exception</tt> thrown is <tt>std::logic_error</tt>.
 *
 * \ingroup TestForException_grp
 */
#define TEUCHOS_ASSERT(assertion_test) TEUCHOS_TEST_FOR_EXCEPT(!(assertion_test))

#ifdef TEUCHOS_DEBUG
/** \brief Behaves as TEUCHOS_ASSERT only if debugging is enabled.
 *
 * \ingroup TestForException_grp
 */
#define TEUCHOS_DEBUG_ASSERT(assertion_test) TEUCHOS_ASSERT(assertion_test)
#else
#define TEUCHOS_DEBUG_ASSERT(assertion_test)
#endif


/** \brief This macro asserts that an integral number fallis in the range
 * <tt>[lower_inclusive,upper_exclusive)</tt>
 *
 * \note <tt>The std::exception</tt> thrown is <tt>std::out_of_range</tt>.
 *
 * WARNING: This assert will evaluate <tt>index</tt>,
 * <tt>lower_inclusive</tt>, and <tt>upper_inclusive</tt> more than once if
 * there is a failure which will cause the side-effect of an additional
 * evaluation.  This is needed because the return types of these values are
 * unknown.  Therefore, only pass in arguments that are objects or function
 * calls that have not side-effects!
 *
 * \ingroup TestForException_grp
 */
#define TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( index, lower_inclusive, upper_exclusive ) \
  { \
    TEUCHOS_TEST_FOR_EXCEPTION( \
      !( (lower_inclusive) <= (index) && (index) < (upper_exclusive) ), \
      std::out_of_range, \
      "Error, the index " #index " = " << (index) << " does not fall in the range" \
      "["<<(lower_inclusive)<<","<<(upper_exclusive)<<")!" ); \
  }


/** \brief This macro is checks that to numbers are equal and if not then
 * throws an exception with a good error message.
 *
 * \note The <tt>std::exception</tt> thrown is <tt>std::out_of_range</tt>.
 *
 * WARNING: This assert will evaluate <tt>val1</tt> and <tt>val2</tt> more
 * than once if there is a failure which will cause the side-effect of an
 * additional evaluation.  This is needed because the return types of
 * <tt>val1</tt> and <tt>val2</tt> are unknown.  Therefore, only pass in
 * arguments that are objects or function calls that have not side-effects!
 *
 * \ingroup TestForException_grp
 */
#define TEUCHOS_ASSERT_EQUALITY( val1, val2 ) \
  { \
    TEUCHOS_TEST_FOR_EXCEPTION( \
      (val1) != (val2), std::out_of_range, \
      "Error, (" #val1 " = " << (val1) << ") != (" #val2 " = " << (val2) << ")!" ); \
  }


/** \brief This macro is checks that an inequality between two numbers is
 * satisified and if not then throws a good exception message.
 *
 * \note The <tt>std::exception</tt> thrown is <tt>std::out_of_range</tt>.
 *
 * WARNING: This assert will evaluate <tt>val1</tt> and <tt>val2</tt> more
 * than once if there is a failure which will cause the side-effect of an
 * additional evaluation.  This is needed because the return types of
 * <tt>val1</tt> and <tt>val2</tt> are unknown.  Therefore, only pass in
 * arguments that are objects or function calls that have not side-effects!
 *
 * \ingroup TestForException_grp
 */
#define TEUCHOS_ASSERT_INEQUALITY( val1, comp, val2 ) \
  { \
    TEUCHOS_TEST_FOR_EXCEPTION( \
      !( (val1) comp (val2) ), std::out_of_range, \
      "Error, (" #val1 " = " << (val1) << ") " \
      #comp " (" #val2 " = " << (val2) << ")! FAILED!" ); \
  }


#endif // TEUCHOS_ASSERT_HPP
