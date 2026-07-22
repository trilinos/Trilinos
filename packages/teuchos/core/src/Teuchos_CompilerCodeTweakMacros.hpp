// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP
#define TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP

#include "TeuchosCore_config.h"


//
// Implementations of below macros
//


#ifdef __NVCC__
  #define TEUCHOS_UNREACHABLE_RETURN_IMPL(dummyReturnVal)
#else
  #define TEUCHOS_UNREACHABLE_RETURN_IMPL(dummyReturnVal) \
    return dummyReturnVal
#endif
// Above, the unreachable return statement is only removed for __NVCC__.  This
// is because only __NVCC__ warns about the unreachable return.  Having that
// return statement added for other compilers that don't warn leads to safer
// code by avoiding undefined behavior in case a function is modified in such
// a way that this return is no longer unreachable.  In that case, we want it
// to return something to avoid undefined behavior in case it ever gets
// executed.  That is just good defensive programming.  Also, by leaving this
// return statement by default, we avoid the (incorrect) warning by some
// compilers that the function may not return a value.


//
// User interface macros
//


/** \brief Avoid warning about unreachable or missing return from function.
 *
 * Consider a function like:
 *
 * \code
  int func(const ESomeEnum val)
  {
    switch (val) {
      case VAL1: return 1;
      case VAL2: return 2;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
 * \endcode
 *
 * That code will never execute out of the switch statement.  However, some
 * compilers will provide a warning that the function may not return a value.
 * Therefore, one can remove this warning by adding a dummy return value like:
 *
 * \code
  int func(const ESomeEnum val)
  {
    switch (val) {
      case VAL1: return 1;
      case VAL2: return 2;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    return -1; // Will never get called!
  }
 * \endcode
 *
 * That removes the "may not return value" warning on those compilers.  But
 * other compilers will correctly warn that <tt>return -1;</tt> will never be
 * executed with a warning like "statement is unreachable".  Therefore, to
 * address warnings like this, this macro is used like:
 *
 * \code
  int func(const ESomeEnum val)
  {
    switch (val) {
      case VAL1: return 1;
      case VAL2: return 2;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }
 * \endcode
 *
 * On compilers that warn about the return being unreachable the return
 * statement is skipped.  On every other compiler, the return statement is
 * kept which results in safer code under refactoring (by avoiding undefined
 * behavior when returning from a function by fall-through without returning
 * an explicit value.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_UNREACHABLE_RETURN(dummyReturnVal) \
  TEUCHOS_UNREACHABLE_RETURN_IMPL(dummyReturnVal)


#endif // TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP
