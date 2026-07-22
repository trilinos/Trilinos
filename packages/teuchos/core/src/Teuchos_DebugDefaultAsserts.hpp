// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP
#define TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP


#include "Teuchos_Assert.hpp"


#ifdef TEUCHOS_DEBUG
#define TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT() default: TEUCHOS_TEST_FOR_EXCEPT(true)
#else
/** \brief Macro to insert switch default that throws in a debug build.
 *
 * In a non-debug build, however, this does nothing.  This is also helpful for
 * removing code that would otherwise show as not being covered.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT()
#endif

// NOTE: Some explaination for the above TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT()
// macro:
//
// In a debug build where default: throws always, we can't have a follow-up
// break statement or some compilers (e.g. NVCC) will (correctly) complain
// that the 'break' statement is unreachable.  Older compilers did not do
// this.  However, in a non-debug build, you want to not even put in a
// 'default:' block.  That way, if all of the enum values are not covered,
// then must compilers will issue a warning about that and we want to see that
// warning in a non-debug build.

#ifdef TEUCHOS_DEBUG
#define TEUCHOS_IF_ELSE_DEBUG_ASSERT() else { TEUCHOS_TEST_FOR_EXCEPT(true); }
#else
/** \brief Macro to insert else block that throws in a debug build.
 *
 * In a non-debug build, however, this does nothing.  This is also helpful for
 * removing code that would otherwise show as not being covered.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_IF_ELSE_DEBUG_ASSERT() else {}
#endif


#endif // TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP
