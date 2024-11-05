// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_COMPILE_TIME_ASSERT_HPP
#define TEUCHOS_COMPILE_TIME_ASSERT_HPP

/*! \file Teuchos_CompileTimeAssert.hpp
    \brief Template classes for testing assertions at compile time.
*/

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

/*! \defgroup CompileTimeAssert_grp  Template classes for testing assertions at compile time.
 \ingroup teuchos_language_support_grp
*/
///@{

/// If instantiated (for Test!=0) then this should not compile!
template <int Test>
class CompileTimeAssert {
	int compile_time_assert_failed[Test-1000]; // Should not compile if instantiated!
};

/// If instantiated (i.e. Test==0) then this will compile!
template <>
class CompileTimeAssert<0> {};

///@}

} // namespace Teuchos

#endif // TEUCHOS_COMPILE_TIME_ASSERT_HPP
