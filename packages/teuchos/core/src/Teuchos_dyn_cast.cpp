// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

/** We throw a m_bad_cast, which is a subclass of bad_cast.
	This is necessary, since bad_cast lacks the appropriate
	constructor for use with the TEUCHOS_TEST_FOR_EXCEPTION macro.
*/
void Teuchos::dyn_cast_throw_exception(
  const std::string &T_from,
  const std::string &T_from_concr,
  const std::string &T_to
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, m_bad_cast
    ,"dyn_cast<" << T_to << ">(" << T_from
    << ") : Error, the object with the concrete type \'"
    << T_from_concr << "\' (passed in through the interface type \'" << T_from <<  "\') "
    " does not support the interface \'"
    << T_to << "\' and the dynamic cast failed!" );
}
