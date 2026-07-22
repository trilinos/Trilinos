// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Ptr.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Exceptions.hpp"


void Teuchos::PtrPrivateUtilityPack::throw_null( const std::string &type_name )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, NullReferenceError,
    "Ptr<"<<type_name<<">::assert_not_null() : You can not"
    " call operator->() or operator*() if get()==NULL!" );
}
