// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestBase.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


namespace Teuchos {


UnitTestBase::UnitTestBase(const std::string groupName, std::string testName)
{
  Teuchos::UnitTestRepository::addUnitTest(this, groupName, testName);
}


bool UnitTestBase::runUnitTest( FancyOStream &out ) const
{
  bool success = true;
  try {
    runUnitTestImpl(out, success);
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, out, success)
  return success;
}


} // namespace Teuchos
