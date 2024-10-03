// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_StatusTest_Wrapper.H"
#include "LOCA_Solver_Wrapper.H"

LOCA::StatusTest::Wrapper::Wrapper(
            const Teuchos::RCP<NOX::StatusTest::Generic>& s) :
  statusTestPtr(s)
{
}

LOCA::StatusTest::Wrapper::~Wrapper()
{
}

NOX::StatusTest::StatusType
LOCA::StatusTest::Wrapper::checkStatus(const NOX::Solver::Generic& problem,
                       NOX::StatusTest::CheckType checkType)
{
  LOCA::Solver::Wrapper problemWrapper(Teuchos::rcp(&problem,false));

  return statusTestPtr->checkStatus(problemWrapper, checkType);
}

NOX::StatusTest::StatusType
LOCA::StatusTest::Wrapper::getStatus() const
{
  return statusTestPtr->getStatus();
}

std::ostream& LOCA::StatusTest::Wrapper::print(std::ostream& stream, int indent) const
{
  return statusTestPtr->print(stream, indent);
}

Teuchos::RCP<NOX::StatusTest::Generic>
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest()
{
  return statusTestPtr;
}

Teuchos::RCP<const NOX::StatusTest::Generic>
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest() const
{
  return statusTestPtr;
}
