// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_StatusTest_Generic.H"
#include "NOX_Common.H"
#include <iomanip>

std::ostream&
NOX::StatusTest::operator<<(std::ostream& os, NOX::StatusTest::StatusType type)
{
  os << std::setiosflags(std::ios::left) << std::setw(13) << std::setfill('.');
  switch (type) {
  case  NOX::StatusTest::Failed:
    os << "Failed";
    break;
  case  NOX::StatusTest::Converged:
    os << "Converged";
    break;
  case NOX::StatusTest::Unevaluated:
    os << "??";
    break;
  case  NOX::StatusTest::Unconverged:
  default:
    os << "**";
    break;
  }
  os << std::resetiosflags(std::ios::adjustfield) << std::setfill(' ');
  return os;
}
