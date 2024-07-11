// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_StatusTest_Abstract.H"
#include <iomanip>

std::ostream&
LOCA::StatusTest::operator<<(std::ostream& os, LOCA::StatusTest::StatusType status)
{
  os << std::setiosflags(std::ios::left) << std::setw(13) << std::setfill('.');
  switch (status) {
  case  LOCA::StatusTest::Finished:
    os << "Finished";
    break;
  case  LOCA::StatusTest::Failed:
    os << "Failed";
    break;
  case LOCA::StatusTest::NotFinished:
    os << "Not finished";
    break;
  default:
    os << "**";
    break;
  }
  os << std::resetiosflags(std::ios::adjustfield) << std::setfill(' ');
  return os;
}
