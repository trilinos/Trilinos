// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_toString.hpp"


namespace {


// Make sure this is initialized whenever needed before main starts!
bool& showTestFailureLocationImpl()
{
  static bool showTFL = false;
  return showTFL;
}


} // namespace


const std::string
Teuchos::passfail_with_location(const bool result,
  const std::string &file, const int lineNumber)
{
  std::string rtn = passfail(result);
  if (!result && showTestFailureLocation()) {
    rtn += " ==> "+file+":"+toString(lineNumber);
  }
  return rtn;
}


void Teuchos::showTestFailureLocation(bool showTFL)
{
  showTestFailureLocationImpl() = showTFL;
}


bool Teuchos::showTestFailureLocation()
{
  return showTestFailureLocationImpl();
}
