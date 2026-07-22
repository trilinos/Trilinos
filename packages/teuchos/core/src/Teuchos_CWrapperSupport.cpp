// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CWrapperSupport_Cpp.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace {


Teuchos::RCP<Teuchos::FancyOStream>& printErrorOStreamImpl()
{
  static Teuchos::RCP<Teuchos::FancyOStream> printErrorOStream;
  if (is_null(printErrorOStream)) {
    printErrorOStream = Teuchos::VerboseObjectBase::getDefaultOStream();
  }
  return printErrorOStream;
}


bool& showStackTraceOnExceptionImpl()
{
  static bool showStackTraceOnException = false;
  return showStackTraceOnException;
}


} // namespace


namespace Teuchos {


//
// CWrapperErrorHandling
//


void CWrapperErrorHandling::setPrintErrorOStream(
  const RCP<FancyOStream> &errorOStream
  )
{
  printErrorOStreamImpl() = errorOStream;
}


RCP<FancyOStream> CWrapperErrorHandling::getPrintErrorOStream()
{
  return printErrorOStreamImpl();
}


void CWrapperErrorHandling::setShowStackTraceOnException(
  const bool showStrackTraceOnException)
{
  showStackTraceOnExceptionImpl() = showStrackTraceOnException;
}


bool CWrapperErrorHandling::getShowStackTraceOnException()
{
  return showStackTraceOnExceptionImpl();
}


} // namespace Teuchos
