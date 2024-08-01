// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TestForException.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

//
// ToDo: Make these functions thread-safe!
//


namespace {


int throwNumber = 0;


bool& loc_enableStackTrace()
{
  static bool static_enableStackTrace =
#ifdef HAVE_TEUCHOS_DEFAULT_STACKTRACE
    true
#else
    false
#endif
    ;
  return static_enableStackTrace;
}


} // namespace


void Teuchos::TestForException_incrThrowNumber()
{
  ++throwNumber;
}


int Teuchos::TestForException_getThrowNumber()
{
  return throwNumber;
}


void Teuchos::TestForException_break(const std::string &errorMsg, int throwNumber)
{
  (void)throwNumber; // Ignore unused arg
  // Provide a statement to break on
  size_t break_on_me;
  break_on_me = errorMsg.length(); // Use errMsg to avoid compiler warning.
  if (break_on_me) {} // Avoid warning
  // Above is just some statement for the debugger to break on.  Note: now is
  // a good time to examine the stack trace and look at the error message in
  // 'errorMsg' to see what happened.  In GDB just type 'where' or you can go
  // up by typing 'up' and moving up in the stack trace to see where you are
  // and how you got to this point in the code where you are throwing this
  // exception!  Typing in a 'p errorMsg' will show you what the error message
  // is.  Also, you should consider adding a conditional break-point in this
  // function based on a specific value of 'throwNumber' if the exception you
  // want to examine is not the first exception thrown.
}


void Teuchos::TestForException_setEnableStacktrace(bool enableStrackTrace)
{
  loc_enableStackTrace() = enableStrackTrace;
}


bool Teuchos::TestForException_getEnableStacktrace()
{
  return loc_enableStackTrace();
}

void Teuchos::TestForTermination_terminate(const std::string &msg) {
  std::ostringstream omsg;
  if (GlobalMPISession::getNProc() > 1) {
    omsg << "p="<<GlobalMPISession::getRank()<<": ";
  }
  omsg << msg << "\n";
  std::cerr << omsg.str();
  std::terminate();
}
// NOTE: The above usage of ostringstream is so that the output to std::cerr
// is done as one string.  This should help to avoid jumbled output like is
// occurring in tests that grep for this output (see #3163).
