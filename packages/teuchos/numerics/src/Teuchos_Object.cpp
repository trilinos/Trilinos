// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Object.hpp"

namespace Teuchos {

// Set TracebackMode value to default.
int Object::tracebackMode = -1;

Object::Object (int tracebackModeIn)
{
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}

Object::Object (const char* label_in, int tracebackModeIn) :
  label_ (label_in) // this does a deep copy of the input string
{
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}

Object::Object (const std::string& label_in, int tracebackModeIn) :
  label_ (label_in)
{
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}

void Object::setLabel (const char* theLabel) {
  label_ = std::string (theLabel);
}

void Object::setTracebackMode (int tracebackModeValue)
{
  if (tracebackModeValue < 0) {
    tracebackModeValue = 0;
  }
  Object tempObject (tracebackModeValue);
}

int Object::getTracebackMode()
{
  int temp = Object::tracebackMode;
  if (temp == -1) {
    temp = Teuchos_DefaultTracebackMode;
  }
  return(temp);
}

void Object::print (std::ostream& /* os */) const
{
  // os << label_; // No need to print label, since std::ostream does it already
}

int Object::reportError (const std::string message, int errorCode) const
{
  using std::cerr;
  using std::endl;

  // mfh 23 Nov 2014: I found the following message here:
  //
  // NOTE:  We are extracting a C-style std::string from Message because
  //        the SGI compiler does not have a real std::string class with
  //        the << operator.  Some day we should get rid of ".c_str()"
  //
  // All the compilers we support now have a correct implementation of
  // std::string, so I've corrected the code below to use std::string
  // instead.

  if (tracebackMode == 1 && errorCode < 0) {
    // Report fatal error
    cerr << endl << "Error in Teuchos Object with label: " << label_
         << endl << "Teuchos Error:  " << message << "  Error Code:  "
         << errorCode << endl;
    return errorCode;
  }
  if (tracebackMode == 2 && errorCode != 0) {
    cerr << endl << "Error in Teuchos Object with label: " << label_
         << endl << "Teuchos Error:  " << message << "  Error Code:  "
         << errorCode << endl;
    return errorCode;
  }
  return errorCode;
}

const char* Object::label () const
{
  return label_.c_str ();
}

std::ostream& operator<< (std::ostream& os, const Teuchos::Object& obj)
{
  os << obj.label () << std::endl;
  obj.print (os);
  return os;
}

} // namespace Teuchos
