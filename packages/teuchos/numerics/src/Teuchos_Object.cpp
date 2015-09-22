// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

void Object::print (std::ostream& os) const
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
