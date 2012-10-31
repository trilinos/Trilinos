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

namespace Teuchos
{
//=============================================================================
Object::Object(int tracebackModeIn) : label_(0)
{
  setLabel("Teuchos::Object");
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}
//=============================================================================
Object::Object(const char* label_in, int tracebackModeIn) : label_(0)
{
  setLabel(label_in);
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}
//=============================================================================
Object::Object(const Object& Obj) : label_(0)
{
  setLabel(Obj.label());
}
// Set TracebackMode value to default
int Object::tracebackMode(-1);

void Object::setTracebackMode(int tracebackModeValue)
{
  if (tracebackModeValue < 0)
    tracebackModeValue = 0;
  Object tempObject(tracebackModeValue);
}

int Object::getTracebackMode()
{
  int temp = Object::tracebackMode;
  if (temp == -1)
    temp = Teuchos_DefaultTracebackMode;
  return(temp);
}
//=============================================================================
void Object::print(std::ostream& os) const
{
  // os << label_; // No need to print label, since std::ostream does it already
}
//=============================================================================
Object::~Object()  
{
  if (label_!=0) {
    delete [] label_;
		label_ = 0;
	}
}
//=============================================================================
char* Object::label() const
{
  return(label_);
}
//=============================================================================
void Object::setLabel(const char* label_in)
{ 
  if (label_ != 0)
    delete [] label_;
  label_ = new char[std::strlen(label_in) + 1];
  std::strcpy(label_, label_in);
}
} // namespace Teuchos
