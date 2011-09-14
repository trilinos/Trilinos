
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

#include "Epetra_Object.h"

#ifdef HAVE_EPETRA_TEUCHOS
#include "Teuchos_VerboseObject.hpp"
#endif

//=============================================================================
Epetra_Object::Epetra_Object(int TracebackModeIn, bool set_label) 
  : Label_(0)
{
  if (set_label) {
    SetLabel("Epetra::Object");
  }
  TracebackMode = (TracebackModeIn != -1) ? TracebackModeIn : TracebackMode;
}

//=============================================================================
Epetra_Object::Epetra_Object(const char * const Label_in, 
			     int TracebackModeIn) 
  : Label_(0)
{
  SetLabel(Label_in);
  TracebackMode = (TracebackModeIn != -1) ? TracebackModeIn : TracebackMode;
}
//=============================================================================
Epetra_Object::Epetra_Object(const Epetra_Object& Object)
  : Label_(0)
{
  SetLabel(Object.Label_);
}

// Set TracebackMode value to default
int Epetra_Object::TracebackMode(-1);

void Epetra_Object::SetTracebackMode(int TracebackModeValue) {
  if (TracebackModeValue < 0) TracebackModeValue = 0;
  Epetra_Object TempObject(TracebackModeValue);
}

int Epetra_Object::GetTracebackMode() {
  int temp = Epetra_Object::TracebackMode;
  if (temp==-1) temp = DefaultTracebackMode;
  return(temp);
}

std::ostream& Epetra_Object::GetTracebackStream() {
#ifdef HAVE_EPETRA_TEUCHOS
  return (*Teuchos::VerboseObjectBase::getDefaultOStream());
#else
  return cerr;
#endif
}

//=============================================================================
void Epetra_Object::Print(ostream & os) const {
  (void)os;//prevents unused variable compiler-warning
  // os << Label_; // No need to print label, since ostream does it already
  return;
}
//=============================================================================
int Epetra_Object::ReportError(const string Message, int ErrorCode) const {
#ifndef EPETRA_NO_ERROR_REPORTS
  // NOTE:  We are extracting a C-style string from Message because
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"
  if (
    ( ErrorCode < 0 && Epetra_Object::GetTracebackMode() > 0 )
    ||
    ( ErrorCode > 0 && Epetra_Object::GetTracebackMode() > 1 )
    )
  {
    GetTracebackStream()
      << endl << "Error in Epetra Object with label:  " << Label_ << endl
      << "Epetra Error:  " << Message.c_str() << "  Error Code:  " << ErrorCode << endl;
  }
#endif
  return(ErrorCode);
}
//=============================================================================
Epetra_Object::~Epetra_Object()  
{
  if (Label_!=0) {
    delete [] Label_;
    Label_ = 0;
  }
}
//=============================================================================
const char * Epetra_Object::Label() const {
  return(Label_);
}
//=============================================================================
void Epetra_Object::SetLabel(const char * const Label_in)
{ 
  if (Label_!=0) {
    delete [] Label_;
    Label_ = 0;
  }
  if (Label_in==0) return;
  Label_ = new char[std::strlen(Label_in)+1];
  std::strcpy(Label_,Label_in);
  return;
}
