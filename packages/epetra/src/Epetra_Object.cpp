
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
