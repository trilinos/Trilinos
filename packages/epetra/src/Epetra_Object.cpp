
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_Object.h"


//=============================================================================
Epetra_Object::Epetra_Object(int TracebackModeIn) 
  : Label_(0)
{
  SetLabel("Epetra::Object");
  TracebackMode = (TracebackModeIn != -1) ? TracebackModeIn : TracebackMode;
}
//=============================================================================
Epetra_Object::Epetra_Object(const char * const Label, 
			     int TracebackModeIn) 
  : Label_(0)
{
  SetLabel(Label);
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
//=============================================================================
void Epetra_Object::Print(ostream & os) const {
  // os << Label_; // No need to print label, since ostream does it already
  return;
}
//=============================================================================
int Epetra_Object::ReportError(const string Message, int ErrorCode) const {

#ifndef EPETRA_NO_ERROR_REPORTS
  // NOTE:  We are extracting a C-style string from Message because
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"
                              if ((ErrorCode < 0 && Epetra_Object::GetTracebackMode() > 0) || \
                                  (ErrorCode > 0 && Epetra_Object::GetTracebackMode() > 1)) { \
  cerr << endl << "Error in Epetra Object with label:  " << Label_ << endl
       << "Epetra Error:  " << Message.c_str() << "  Error Code:  " << ErrorCode << endl;
}
#endif
  return(ErrorCode);
}
//=============================================================================
Epetra_Object::~Epetra_Object()  
{
  if (Label_!=0) delete [] Label_;
}
//=============================================================================
char * Epetra_Object::Label() const {
  return(Label_);
}
//=============================================================================
void Epetra_Object::SetLabel(const char * const Label) { 
  if (Label_!=0) delete [] Label_;
  Label_ = new char[strlen(Label)+1];
  strcpy(Label_,Label);
  return;
}
