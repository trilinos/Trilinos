
//@HEADER
// ***********************************************************************
// 
//                     New_Package Example Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "Newp_Jambo.h"
#include "Epetra_Comm.h"
//=============================================================================
Newp_Jambo::Newp_Jambo(const Epetra_Comm& Comm):Comm_(Comm) {} 

//=======================================================================
void Newp_Jambo::Print(ostream& os) const {
  int MyPID = Comm_.MyPID();
  int NumProc = Comm_.NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {

      if (MyPID==0) {
	os <<  "This will print out one line for each of the " << NumProc << " processes " << endl ; 
	os << endl;
      }
      os << "Jambo.  I am process " << MyPID << endl ; 
      os << flush; 
    }

    // Do a few global ops to give I/O a chance to complete
    Comm_.Barrier();
    Comm_.Barrier();
    Comm_.Barrier();
  }
  return; }


