
/* Copyright (2003) Sandia Corportation. Under the terms of Contract 
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


