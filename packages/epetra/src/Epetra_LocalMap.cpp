
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

#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"


//============================================================================
Epetra_LocalMap::Epetra_LocalMap(int NumMyElements, int IndexBase, 
			       const Epetra_Comm& Comm)
  // LocalMap is just a special case of Map
: Epetra_Map(NumMyElements, NumMyElements, IndexBase, Comm) {
  SetLabel("Epetra::LocalMap");
  if (CheckInput()!=0)
    throw ReportError("Replicated Local Map not the same size on all PEs",-1);
}
//============================================================================
Epetra_LocalMap::Epetra_LocalMap(const Epetra_LocalMap& map)
: Epetra_Map(map) {
  if (CheckInput()!=0)
    throw ReportError("Replicated Local Map not the same size on all PEs",-1);
}
 
//============================================================================
int Epetra_LocalMap::CheckInput() {
  int * tmp = new int[4];
  tmp[0] = NumMyElements_;
  tmp[1] = - NumMyElements_;
  Comm().MaxAll(tmp, tmp+2, 2);

  int tmp1 = tmp[2]; // Max of all NumMyElements across all processors
  int tmp2 = - tmp[3]; // Min of all ...
  delete [] tmp;

  if (tmp1==tmp2) return(0);
  else return(-1);
}
//=========================================================================
Epetra_LocalMap::~Epetra_LocalMap(){}
