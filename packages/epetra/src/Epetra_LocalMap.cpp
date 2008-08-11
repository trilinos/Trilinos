
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

#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

//============================================================================
Epetra_LocalMap::Epetra_LocalMap(int NumMyElements, int IndexBase, 
																 const Epetra_Comm& Comm)
  // LocalMap is just a special case of Map
	: Epetra_Map(NumMyElements, NumMyElements, IndexBase, Comm) 
{
  SetLabel("Epetra::LocalMap");
  if (CheckInput()!=0)
    throw ReportError("Replicated Local Map not the same size on all PEs",-1);
}
//============================================================================
Epetra_LocalMap::Epetra_LocalMap(const Epetra_LocalMap& map)
	: Epetra_Map(map) 
{
  if (CheckInput()!=0)
    throw ReportError("Replicated Local Map not the same size on all PEs",-1);
}
 
//============================================================================
int Epetra_LocalMap::CheckInput() {
  int * tmp = new int[4];
  tmp[0] = NumMyElements();
  tmp[1] = - NumMyElements();
  Comm().MaxAll(tmp, tmp+2, 2);

  int tmp1 = tmp[2]; // Max of all NumMyElements across all processors
  int tmp2 = - tmp[3]; // Min of all ...
  delete [] tmp;

  if (tmp1==tmp2) 
		return(0);
  else 
		return(-1);
}
//=========================================================================
Epetra_LocalMap::~Epetra_LocalMap()
{
}
//=============================================================================
Epetra_LocalMap & Epetra_LocalMap::operator= (const Epetra_LocalMap& map) {
	if(this != &map)
		Epetra_BlockMap::operator=(map); // call this->Epetra_BlockMap::operator=
	return(*this);
}
