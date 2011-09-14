
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

#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

//============================================================================
Epetra_LocalMap::Epetra_LocalMap(int numMyElements, int indexBase, 
				 const Epetra_Comm& comm)
  // LocalMap is just a special case of Map
	: Epetra_Map(numMyElements, numMyElements, indexBase, comm) 
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
