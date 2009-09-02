// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"
#include "Galeri_VbrMatrices.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include <vector>

namespace Galeri {

Epetra_VbrMatrix* 
CreateVbrMatrix(const Epetra_CrsMatrix* CrsMatrix, const int NumPDEs)
{
  const Epetra_Comm& Comm = CrsMatrix->Comm();
  const Epetra_Map& Map = CrsMatrix->RowMatrixRowMap();

  int NumGlobalElements = Map.NumGlobalElements();
  int NumMyElements = Map.NumMyElements();
  int* MyGlobalElements = Map.MyGlobalElements();

  Epetra_BlockMap* BlockMap = new Epetra_BlockMap(NumGlobalElements,
                                                  NumMyElements,
                                                  MyGlobalElements,
                                                  NumPDEs, 0, Comm);

  int MaxNnzPerRow = CrsMatrix->MaxNumEntries();
  // create a VBR matrix based on BlockMap
  Epetra_VbrMatrix* VbrMatrix = new Epetra_VbrMatrix(Copy, *BlockMap,
                                                     MaxNnzPerRow);

  delete BlockMap;

  // size of each VBR block
  int MaxBlockSize = NumPDEs * NumPDEs;

  int CrsNumEntries;
  int* CrsIndices;
  double* CrsValues;

  vector<int>    VbrIndices(MaxNnzPerRow);
  vector<double> VbrValues(MaxBlockSize);
  int BlockRows = NumPDEs;
  int ierr;

  // cycle over all the local rows. 

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int GlobalNode = MyGlobalElements[i];

    ierr = CrsMatrix->ExtractMyRowView(i, CrsNumEntries, CrsValues, CrsIndices);

    for (int kk = 0 ; kk < CrsNumEntries ; ++kk)
      VbrIndices[kk] = CrsMatrix->GCID(CrsIndices[kk]);

    VbrMatrix->BeginInsertGlobalValues(GlobalNode, CrsNumEntries, &VbrIndices[0]);

    for (int i = 0 ; i < CrsNumEntries ; ++i) 
    {
      for (int k = 0 ; k < BlockRows ; ++k) 
      {
        for (int h = 0 ; h < BlockRows ; ++h)
        {
          if (k == h) VbrValues[k + h * BlockRows] = CrsValues[i];
          else        VbrValues[k + h * BlockRows] = 0.0;
        }
      }
      VbrMatrix->SubmitBlockEntry(&VbrValues[0], BlockRows, BlockRows, BlockRows);
    }

    VbrMatrix->EndSubmitEntries();
  }

  VbrMatrix->FillComplete();

  return(VbrMatrix);

} // CreateVbrMatrix()

} // namespace Galeri
