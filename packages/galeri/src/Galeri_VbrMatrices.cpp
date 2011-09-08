// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
