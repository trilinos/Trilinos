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
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES // FIXME
  throw "Galeri::CreateVbrMatrix: VbrMatrix support not available in 64 bit global indices.";
  return 0;
#else

  const Epetra_Comm& Comm = CrsMatrix->Comm();
  const Epetra_Map& Map = CrsMatrix->RowMatrixRowMap();
  long long NumGlobalElements = 0;
  int NumMyElements = Map.NumMyElements();
  const int* MyGlobalElements_int = 0;
  const long long* MyGlobalElements_LL = 0;
  Epetra_BlockMap* BlockMap = 0;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map.GlobalIndicesInt()) {
    NumGlobalElements = Map.NumGlobalElements();
    MyGlobalElements_int = Map.MyGlobalElements();

    BlockMap = new Epetra_BlockMap((int) NumGlobalElements,
                                   NumMyElements,
                                   MyGlobalElements_int,
                                   NumPDEs, 0, Comm);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map.GlobalIndicesLongLong()) {
    NumGlobalElements = Map.NumGlobalElements64();
    MyGlobalElements_LL = Map.MyGlobalElements64();

    BlockMap = new Epetra_BlockMap(NumGlobalElements,
                                   NumMyElements,
                                   MyGlobalElements_LL,
                                   NumPDEs, 0, Comm);
  }
  else
#endif
    throw "Galeri::CreateVbrMatrix: GlobalIndices type unknown";

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

  std::vector<int>    VbrIndices(MaxNnzPerRow);
  std::vector<double> VbrValues(MaxBlockSize);
  int BlockRows = NumPDEs;
  int ierr;

  // cycle over all the local rows. 

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    long long GlobalNode = MyGlobalElements_int ? MyGlobalElements_int[i] : MyGlobalElements_LL[i];

    ierr = CrsMatrix->ExtractMyRowView(i, CrsNumEntries, CrsValues, CrsIndices);

    // CJ TODO FIXME : Change this when Epetra VBR matrix is changed for long long
	// Commenting out right now so the package builds when 32 bit GIDs are disabled.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    for (int kk = 0 ; kk < CrsNumEntries ; ++kk)
      VbrIndices[kk] = CrsMatrix->GCID(CrsIndices[kk]);

    VbrMatrix->BeginInsertGlobalValues(GlobalNode, CrsNumEntries, &VbrIndices[0]);
#else
    throw "Galeri::CreateVbrMatrix: not implemented fully";
#endif

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
#endif

} // CreateVbrMatrix()

} // namespace Galeri
