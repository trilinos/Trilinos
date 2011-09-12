//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "EpetraExt_BlockMultiVector.h"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"

using std::vector;

namespace EpetraExt {

//=============================================================================
// EpetraExt::BlockMultiVector Constructor
BlockMultiVector::BlockMultiVector(
      const Epetra_BlockMap & BaseMap,
      const Epetra_BlockMap & GlobalMap,
      int NumVectors )
  : Epetra_MultiVector( GlobalMap, NumVectors ),
    BaseMap_( BaseMap ),
    Offset_( BlockUtility::CalculateOffset( BaseMap ) )
{
}

//=============================================================================
// EpetraExt::BlockMultiVector Constructor
BlockMultiVector::BlockMultiVector(
      Epetra_DataAccess CV, 
      const Epetra_BlockMap & BaseMap, 
      const Epetra_MultiVector & BlockVec)
  : Epetra_MultiVector( CV, BlockVec, 0, BlockVec.NumVectors() ),
    BaseMap_( BaseMap ),
    Offset_( BlockUtility::CalculateOffset( BaseMap ) )
{
}

//==========================================================================
// Copy Constructor
BlockMultiVector::BlockMultiVector(const BlockMultiVector& Source)
  : Epetra_MultiVector( dynamic_cast<const Epetra_MultiVector &>(Source) ),
    BaseMap_( Source.BaseMap_ ),
    Offset_( Source.Offset_ )
{
}

//=========================================================================
BlockMultiVector::~BlockMultiVector()
{
}

//=========================================================================
int BlockMultiVector::ExtractBlockValues(Epetra_MultiVector & BaseVector, int GlobalBlockRow) const
{
   int IndexOffset = GlobalBlockRow * Offset_;
   int localIndex=0;

   // For each entry in the base vector, translate its global ID
   // by the IndexOffset and extract the value from this blockVector
   for (int i=0; i<BaseMap_.NumMyElements(); i++) {
      localIndex = this->Map().LID((IndexOffset + BaseMap_.GID(i)));
      if (localIndex==-1) { 
	     cout << "Error in  BlockMultiVector::GetBlock: " << i << " " 
		  << IndexOffset << " " << BaseMap_.GID(i) << endl;
	     return -1;
      }
      for (int j=0; j<NumVectors(); j++)
	BaseVector[j][i] = (*this)[j][localIndex]; 
   }

   return 0;
}

//=========================================================================
int BlockMultiVector::LoadBlockValues(const Epetra_MultiVector & BaseVector, int GlobalBlockRow) 
{
   int IndexOffset = GlobalBlockRow * Offset_;
   int localIndex=0;

   // For each entry in the base vector, translate its global ID
   // by the IndexOffset and load into this blockVector
   for (int i=0; i<BaseMap_.NumMyElements(); i++) {
      localIndex = this->Map().LID((IndexOffset + BaseMap_.GID(i)));
      if (localIndex==-1) { 
	     cout << "Error in  BlockMultiVector::GetBlock: " << i << " " 
		  << IndexOffset << " " << BaseMap_.GID(i) << endl;
	     return -1;
      }
      for (int j=0; j<NumVectors(); j++)
	(*this)[j][localIndex] = BaseVector[j][i];
   }

   return 0;
}

//=========================================================================
Teuchos::RCP<const Epetra_MultiVector>
BlockMultiVector::GetBlock(int GlobalBlockRow) const
{
  int offset = GlobalBlockRow * BaseMap_.NumMyElements();
  int numVecs = NumVectors();
  double **pointers = Pointers();
  double **block_pointers = new double*[numVecs];
  for (int i=0; i<numVecs; i++)
    block_pointers[i] = pointers[i]+offset;
  Teuchos::RCP<Epetra_MultiVector> block = 
    Teuchos::rcp(new Epetra_MultiVector(View, BaseMap_, block_pointers,
					numVecs));
  delete [] block_pointers;
  return block;
}

//=========================================================================
Teuchos::RCP<Epetra_MultiVector>
BlockMultiVector::GetBlock(int GlobalBlockRow)
{
  int offset = GlobalBlockRow * BaseMap_.NumMyElements();
  int numVecs = NumVectors();
  double **pointers = Pointers();
  double **block_pointers = new double*[numVecs];
  for (int i=0; i<numVecs; i++)
    block_pointers[i] = pointers[i]+offset;
  Teuchos::RCP<Epetra_MultiVector> block = 
    Teuchos::rcp(new Epetra_MultiVector(View, BaseMap_, block_pointers,
					numVecs));
  delete [] block_pointers;
  return block;
}

//=========================================================================
const Epetra_BlockMap&
BlockMultiVector::GetBaseMap() const
{
  return BaseMap_;
}

} //namespace EpetraExt
