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

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

namespace EpetraExt {

//=============================================================================
// EpetraExt::BlockVector Constructor
BlockVector::BlockVector(
      const Epetra_BlockMap & BaseMap,
      const Epetra_BlockMap & GlobalMap)
  : Epetra_Vector( GlobalMap ),
    BaseMap_( BaseMap ),
    Offset_( BlockUtility::CalculateOffset( BaseMap ) )
{
}

//=============================================================================
// EpetraExt::BlockVector Constructor
BlockVector::BlockVector(
      Epetra_DataAccess CV, 
      const Epetra_BlockMap & BaseMap, 
      const Epetra_Vector & BlockVec)
  : Epetra_Vector( CV, BlockVec, 0 ),
    BaseMap_( BaseMap ),
    Offset_( BlockUtility::CalculateOffset( BaseMap ) )
{
}

//==========================================================================
// Copy Constructor
BlockVector::BlockVector(const BlockVector& Source)
  : Epetra_Vector( dynamic_cast<const Epetra_Vector &>(Source) ),
    BaseMap_( Source.BaseMap_ ),
    Offset_( Source.Offset_ )
{
}

//=========================================================================
BlockVector::~BlockVector()
{
}

//=========================================================================
int BlockVector::ExtractBlockValues(Epetra_Vector & BaseVector, int GlobalBlockRow) const
{
   int IndexOffset = GlobalBlockRow * Offset_;
   int localIndex=0;

   // For each entry in the base vector, translate its global ID
   // by the IndexOffset and extract the value from this blockVector
   for (int i=0; i<BaseMap_.NumMyElements(); i++) {
      localIndex = this->Map().LID((IndexOffset + BaseMap_.GID(i)));
      if (localIndex==-1) { 
	     cout << "Error in  BlockVector::GetBlock: " << i << " " 
		  << IndexOffset << " " << BaseMap_.GID(i) << endl;
	     return -1;
      }
      BaseVector[i] = Values_[localIndex]; 
   }

   return 0;
}

//=========================================================================
int BlockVector::LoadBlockValues(const Epetra_Vector & BaseVector, int GlobalBlockRow) 
{
   int IndexOffset = GlobalBlockRow * Offset_;
   int localIndex=0;

   // For each entry in the base vector, translate its global ID
   // by the IndexOffset and load into this blockVector
   for (int i=0; i<BaseMap_.NumMyElements(); i++) {
      localIndex = this->Map().LID((IndexOffset + BaseMap_.GID(i)));
      if (localIndex==-1) { 
	     cout << "Error in  BlockVector::GetBlock: " << i << " " 
		  << IndexOffset << " " << BaseMap_.GID(i) << endl;
	     return -1;
      }
      (*this)[localIndex] = BaseVector[i];
   }

   return 0;
}   
//=========================================================================
int BlockVector::BlockSumIntoGlobalValues(int NumIndices, double* Values,
                                          int* Indices, int GlobalBlockRow)
{
   int IndexOffset = GlobalBlockRow * Offset_;
   int localIndex=0;

   // For each entry in the base vector, translate its global ID
   // by the IndexOffset and load into this blockVector
   for (int i=0; i<NumIndices; i++) {
      localIndex = this->Map().LID((IndexOffset + Indices[i]));
      if (localIndex==-1) { 
	     cout << "Error in  BlockVector::BlockSumIntoGlobalValues: " << i
                  << " " << IndexOffset << " " << Indices[i] << endl;
	     return -1;
      }
      (*this)[localIndex] += Values[i];
   }

   return 0;
}   
//=========================================================================
int BlockVector::BlockReplaceGlobalValues(int NumIndices, double* Values,
                                          int* Indices, int GlobalBlockRow)
{
   int IndexOffset = GlobalBlockRow * Offset_;
   int localIndex=0;

   // For each entry in the base vector, translate its global ID
   // by the IndexOffset and load into this blockVector
   for (int i=0; i<NumIndices; i++) {
      localIndex = this->Map().LID((IndexOffset + Indices[i]));
      if (localIndex==-1) { 
	     cout << "Error in  BlockVector::BlockReplaceGlobalValues: " << i
                  << " " << IndexOffset << " " << Indices[i] << endl;
	     return -1;
      }
      (*this)[localIndex] = Values[i];
   }

   return 0;
}   

//=========================================================================
Teuchos::RCP<const Epetra_Vector>
BlockVector::GetBlock(int GlobalBlockRow) const
{
  int offset = GlobalBlockRow * BaseMap_.NumMyElements();
  return Teuchos::rcp(new Epetra_Vector(View, BaseMap_, Values_+offset));
}

//=========================================================================
Teuchos::RCP<Epetra_Vector>
BlockVector::GetBlock(int GlobalBlockRow)
{
  int offset = GlobalBlockRow * BaseMap_.NumMyElements();
  return Teuchos::rcp(new Epetra_Vector(View, BaseMap_, Values_+offset));
}

//=========================================================================
const Epetra_BlockMap&
BlockVector::GetBaseMap() const
{
  return BaseMap_;
}


} //namespace EpetraExt
