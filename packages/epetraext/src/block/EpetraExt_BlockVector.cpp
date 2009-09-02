//@HEADER
// ************************************************************************
// 
//               EpetraExt: Extended Linear Algebra Services Package 
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
