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

} //namespace EpetraExt
