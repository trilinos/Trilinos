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
      const Epetra_BlockMap & GlobalMap,
      int NumBlocks )
  : Epetra_Vector( GlobalMap ),
    BaseMap_( BaseMap ),
    NumBlocks_( NumBlocks ),
    Offset_( BlockUtility::CalculateOffset( BaseMap ) )
{
  AllocateBlocks_();
}

//==========================================================================
// Copy Constructor
BlockVector::BlockVector(const BlockVector& Source)
  : Epetra_Vector( dynamic_cast<const Epetra_Vector &>(Source) ),
    BaseMap_( Source.BaseMap_ ),
    NumBlocks_( Source.NumBlocks_ ),
    Offset_( Source.Offset_ )
{
  AllocateBlocks_();
}

//=========================================================================
BlockVector::~BlockVector()
{
  DeleteBlocks_();
}

//=========================================================================
void BlockVector::AllocateBlocks_(void)
{
  if (BaseMap_.Comm().NumProc() > 1 && NumBlocks_ > 1) {
     if (BaseMap_.Comm().MyPID()==0) 
     cout << "Warning in BlockVector::AllocateBlocks_: This routine does not work\n" 
	  << "\tfor multi-proc base vectors becasue of re-ordering of externals" <<endl;
  }
  
  double * Ptrs;
  ExtractView( &Ptrs );

  Blocks_.resize( NumBlocks_ );
  int NumElements = BaseMap_.NumMyElements();
  for( int i = 0; i < NumBlocks_; ++i )
    Blocks_[i] = new Epetra_Vector( View, BaseMap_, Ptrs+(i*NumElements) );
}

//=========================================================================
void BlockVector::DeleteBlocks_(void)
{
  for( int i = 0; i < NumBlocks_; ++i )
  {
    delete Blocks_[i];
    Blocks_[i] = 0;
  }
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
int BlockVector::LoadBlockValues(Epetra_Vector & BaseVector, int GlobalBlockRow) 
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


} //namespace EpetraExt
