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
#include "Epetra_Map.h"

namespace EpetraExt {

//=============================================================================
// EpetraExt::BlockVector Constructor
BlockVector::BlockVector(
      const Epetra_BlockMap & BaseMap,
      const Epetra_BlockMap & GlobalMap,
      int NumBlocks )
  : Epetra_Vector( GlobalMap ),
    BaseMap_( BaseMap ),
    NumBlocks_( NumBlocks )
{
  AllocateBlocks_();
}

//==========================================================================
// Copy Constructor
BlockVector::BlockVector(const BlockVector& Source)
  : Epetra_Vector( dynamic_cast<const Epetra_Vector &>(Source) ),
    BaseMap_( Source.BaseMap_ ),
    NumBlocks_( Source.NumBlocks_ )
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

} //namespace EpetraExt
