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
    NumBlocks_(1)
{
  AllocateBlocks_();
}

//=============================================================================
// EpetraExt::BlockMultiVector Constructor
BlockMultiVector::BlockMultiVector(
      const Epetra_BlockMap & BaseMap,
      const Epetra_BlockMap & GlobalMap,
      int NumBlocks,
      int NumVectors )
  : Epetra_MultiVector( GlobalMap, NumVectors ),
    BaseMap_( BaseMap ),
    NumBlocks_(NumBlocks)
{
  AllocateBlocks_();
}

//==========================================================================
// Copy Constructor
BlockMultiVector::BlockMultiVector(const BlockMultiVector& Source)
  : Epetra_MultiVector( dynamic_cast<const Epetra_MultiVector &>(Source) ),
    BaseMap_( Source.BaseMap_ ),
    NumBlocks_( Source.NumBlocks_ )
{
  AllocateBlocks_();
}

//=========================================================================
BlockMultiVector::~BlockMultiVector()
{
  DeleteBlocks_();
}

//=========================================================================
void BlockMultiVector::AllocateBlocks_(void)
{
  int NumElements = BaseMap_.NumMyElements();
  Ptrs_.resize(NumBlocks_);
  for( int i = 0; i < NumBlocks_; ++i ) Ptrs_[i] = new double*[NumVectors()];

  double ** OrigPtrs;
  ExtractView( &OrigPtrs );

  for( int i = 0; i < NumBlocks_; ++i )
  {
    for( int j = 0; j < NumVectors(); ++j )
      Ptrs_[i][j] = OrigPtrs[j]+(i*NumElements);

    Blocks_[i] = new Epetra_MultiVector( View, BaseMap_, Ptrs_[i], NumVectors() );
  }
}

//=========================================================================
void BlockMultiVector::DeleteBlocks_(void)
{
  for( int i = 0; i < NumBlocks_; ++i )
  {
    delete Blocks_[i];
    Blocks_[i] = 0;
    delete [] Ptrs_[i];
    Ptrs_[i] = 0;
  }
}

} //namespace EpetraExt
