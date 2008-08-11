// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// ***********************************************************************
// @HEADER

#include <EpetraExt_View_MultiVector.h>

#include <Epetra_MultiVector.h>
#include <Epetra_BlockMap.h>

namespace EpetraExt {

MultiVector_View::
~MultiVector_View()
{
  if( newObj_ ) delete newObj_;
}

MultiVector_View::NewTypeRef
MultiVector_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int numVec = NumVec_;
  if( numVec == -1 ) numVec = orig.NumVectors();

  double ** ptrArray;
  orig.ExtractView( &ptrArray );

  Epetra_MultiVector * newMV = new Epetra_MultiVector( View, NewMap_, ptrArray, numVec );

  newObj_ = newMV;

  return *newMV;
}

} // namespace EpetraExt

