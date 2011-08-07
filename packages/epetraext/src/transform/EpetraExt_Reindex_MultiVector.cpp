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

#include <EpetraExt_Reindex_MultiVector.h>

#include <vector>

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>

namespace EpetraExt {

MultiVector_Reindex::
~MultiVector_Reindex()
{
  if( newObj_ ) delete newObj_;
}

MultiVector_Reindex::NewTypeRef
MultiVector_Reindex::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //test std::map, must have same number of local and global elements as original row std::map
  assert( orig.Map().NumMyElements() == NewRowMap_.NumMyElements() );

  std::vector<double*> MyValues(1);
  int MyLDA;
  int NumVectors = orig.NumVectors();
  orig.ExtractView( &MyValues[0], &MyLDA );

  Epetra_MultiVector * NewMV = new Epetra_MultiVector( View, NewRowMap_, MyValues[0], MyLDA, NumVectors );

  newObj_ = NewMV;

  return *NewMV;
}

MultiVector_Reindex::NewTypeRCP
MultiVector_Reindex::transform(OriginalTypeRef orig)
{
  //test std::map, must have same number of local and global elements as original row std::map
  assert( orig.Map().NumMyElements() == NewRowMap_.NumMyElements() );

  std::vector<double*> MyValues(1);
  int MyLDA;
  int NumVectors = orig.NumVectors();
  orig.ExtractView( &MyValues[0], &MyLDA );

  return Teuchos::rcp(new Epetra_MultiVector( View, NewRowMap_, MyValues[0], MyLDA, NumVectors ));
}

} // namespace EpetraExt

