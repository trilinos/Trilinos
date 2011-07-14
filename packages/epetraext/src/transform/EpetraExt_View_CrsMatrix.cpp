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

#include <EpetraExt_View_CrsMatrix.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

namespace EpetraExt {

CrsMatrix_View::
~CrsMatrix_View()
{
  if( newObj_ ) delete newObj_;
}

CrsMatrix_View::NewTypeRef
CrsMatrix_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  if( orig.IndicesAreGlobal() ) cout << "EDT_CrsMatrix_View: Indices must be LOCAL!\n";
  assert( !orig.IndicesAreGlobal() );

  //test graph, new graph must be contiguous subset of old

  //intial construction of matrix 
  Epetra_CrsMatrix * newMatrix( new Epetra_CrsMatrix( View, NewGraph_ ) );

  //insert views of row values
  int * myIndices;
  double * myValues;
  int indicesCnt;
  int numMyRows = newMatrix->NumMyRows();
  for( int i = 0; i < numMyRows; ++i )
  {
    orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );

    int newIndicesCnt = indicesCnt;
    bool done = false;
    for( int j = 0; j < indicesCnt; ++j )
      if( !done && NewGraph_.GCID( myIndices[j] ) == -1 )
      {
        newIndicesCnt = j;
        done = true;
      }

    newMatrix->InsertMyValues( i, newIndicesCnt, myValues, myIndices );
  }

  newMatrix->FillComplete(false);

  newObj_ = newMatrix;

  return *newMatrix;
}

} // namespace EpetraExt

