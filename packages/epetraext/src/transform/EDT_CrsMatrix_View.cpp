//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

#include <EDT_CrsMatrix_View.h>

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

  newMatrix->TransformToLocal();

  newObj_ = newMatrix;

  return *newMatrix;
}

} // namespace EpetraExt

