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

#include <EpetraExt_Reindex_CrsMatrix.h>

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_IntVector.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace EpetraExt {

CrsMatrix_Reindex::
~CrsMatrix_Reindex()
{
  if( newObj_ ) delete newObj_;
  if( NewColMap_ ) delete NewColMap_;
}

CrsMatrix_Reindex::NewTypeRef
CrsMatrix_Reindex::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //test map, must have same number of local and global elements as original row map
  Epetra_Map & OldRowMap = const_cast<Epetra_Map&>(orig.RowMap());
  Epetra_Map & OldColMap = const_cast<Epetra_Map&>(orig.ColMap());
  int NumMyElements = OldRowMap.NumMyElements();
  assert( OldRowMap.NumMyElements() == NewRowMap_.NumMyElements() );

  //Construct new Column Map
  Epetra_IntVector Cols( OldRowMap );
  Epetra_IntVector NewCols( OldColMap );
  Epetra_Import Importer( OldColMap, OldRowMap );

  for( int i = 0; i < NumMyElements; ++i )
    Cols[i] = NewRowMap_.GID(i);

  NewCols.Import( Cols, Importer, Insert );

  vector<int*> NewColIndices(1);
  NewCols.ExtractView( &NewColIndices[0] );

  int NumMyColElements = OldColMap.NumMyElements();
  int NumGlobalColElements = OldColMap.NumGlobalElements();

  NewColMap_ = new Epetra_Map( NumGlobalColElements, NumMyColElements, NewColIndices[0], OldColMap.IndexBase(), OldColMap.Comm() );

  //intial construction of matrix 
  Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, NewRowMap_, *NewColMap_, 0 );

  //insert views of row values
  int * myIndices;
  double * myValues;
  int indicesCnt;
  int numMyRows = NewMatrix->NumMyRows();
  for( int i = 0; i < numMyRows; ++i )
  {
    orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
    NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
  }

  NewMatrix->FillComplete();

  newObj_ = NewMatrix;

  return *NewMatrix;
}

} // namespace EpetraExt

