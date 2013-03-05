//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

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

  //test std::map, must have same number of local and global elements as original row std::map
  //Epetra_Map & OldRowMap = const_cast<Epetra_Map&>(orig.RowMap());
  Epetra_Map & OldDomainMap = const_cast<Epetra_Map&>(orig.OperatorDomainMap());
  Epetra_Map & OldColMap = const_cast<Epetra_Map&>(orig.ColMap());
  int NumMyElements = OldDomainMap.NumMyElements();
  int NumGlobalElements = OldDomainMap.NumGlobalElements();
  assert( orig.RowMap().NumMyElements() == NewRowMap_.NumMyElements() );

  if (NumGlobalElements == 0 && orig.RowMap().NumGlobalElements() == 0 )
  {
    //construct a zero matrix as a placeholder, don't do reindexing analysis.
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, orig.RowMap(), orig.ColMap(), 0 );
    newObj_ = NewMatrix;
  }
  else {

    //Construct new Column Map
    Epetra_IntVector Cols( OldDomainMap );
    Epetra_IntVector NewCols( OldColMap );
    Epetra_Import Importer( OldColMap, OldDomainMap );
 
    Epetra_Map tmpColMap( NumGlobalElements, NumMyElements, 0, OldDomainMap.Comm() );
 
    for( int i = 0; i < NumMyElements; ++i )
      Cols[i] = tmpColMap.GID(i);

    NewCols.Import( Cols, Importer, Insert );

    std::vector<int*> NewColIndices(1);
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

  }

  return *newObj_;
}

} // namespace EpetraExt

