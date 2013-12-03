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

#include <Epetra_ConfigDefs.h>
#include <EpetraExt_Reindex_CrsMatrix.h>

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_GIDTypeVector.h>

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
#include <Epetra_LongLongVector.h>
#endif

#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace EpetraExt {

CrsMatrix_Reindex::
~CrsMatrix_Reindex()
{
  if( newObj_ ) delete newObj_;
  if( NewColMap_ ) delete NewColMap_;
}

template<typename int_type>
CrsMatrix_Reindex::NewTypeRef
CrsMatrix_Reindex::
Toperator( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //test std::map, must have same number of local and global elements as original row std::map
  //Epetra_Map & OldRowMap = const_cast<Epetra_Map&>(orig.RowMap());
  Epetra_Map & OldDomainMap = const_cast<Epetra_Map&>(orig.OperatorDomainMap());
  Epetra_Map & OldColMap = const_cast<Epetra_Map&>(orig.ColMap());
  int NumMyElements = OldDomainMap.NumMyElements();
  int_type NumGlobalElements = (int_type) OldDomainMap.NumGlobalElements64();
  assert( orig.RowMap().NumMyElements() == NewRowMap_.NumMyElements() );

  if (NumGlobalElements == 0 && orig.RowMap().NumGlobalElements64() == 0 )
  {
    //construct a zero matrix as a placeholder, don't do reindexing analysis.
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, orig.RowMap(), orig.ColMap(), 0 );
    newObj_ = NewMatrix;
  }
  else {

    //Construct new Column Map
    typename Epetra_GIDTypeVector<int_type>::impl Cols( OldDomainMap );
    typename Epetra_GIDTypeVector<int_type>::impl NewCols( OldColMap );
    Epetra_Import Importer( OldColMap, OldDomainMap );
 
    Epetra_Map tmpColMap( NumGlobalElements, NumMyElements, 0, OldDomainMap.Comm() );
 
    for( int i = 0; i < NumMyElements; ++i )
      Cols[i] = (int_type) tmpColMap.GID64(i);

    NewCols.Import( Cols, Importer, Insert );

    std::vector<int_type*> NewColIndices(1);
    NewCols.ExtractView( &NewColIndices[0] );

    int NumMyColElements = OldColMap.NumMyElements();
    int_type NumGlobalColElements = (int_type) OldColMap.NumGlobalElements64();

    NewColMap_ = new Epetra_Map( NumGlobalColElements, NumMyColElements, NewColIndices[0], (int_type) OldColMap.IndexBase64(), OldColMap.Comm() );

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

CrsMatrix_Reindex::NewTypeRef
CrsMatrix_Reindex::
operator()( OriginalTypeRef orig )
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(orig.RowMatrixRowMap().GlobalIndicesInt())
    return Toperator<int>(orig);
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(orig.RowMatrixRowMap().GlobalIndicesLongLong())
    return Toperator<long long>(orig);
  else
#endif
    throw "EpetraExt::CrsMatrix_Reindex::operator(): Global indices unknown.";
}

} // namespace EpetraExt

