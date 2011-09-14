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

#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

namespace EpetraExt {

using std::vector;

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
	const Epetra_BlockMap & BaseMap,
        const int * RowIndices,
	int NumBlockRows,
        const Epetra_Comm & GlobalComm ) 
{
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map *GlobalMap = 
    new Epetra_Map( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, 
		    GlobalComm );

  return GlobalMap;
}

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
	const Epetra_BlockMap & BaseMap,
        const std::vector<int>& RowIndices,
	const Epetra_Comm & GlobalComm ) 
{
  return GenerateBlockMap(BaseMap, &RowIndices[0], RowIndices.size(), 
			  GlobalComm);
}

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
	const Epetra_BlockMap & BaseMap,
        const Epetra_BlockMap& BlockMap,
	const Epetra_Comm & GlobalComm ) 
{
  return GenerateBlockMap(BaseMap, 
			  BlockMap.MyGlobalElements(), 
			  BlockMap.NumMyElements(), 
			  GlobalComm);
}

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm ) 
{

  const Epetra_BlockMap & BaseMap = BaseGraph.RowMap();
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int> Indices(MaxIndices);
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    int StencilSize = RowStencil[i].size();
    for( int j = 0; j < Size; ++j )
    {
      int BaseRow = BaseMap.GID(j);
      int GlobalRow = GlobalMap.GID(j+i*Size);

      BaseGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      for( int k = 0; k < StencilSize; ++k )
      {
        int ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
        if( k > 0 ) ColOffset -= (RowIndices[i]+RowStencil[i][k-1]) * Offset;

        for( int l = 0; l < NumIndices; ++l )
          Indices[l] += ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices[0] );
      }
    }
  }

  GlobalGraph->FillComplete();

  return GlobalGraph;
}

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm ) 
{

  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();
  int BaseIndex = BaseMap.IndexBase();
  int Offset = BlockUtility::CalculateOffset(BaseMap);

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices(MaxIndices);
  vector<double> Values(MaxIndices); 
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    int StencilSize = RowStencil[i].size();
    for( int j = 0; j < Size; ++j )
    {
      int GlobalRow = GlobalMap.GID(j+i*Size);

      BaseMatrix.ExtractMyRowCopy( j, MaxIndices, NumIndices, &Values[0], &Indices[0] );
      for( int l = 0; l < NumIndices; ++l ) Indices[l] = BaseColMap.GID(Indices[l]);

      for( int k = 0; k < StencilSize; ++k )
      {
        int ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
        if( k > 0 ) ColOffset -= (RowIndices[i]+RowStencil[i][k-1]) * Offset;

        for( int l = 0; l < NumIndices; ++l )
          Indices[l] += ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices[0] );
      }
    }
  }

  GlobalGraph->FillComplete();

  return GlobalGraph;
}

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const Epetra_CrsGraph & LocalBlockGraph,
        const Epetra_Comm & GlobalComm ) 
{
  const Epetra_BlockMap & BaseRowMap = BaseGraph.RowMap();
  const Epetra_BlockMap & BaseColMap = BaseGraph.ColMap();
  int ROffset = BlockUtility::CalculateOffset(BaseRowMap);
  int COffset = BlockUtility::CalculateOffset(BaseColMap);

  //Get Base Global IDs
  const Epetra_BlockMap & BlockRowMap = LocalBlockGraph.RowMap();
  const Epetra_BlockMap & BlockColMap = LocalBlockGraph.ColMap();
  
  int NumBlockRows = BlockRowMap.NumMyElements();
  vector<int> RowIndices(NumBlockRows);
  BlockRowMap.MyGlobalElements(&RowIndices[0]);

  int Size = BaseRowMap.NumMyElements();
  
  Epetra_Map *GlobalRowMap = 
    GenerateBlockMap(BaseRowMap, BlockRowMap, GlobalComm);
  

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int> Indices(MaxIndices);

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(*GlobalRowMap),
                               0 );

  int NumBlockIndices, NumBaseIndices;
  int *BlockIndices, *BaseIndices;
  for( int i = 0; i < NumBlockRows; ++i )
  {
    LocalBlockGraph.ExtractMyRowView(i, NumBlockIndices, BlockIndices);
    
    for( int j = 0; j < Size; ++j )
    {
      int GlobalRow = GlobalRowMap->GID(j+i*Size);

      BaseGraph.ExtractMyRowView( j, NumBaseIndices, BaseIndices );
      for( int k = 0; k < NumBlockIndices; ++k )
      {
        int ColOffset = BlockColMap.GID(BlockIndices[k]) * COffset;

        for( int l = 0; l < NumBaseIndices; ++l )
          Indices[l] = BaseGraph.GCID(BaseIndices[l]) + ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumBaseIndices, &Indices[0] );
      }
    }
  }

  const Epetra_BlockMap & BaseDomainMap = BaseGraph.DomainMap();
  const Epetra_BlockMap & BaseRangeMap = BaseGraph.RangeMap();
  const Epetra_BlockMap & BlockDomainMap = LocalBlockGraph.DomainMap();
  const Epetra_BlockMap & BlockRangeMap = LocalBlockGraph.RangeMap();

  Epetra_Map *GlobalDomainMap = 
    GenerateBlockMap(BaseDomainMap, BlockDomainMap, GlobalComm);
  Epetra_Map *GlobalRangeMap = 
    GenerateBlockMap(BaseRangeMap, BlockRangeMap, GlobalComm);

  GlobalGraph->FillComplete(*GlobalDomainMap, *GlobalRangeMap);

  delete GlobalDomainMap;
  delete GlobalRangeMap;
  delete GlobalRowMap;

  return GlobalGraph;
}

//==============================================================================
void BlockUtility::GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph, 
				      std::vector<int> RowIndices, 
				      std::vector< std::vector<int> >& RowStencil)
{
  // Get row indices
  int NumMyRows = LocalBlockGraph.NumMyRows();
  RowIndices.resize(NumMyRows);
  const Epetra_BlockMap& RowMap = LocalBlockGraph.RowMap();
  RowMap.MyGlobalElements(&RowIndices[0]);

  // Get stencil
  RowStencil.resize(NumMyRows);
  if (LocalBlockGraph.IndicesAreGlobal()) {
    for (int i=0; i<NumMyRows; i++) {
      int Row = RowIndices[i];
      int NumCols = LocalBlockGraph.NumGlobalIndices(Row);
      RowStencil[i].resize(NumCols);
      LocalBlockGraph.ExtractGlobalRowCopy(Row, NumCols, NumCols, 
					   &RowStencil[i][0]);
      for (int k=0; k<NumCols; k++)
	RowStencil[i][k] -= Row;
    }
  }
  else {
    for (int i=0; i<NumMyRows; i++) {
      int NumCols = LocalBlockGraph.NumMyIndices(i);
      RowStencil[i].resize(NumCols);
      LocalBlockGraph.ExtractMyRowCopy(i, NumCols, NumCols, 
				       &RowStencil[i][0]);
      for (int k=0; k<NumCols; k++)
	RowStencil[i][k] = LocalBlockGraph.GCID(RowStencil[i][k]) - 
	  RowIndices[i];
    }
  }
}

//==============================================================================
int BlockUtility::CalculateOffset(const Epetra_BlockMap & BaseMap)
{
  int MaxGID = BaseMap.MaxAllGID();

//   int Offset = 1;
//   while( Offset <= MaxGID ) Offset *= 10;

//   return Offset;

  return MaxGID+1;
}


} //namespace EpetraExt
