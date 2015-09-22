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

template<typename int_type>
Epetra_Map * BlockUtility::TGenerateBlockMap(
        const Epetra_BlockMap & BaseMap,
        const int_type * RowIndices,
        int NumBlockRows,
        const Epetra_Comm & GlobalComm,
        int_type Offset)
{
  int_type BaseIndex = (int_type) BaseMap.IndexBase64();
  if (Offset == 0)
    Offset = BlockUtility::TCalculateOffset<int_type>(BaseMap);

  //Get Base Global IDs
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int_type> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int_type> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int_type GlobalSize;
  int_type TotalSize_int_type = TotalSize;
  GlobalComm.SumAll( &TotalSize_int_type, &GlobalSize, 1 );

  Epetra_Map *GlobalMap =
    new Epetra_Map( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex,
                    GlobalComm );

  return GlobalMap;
}

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map * BlockUtility::GenerateBlockMap(
        const Epetra_BlockMap & BaseMap,
        const int * RowIndices,
        int NumBlockRows,
        const Epetra_Comm & GlobalComm,
        int Offset)
{
  if(BaseMap.GlobalIndicesInt())
    return TGenerateBlockMap<int>(BaseMap, RowIndices, NumBlockRows, GlobalComm, Offset);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockMap: Global Indices not int.";
}
#endif
//==============================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map * BlockUtility::GenerateBlockMap(
        const Epetra_BlockMap & BaseMap,
        const long long * RowIndices,
        int NumBlockRows,
        const Epetra_Comm & GlobalComm,
        long long Offset)
{
  if(BaseMap.GlobalIndicesLongLong())
    return TGenerateBlockMap<long long>(BaseMap, RowIndices, NumBlockRows, GlobalComm, Offset);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockMap: Global Indices not long long.";
}
#endif
//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map * BlockUtility::GenerateBlockMap(
        const Epetra_BlockMap & BaseMap,
        const std::vector<int>& RowIndices,
        const Epetra_Comm & GlobalComm,
        int Offset )
{
  return GenerateBlockMap(BaseMap, &RowIndices[0], RowIndices.size(),
                          GlobalComm, Offset);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map * BlockUtility::GenerateBlockMap(
        const Epetra_BlockMap & BaseMap,
        const std::vector<long long>& RowIndices,
        const Epetra_Comm & GlobalComm,
        long long Offset )
{
  return GenerateBlockMap(BaseMap, &RowIndices[0], RowIndices.size(),
                          GlobalComm, Offset);
}
#endif

//==============================================================================
Epetra_Map * BlockUtility::GenerateBlockMap(
        const Epetra_BlockMap & BaseMap,
        const Epetra_BlockMap& BlockMap,
        const Epetra_Comm & GlobalComm,
        int Offset)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(BaseMap.GlobalIndicesInt() && BlockMap.GlobalIndicesInt())
    return GenerateBlockMap(BaseMap,
                          BlockMap.MyGlobalElements(),
                          BlockMap.NumMyElements(),
                          GlobalComm,
                          Offset);
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
	  if(BaseMap.GlobalIndicesLongLong() && BlockMap.GlobalIndicesLongLong())
        return GenerateBlockMap(BaseMap,
                          BlockMap.MyGlobalElements64(),
                          BlockMap.NumMyElements(),
                          GlobalComm,
                          Offset);
  else
#endif
    throw "EpetraExt::BlockUtility::GenerateBlockMap: Error Global Indices unknown.";
}

//==============================================================================
template<typename int_type>
Epetra_CrsGraph * BlockUtility::TGenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int_type> > & RowStencil,
        const vector<int_type> & RowIndices,
        const Epetra_Comm & GlobalComm )
{

  const Epetra_BlockMap & BaseMap = BaseGraph.RowMap();
  int_type BaseIndex = (int_type) BaseMap.IndexBase64();
  int_type Offset = BlockUtility::TCalculateOffset<int_type>(BaseMap);

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int_type> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int_type> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int_type GlobalSize;
  int_type TotalSize_int_type = TotalSize;
  GlobalComm.SumAll( &TotalSize_int_type, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int_type> Indices(MaxIndices);
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy,
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    int StencilSize = RowStencil[i].size();
    for( int j = 0; j < Size; ++j )
    {
      int_type BaseRow = (int_type) BaseMap.GID64(j);
      int_type GlobalRow = (int_type) GlobalMap.GID64(j+i*Size);

      BaseGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      for( int k = 0; k < StencilSize; ++k )
      {
        int_type ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm )
{
  if(BaseGraph.RowMap().GlobalIndicesInt())
    return TGenerateBlockGraph<int>(BaseGraph, RowStencil, RowIndices, GlobalComm);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockGraph: Error Global Indices not int.";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<long long> > & RowStencil,
        const vector<long long> & RowIndices,
        const Epetra_Comm & GlobalComm )
{
  if(BaseGraph.RowMap().GlobalIndicesLongLong())
    return TGenerateBlockGraph<long long>(BaseGraph, RowStencil, RowIndices, GlobalComm);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockGraph: Error Global Indices not long long.";
}
#endif

//==============================================================================
template<typename int_type>
Epetra_CrsGraph * BlockUtility::TGenerateBlockGraph(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int_type> > & RowStencil,
        const vector<int_type> & RowIndices,
        const Epetra_Comm & GlobalComm )
{

  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();
  int_type BaseIndex = (int_type) BaseMap.IndexBase64();
  int_type Offset = BlockUtility::TCalculateOffset<int_type>(BaseMap);

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int_type> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int_type> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GlobalGIDs[i*Size+j] = GIDs[j] + RowIndices[i] * Offset;
  }

  int_type GlobalSize;
  int_type TotalSize_int_type = TotalSize;
  GlobalComm.SumAll( &TotalSize_int_type, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices_local(MaxIndices);
  vector<int_type> Indices_global(MaxIndices);
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
      int_type GlobalRow = (int_type) GlobalMap.GID64(j+i*Size);

      BaseMatrix.ExtractMyRowCopy( j, MaxIndices, NumIndices, &Values[0], &Indices_local[0] );
      for( int l = 0; l < NumIndices; ++l ) Indices_global[l] = (int_type) BaseColMap.GID64(Indices_local[l]);

      for( int k = 0; k < StencilSize; ++k )
      {
        int_type ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
        if( k > 0 ) ColOffset -= (RowIndices[i]+RowStencil[i][k-1]) * Offset;

        for( int l = 0; l < NumIndices; ++l )
          Indices_global[l] += ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices_global[0] );
      }
    }
  }

  GlobalGraph->FillComplete();

  return GlobalGraph;
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm )
{
  if(BaseMatrix.RowMatrixRowMap().GlobalIndicesInt())
    return TGenerateBlockGraph<int>(BaseMatrix, RowStencil, RowIndices, GlobalComm);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockGraph: Error Global Indices not int.";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<long long> > & RowStencil,
        const vector<long long> & RowIndices,
        const Epetra_Comm & GlobalComm )
{
  if(BaseMatrix.RowMatrixRowMap().GlobalIndicesLongLong())
    return TGenerateBlockGraph<long long>(BaseMatrix, RowStencil, RowIndices, GlobalComm);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockGraph: Error Global Indices not long long.";
}
#endif

//==============================================================================
template<typename int_type>
Epetra_CrsGraph * BlockUtility::TGenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const Epetra_CrsGraph & LocalBlockGraph,
        const Epetra_Comm & GlobalComm )
{
  const Epetra_BlockMap & BaseRowMap = BaseGraph.RowMap();
  const Epetra_BlockMap & BaseColMap = BaseGraph.ColMap();
  int_type ROffset = BlockUtility::TCalculateOffset<int_type>(BaseRowMap);
  (void) ROffset; // Silence "unused variable" compiler warning.
  int_type COffset = BlockUtility::TCalculateOffset<int_type>(BaseColMap);

  //Get Base Global IDs
  const Epetra_BlockMap & BlockRowMap = LocalBlockGraph.RowMap();
  const Epetra_BlockMap & BlockColMap = LocalBlockGraph.ColMap();

  int NumBlockRows = BlockRowMap.NumMyElements();
  vector<int_type> RowIndices(NumBlockRows);
  BlockRowMap.MyGlobalElements(&RowIndices[0]);

  int Size = BaseRowMap.NumMyElements();

  Epetra_Map *GlobalRowMap =
    GenerateBlockMap(BaseRowMap, BlockRowMap, GlobalComm);


  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int_type> Indices(MaxIndices);

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
      int_type GlobalRow = (int_type) GlobalRowMap->GID64(j+i*Size);

      BaseGraph.ExtractMyRowView( j, NumBaseIndices, BaseIndices );
      for( int k = 0; k < NumBlockIndices; ++k )
      {
        int_type ColOffset = (int_type) BlockColMap.GID64(BlockIndices[k]) * COffset;

        for( int l = 0; l < NumBaseIndices; ++l )
          Indices[l] = (int_type) BaseGraph.GCID64(BaseIndices[l]) + ColOffset;

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

Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const Epetra_CrsGraph & LocalBlockGraph,
        const Epetra_Comm & GlobalComm )
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(BaseGraph.RowMap().GlobalIndicesInt() && LocalBlockGraph.RowMap().GlobalIndicesInt())
    return TGenerateBlockGraph<int>(BaseGraph, LocalBlockGraph, GlobalComm);
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(BaseGraph.RowMap().GlobalIndicesLongLong() && LocalBlockGraph.RowMap().GlobalIndicesLongLong())
      return TGenerateBlockGraph<long long>(BaseGraph, LocalBlockGraph, GlobalComm);
  else
#endif
    throw "EpetraExt::BlockUtility::GenerateBlockGraph: Error Global Indices unknown.";
}

//==============================================================================
template<typename int_type>
void BlockUtility::TGenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph,
                                      std::vector<int_type> RowIndices,
                                      std::vector< std::vector<int_type> >& RowStencil)
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
      int_type Row = RowIndices[i];
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
	  std::vector<int> RowStencil_local(NumCols);
      RowStencil[i].resize(NumCols);
      LocalBlockGraph.ExtractMyRowCopy(i, NumCols, NumCols,
                                       &RowStencil_local[0]);
      for (int k=0; k<NumCols; k++)
        RowStencil[i][k] = (int_type) ((int) (LocalBlockGraph.GCID64(RowStencil_local[k]) - RowIndices[i]));
    }
  }
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
void BlockUtility::GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph,
                                      std::vector<int> RowIndices,
                                      std::vector< std::vector<int> >& RowStencil)
{
  if(LocalBlockGraph.RowMap().GlobalIndicesInt())
    BlockUtility::TGenerateRowStencil<int>(LocalBlockGraph, RowIndices, RowStencil);
  else
    throw "EpetraExt::BlockUtility::GenerateRowStencil: Global Indices not int.";
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
void BlockUtility::GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph,
                                      std::vector<long long> RowIndices,
                                      std::vector< std::vector<long long> >& RowStencil)
{
  if(LocalBlockGraph.RowMap().GlobalIndicesLongLong())
    BlockUtility::TGenerateRowStencil<long long>(LocalBlockGraph, RowIndices, RowStencil);
  else
    throw "EpetraExt::BlockUtility::GenerateRowStencil: Global Indices not long long.";
}
#endif


//==============================================================================
template<typename int_type>
int_type BlockUtility::TCalculateOffset(const Epetra_BlockMap & BaseMap)
{
  int_type MaxGID = (int_type) BaseMap.MaxAllGID64();

//   int Offset = 1;
//   while( Offset <= MaxGID ) Offset *= 10;

//   return Offset;

  return MaxGID+1;
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int BlockUtility::CalculateOffset(const Epetra_BlockMap & BaseMap)
{
  if(BaseMap.GlobalIndicesInt())
    return TCalculateOffset<int>(BaseMap);
  else
    throw "EpetraExt::BlockUtility::GenerateBlockMap: Global Indices not int.";
}
#endif

long long BlockUtility::CalculateOffset64(const Epetra_BlockMap & BaseMap)
{
  return TCalculateOffset<long long>(BaseMap);
}

} //namespace EpetraExt
