//@HEADER
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
//@HEADER

#include <EpetraExt_MapColoring.h>

#include <EpetraExt_Transpose_CrsGraph.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_Map.h>

#ifdef EPETRAEXT_TIMING
#include <Epetra_Time.h>
#endif

#include <vector>
#include <set>
#include <map>

using std::vector;
using std::set;
using std::map;

namespace EpetraExt {

CrsGraph_MapColoring::NewTypeRef
CrsGraph_MapColoring::
operator()( OriginalTypeRef orig  )
{
#ifdef EPETRAEXT_TIMING
  Epetra_Time timer( orig.Comm() );
#endif

  origObj_ = &orig;

  int err;

  const Epetra_BlockMap & RowMap = orig.RowMap();
  int nRows = RowMap.NumMyElements();
  const Epetra_BlockMap & ColMap = orig.ColMap();
  int nCols = ColMap.NumMyElements();

  if( verbose_ ) cout << "RowMap:\n" << RowMap;
  if( verbose_ ) cout << "ColMap:\n" << ColMap;

  Epetra_CrsGraph * base = &( const_cast<Epetra_CrsGraph&>(orig) );

  int NumIndices;
  int * Indices;

#ifdef EPETRAEXT_TIMING
  double wTime1 = timer.WallTime();
#endif

  // For parallel applications, add in boundaries to coloring
  bool distributedGraph = RowMap.DistributedGlobal();
  if( distributedGraph )
  {
    base = new Epetra_CrsGraph( Copy, ColMap, ColMap, 0 );

    for( int i = 0; i < nRows; ++i )
    {
      assert( orig.ExtractMyRowView( i, NumIndices, Indices ) == 0 );
      assert( base->InsertMyIndices( i, NumIndices, Indices ) >= 0 );

      //Do this with a single insert
      //Is this the right thing?
      for( int j = 0; j < NumIndices; ++j )
        if( Indices[j] >= nRows )
          assert( base->InsertMyIndices( Indices[j], 1, &i ) >= 0 );
    } 
    base->TransformToLocal();
  }

  if( verbose_ ) cout << "Base Graph:\n" << *base << endl;

#ifdef EPETRAEXT_TIMING
  double wTime2 = timer.WallTime();
  cout << "EpetraExt::MapColoring [INSERT BOUNDARIES] Time: " << wTime2-wTime1 << endl;
#endif

  //Generate Local Distance-1 Adjacency Graph
  //(Transpose of orig excluding non-local column indices)
  EpetraExt::CrsGraph_Transpose transposeTransform( true );
  Epetra_CrsGraph & Adj1 = transposeTransform( *base );
  if( verbose_ ) cout << "Adjacency 1 Graph!\n" << Adj1;

#ifdef EPETRAEXT_TIMING
  wTime1 = timer.WallTime();
  cout << "EpetraExt::MapColoring [TRANSPOSE GRAPH] Time: " << wTime1-wTime2 << endl;
#endif

  int Delta = Adj1.MaxNumIndices();
  if( verbose_ ) cout << endl << "Delta: " << Delta << endl;

  //Generation of Local Distance-2 Adjacency Graph
  Epetra_CrsGraph Adj2( Copy, ColMap, ColMap, 0 );
  int NumAdj1Indices;
  int * Adj1Indices;
  for( int i = 0; i < nCols; ++i )
  {
    assert( Adj1.ExtractMyRowView( i, NumAdj1Indices, Adj1Indices ) == 0 );

    for( int j = 0; j < NumAdj1Indices; ++j )
    {
      assert( base->ExtractMyRowView( Adj1Indices[j], NumIndices, Indices ) == 0 );
      int NumLocalIndices = 0;
      for( int k = 0; k < NumIndices; ++k )
        if( Indices[k] < nCols ) NumLocalIndices++; 
      assert( Adj2.InsertMyIndices( i, NumLocalIndices, Indices ) >= 0 );
    }
  }
  assert( Adj2.TransformToLocal() == 0 );
  if( verbose_ ) cout << "Adjacency 2 Graph!\n" << Adj2;

#ifdef EPETRAEXT_TIMING
  wTime2 = timer.WallTime();
  cout << "EpetraExt::MapColoring [GEN DIST-2 GRAPH] Time: " << wTime2-wTime1 << endl;
#endif

  vector<int> rowOrder( nCols );
  //Simple reordering
  {
    multimap<int,int> adjMap;
    for( int i = 0; i < nCols; ++i )
#ifdef MUST_CONST_STL_MAP_KEY
      adjMap.insert( pair<const int,int>( Adj2.NumMyIndices(i), i ) );
#else
      adjMap.insert( pair<int,int>( Adj2.NumMyIndices(i), i ) );
#endif
    multimap<int,int>::iterator iter = adjMap.begin();
    multimap<int,int>::iterator end = adjMap.end();
    for( int i = 1; iter != end; ++iter, ++i )
      rowOrder[ nCols - i ] = (*iter).second;
  }

#ifdef EPETRAEXT_TIMING
  wTime1 = timer.WallTime();
  cout << "EpetraExt::MapColoring [REORDERING] Time: " << wTime1-wTime2 << endl;
#endif

  Epetra_MapColoring * ColorMap = new Epetra_MapColoring( ColMap );

  //Application of Greedy Algorithm to generate Color Map
  int Size = Delta * Delta + 1;
  set<int> allowedColors;
  for( int col = 0; col < nCols; ++col )
  {
    for( int i = 0; i < Size; ++i ) allowedColors.insert( i+1 ); 

    Adj2.ExtractMyRowView( rowOrder[col], NumIndices, Indices );

    for( int i = 0; i < NumIndices; ++i )
      if( (*ColorMap)[ Indices[i] ] > 0 ) allowedColors.erase( (*ColorMap)[ Indices[i] ] );

    (*ColorMap)[ rowOrder[col] ] = *(allowedColors.begin());
  }

#ifdef EPETRAEXT_TIMING
  wTime2 = timer.WallTime();
  cout << "EpetraExt::MapColoring [GREEDY COLORING] Time: " << wTime2-wTime1 << endl;
#endif

  if( verbose_ ) cout << "ColorMap!\n" << *ColorMap;

  if( distributedGraph ) delete base;

  newObj_ = ColorMap;

  return *ColorMap;
}

} // namespace EpetraExt
