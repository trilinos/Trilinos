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
#include <EpetraExt_Overlap_CrsGraph.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_IntVector.h>
#include <Epetra_MapColoring.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_Util.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>

#ifdef EPETRAEXT_TIMING
#include <Epetra_Time.h>
#endif

#include <vector>
#include <set>
#include <map>

using std::vector;
using std::set;
using std::map;

#define PRAND(seed) (fmod((double)23*(47594118*seed),(double)10000001))

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

  Epetra_MapColoring * ColorMap = 0;

  if( !colorParallel_ )
  {

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

    set<int> Cols;
    for( int j = 0; j < NumAdj1Indices; ++j )
    {
      assert( base->ExtractMyRowView( Adj1Indices[j], NumIndices, Indices ) == 0 );
      for( int k = 0; k < NumIndices; ++k )
        if( Indices[k] < nCols ) Cols.insert( Indices[k] );
    }
    int nCols2 = Cols.size();
    vector<int> ColVec( nCols2 );
    set<int>::iterator iterIS = Cols.begin();
    set<int>::iterator iendIS = Cols.end();
    for( int j = 0 ; iterIS != iendIS; ++iterIS, ++j ) ColVec[j] = *iterIS;
    assert( Adj2.InsertMyIndices( i, nCols2, &ColVec[0] ) >= 0 );
  }
  assert( Adj2.TransformToLocal() == 0 );
  if( verbose_ ) cout << "Adjacency 2 Graph!\n" << Adj2;

#ifdef EPETRAEXT_TIMING
  wTime2 = timer.WallTime();
  cout << "EpetraExt::MapColoring [GEN DIST-2 GRAPH] Time: " << wTime2-wTime1 << endl;
#endif

  ColorMap = new Epetra_MapColoring( ColMap );

  vector<int> rowOrder( nCols );
  if( reordering_ == 0 || reordering_ == 1 ) 
  {
    multimap<int,int> adjMap;
    typedef multimap<int,int>::value_type adjMapValueType;
    for( int i = 0; i < nCols; ++i )
      adjMap.insert( adjMapValueType( Adj2.NumMyIndices(i), i ) );
    multimap<int,int>::iterator iter = adjMap.begin();
    multimap<int,int>::iterator end = adjMap.end();
    if( reordering_ == 0 ) //largest first (less colors)
    {
      for( int i = 1; iter != end; ++iter, ++i )
        rowOrder[ nCols - i ] = (*iter).second;
    }
    else                  //smallest first (better balance)
    {
      for( int i = 0; iter != end; ++iter, ++i )
        rowOrder[ i ] = (*iter).second;
    }
  }
  else if( reordering_ == 2 ) //random
  {
    for( int i = 0; i < nCols; ++i )
      rowOrder[ i ] = i;
#ifndef TFLOP
    random_shuffle( rowOrder.begin(), rowOrder.end() );
#endif
  }

#ifdef EPETRAEXT_TIMING
  wTime1 = timer.WallTime();
  cout << "EpetraExt::MapColoring [REORDERING] Time: " << wTime1-wTime2 << endl;
#endif

  if( algo_ == ALGO_GREEDY )
  {
    //Application of Greedy Algorithm to generate Color Map
    for( int col = 0; col < nCols; ++col )
    {
      Adj2.ExtractMyRowView( rowOrder[col], NumIndices, Indices );

      set<int> usedColors;
      int color;
      for( int i = 0; i < NumIndices; ++i )
      {
        color = (*ColorMap)[ Indices[i] ];
        if( color > 0 ) usedColors.insert( color );
        color = 0;
        int testcolor = 1;
        while( !color )
        {
          if( !usedColors.count( testcolor ) ) color = testcolor;
          ++testcolor;
        }
      }
      (*ColorMap)[ rowOrder[col] ] = color;
    }

#ifdef EPETRAEXT_TIMING
    wTime2 = timer.WallTime();
    cout << "EpetraExt::MapColoring [GREEDY COLORING] Time: " << wTime2-wTime1 << endl;
    cout << "Num GREEDY Colors: " << ColorMap->NumColors() << endl;
#endif
  }
  else if( algo_ == ALGO_LUBI )
  {
    //Assign Random Keys To Rows
    Epetra_Util util;
    vector<int> Keys(nCols);
    vector<int> State(nCols,-1);

    for( int col = 0; col < nCols; ++col )
      Keys[col] = util.RandomInt();

    int NumRemaining = nCols;
    int CurrentColor = 1;

    while( NumRemaining > 0 )
    {
      //maximal independent set
      while( NumRemaining > 0 )
      {
        NumRemaining = 0;

        //zero out everyone less than neighbor
        for( int i = 0; i < nCols; ++i )
        {
          int col = rowOrder[i];
          if( State[col] < 0 )
          {
#ifdef EPETRAEXT_TEST_LUBI
            cout << "Testing Node: " << col << " " << State[col] << " " << Keys[col] << endl;
#endif
            Adj2.ExtractMyRowView( col, NumIndices, Indices );
            int MyKey = Keys[col];
            for( int j = 0; j < NumIndices; ++j )
              if( col != Indices[j] && State[ Indices[j] ] < 0 )
              {
                if( MyKey > Keys[ Indices[j] ] ) State[ Indices[j] ] = 0;
                else                             State[ col ] = 0;
              }
          }
        }

        //assign -1's to current color
        for( int col = 0; col < nCols; ++col )
        {
#ifdef EPETRAEXT_TEST_LUBI
          cout << "Node: " << col << " State: " << State[col] << endl;
#endif
          if( State[col] < 0 )
          {
            State[col] = CurrentColor;
#ifdef EPETRAEXT_TEST_LUBI
            cout << "Node,Color: " << col << "," << CurrentColor << endl;
#endif
	  }
        }

        //reinstate any zero not neighboring current color
        for( int col = 0; col < nCols; ++col )
          if( State[col] == 0 )
          {
            Adj2.ExtractMyRowView( col, NumIndices, Indices );
	    bool flag = false;
            for( int i = 0; i < NumIndices; ++i )
              if( col != Indices[i] && State[ Indices[i] ] == CurrentColor )
              {
                flag = true;
                break;
              }
            if( !flag )
            {
              State[col] = -1;
              ++NumRemaining;
            }
          }
      }

      //Reset Status for all non-colored nodes
      for( int col = 0; col < nCols; ++col )
        if( State[col] == 0 )
        {
          State[col] = -1;
          ++NumRemaining;
        }

      if( verbose_ )
      {
        cout << "Finished Color: " << CurrentColor << endl;
        cout << "NumRemaining: " << NumRemaining << endl;
      }

      //New color
      ++CurrentColor;
    }

    for( int col = 0; col < nCols; ++col )
      (*ColorMap)[col] = State[col]-1;

#ifdef EPETRAEXT_TIMING
    wTime2 = timer.WallTime();
    cout << "EpetraExt::MapColoring [LUBI COLORING] Time: " << wTime2-wTime1 << endl;
    cout << "Num LUBI Colors: " << ColorMap->NumColors() << endl;
#endif

    if( distributedGraph ) delete base;
  }
  else
    abort(); //UNKNOWN ALGORITHM

  }
  else
  {
    //Generate Parallel Adjacency-2 Graph

    EpetraExt::CrsGraph_Overlap OverlapTrans(1);
    Epetra_CrsGraph & OverlapGraph = OverlapTrans( orig );
    if( verbose_ ) cout << "OverlapGraph:\n" << OverlapGraph;

    int MyPID = RowMap.Comm().MyPID();

    Epetra_CrsGraph Adj2( Copy, RowMap, OverlapGraph.ColMap(), 0 );
    int NumIndices;
    int * Indices;
    int NumAdj1Indices;
    int * Adj1Indices;
    for( int i = 0; i < nRows; ++i )
    {
      assert( OverlapGraph.ExtractMyRowView( i, NumAdj1Indices, Adj1Indices ) == 0 );

      set<int> Cols;
      for( int j = 0; j < NumAdj1Indices; ++j )
      {
        int GID = OverlapGraph.LRID( OverlapGraph.GCID( Adj1Indices[j] ) );
        assert( OverlapGraph.ExtractMyRowView( GID, NumIndices, Indices ) == 0 );
        for( int k = 0; k < NumIndices; ++k ) Cols.insert( Indices[k] );
      }
      int nCols2 = Cols.size();
      vector<int> ColVec( nCols2 );
      set<int>::iterator iterIS = Cols.begin();
      set<int>::iterator iendIS = Cols.end();
      for( int j = 0 ; iterIS != iendIS; ++iterIS, ++j ) ColVec[j] = *iterIS;
      assert( Adj2.InsertMyIndices( i, nCols2, &ColVec[0] ) >= 0 );
    }
    assert( Adj2.TransformToLocal() == 0 );
    RowMap.Comm().Barrier();
    if( verbose_ ) cout << "Adjacency 2 Graph!\n" << Adj2;

    //collect GIDs on boundary
    vector<int> boundaryGIDs;
    vector<int> interiorGIDs;
    for( int row = 0; row < nRows; ++row )
    {
      Adj2.ExtractMyRowView( row, NumIndices, Indices );
      bool testFlag = false;
      for( int i = 0; i < NumIndices; ++i )
        if( Indices[i] >= nRows ) testFlag = true;
      if( testFlag ) boundaryGIDs.push_back( Adj2.GRID(row) );
      else           interiorGIDs.push_back( Adj2.GRID(row) );
    }

    int LocalBoundarySize = boundaryGIDs.size();
    int LocalInteriorSize = interiorGIDs.size();

    Epetra_Map BoundaryMap( -1, boundaryGIDs.size(), &boundaryGIDs[0], 0, RowMap.Comm() );
    if( verbose_ ) cout << "BoundaryMap:\n" << BoundaryMap;
    
    int BoundarySize = BoundaryMap.NumGlobalElements();
    Epetra_MapColoring BoundaryColoring( BoundaryMap );

    if( SerialBoundaryColoring_ )
    {
      Epetra_Map BoundaryIndexMap( BoundarySize, LocalBoundarySize, 0, RowMap.Comm() );
      if( verbose_ ) cout << "BoundaryIndexMap:\n" << BoundaryIndexMap;

      Epetra_IntVector bGIDs( View, BoundaryIndexMap, &boundaryGIDs[0] );
      if( verbose_ ) cout << "BoundaryGIDs:\n" << bGIDs;

      int NumLocalBs = 0;
      if( !RowMap.Comm().MyPID() ) NumLocalBs = BoundarySize;
     
      Epetra_Map LocalBoundaryIndexMap( BoundarySize, NumLocalBs, 0, RowMap.Comm() );
      if( verbose_ ) cout << "LocalBoundaryIndexMap:\n" << LocalBoundaryIndexMap;

      Epetra_IntVector lbGIDs( LocalBoundaryIndexMap );
      Epetra_Import lbImport( LocalBoundaryIndexMap, BoundaryIndexMap );
      lbGIDs.Import( bGIDs, lbImport, Insert );
      if( verbose_ ) cout << "LocalBoundaryGIDs:\n" << lbGIDs;

      Epetra_Map LocalBoundaryMap( BoundarySize, NumLocalBs, lbGIDs.Values(), 0, RowMap.Comm() );
      if( verbose_ ) cout << "LocalBoundaryMap:\n" << LocalBoundaryMap;

      Epetra_CrsGraph LocalBoundaryGraph( Copy, LocalBoundaryMap, LocalBoundaryMap, 0 );
      Epetra_Import LocalBoundaryImport( LocalBoundaryMap, Adj2.RowMap() );
      LocalBoundaryGraph.Import( Adj2, LocalBoundaryImport, Insert );
      LocalBoundaryGraph.TransformToLocal();
      if( verbose_ ) cout << "LocalBoundaryGraph:\n " << LocalBoundaryGraph;

      EpetraExt::CrsGraph_MapColoring BoundaryTrans( algo_, verbose_, reordering_ );
      Epetra_MapColoring & LocalBoundaryColoring = BoundaryTrans( LocalBoundaryGraph );
      if( verbose_ ) cout << "LocalBoundaryColoring:\n " << LocalBoundaryColoring;

      Epetra_Export BoundaryExport( LocalBoundaryMap, BoundaryMap );
      BoundaryColoring.Export( LocalBoundaryColoring, BoundaryExport, Insert );
    }
    else
    {
    /* Alternative Distrib. Memory Coloring of Boundary based on JonesPlassman(sic) paper
     * 1.Random number assignment to all boundary nodes using GID as seed to function
     * (This allows any processor to compute adj. off proc values with a local computation)
     * 2.Find all nodes greater than any neighbor off processor, color them.
     * 3.Send colored node info to neighbors
     * 4.Constrained color all nodes with all off proc neighbors smaller or colored.
     * 5.Goto 3
     */

      vector<int> OverlapBoundaryGIDs( boundaryGIDs );
      for( int i = nRows; i < Adj2.ColMap().NumMyElements(); ++i )
        OverlapBoundaryGIDs.push_back( Adj2.ColMap().GID(i) );

      int OverlapBoundarySize = OverlapBoundaryGIDs.size();
      Epetra_Map BoundaryColMap( -1, OverlapBoundarySize, &OverlapBoundaryGIDs[0], 0, RowMap.Comm() );

      Epetra_CrsGraph BoundaryGraph( Copy, BoundaryMap, BoundaryColMap, 0 );
      Epetra_Import BoundaryImport( BoundaryMap, Adj2.RowMap() );
      BoundaryGraph.Import( Adj2, BoundaryImport, Insert );
      BoundaryGraph.TransformToLocal();
      if( verbose_ ) cout << "BoundaryGraph:\n" << BoundaryGraph;

      Epetra_Import ReverseOverlapBoundaryImport( BoundaryMap, BoundaryColMap );
      Epetra_Import OverlapBoundaryImport( BoundaryColMap, BoundaryMap );

      int Colored = 0;
      int GlobalColored = 0;
      int Level = 0;
      Epetra_MapColoring OverlapBoundaryColoring( BoundaryColMap );

      if( verbose_ )
      {
        cout << MyPID << " PRands: ";
	for( int i = 0; i < LocalBoundarySize; ++i ) cout << PRAND(boundaryGIDs[i]) << " ";
	cout << endl;
      }

      while( GlobalColored < BoundarySize )
      {
        if( verbose_ ) cout << "Boundary Coloring Level: " << Level << endl;

	//Find current "Level" of boundary indices to color
	int NumIndices;
	int * Indices;
	vector<int> LevelIndices;
	for( int i = 0; i < LocalBoundarySize; ++i )
	{
          if( !OverlapBoundaryColoring[i] )
          {
            int MyVal = PRAND(BoundaryColMap.GID(i));
            BoundaryGraph.ExtractMyRowView( i, NumIndices, Indices );
	    bool ColorFlag = true;
	    int Loc = 0;
	    while( Loc<NumIndices && Indices[Loc]<LocalBoundarySize ) ++Loc;
	    for( int j = Loc; j < NumIndices; ++j )
              if( PRAND(BoundaryColMap.GID(Indices[j])) > MyVal
	          && !OverlapBoundaryColoring[Indices[j]] )
              {
                ColorFlag = false;
		break;
              }
            if( ColorFlag ) LevelIndices.push_back(i);
          }
        }

	if( verbose_ )
        {
          cout << MyPID << " Level Indices: ";
	  for( int i = 0; i < LevelIndices.size(); ++i ) cout << LevelIndices[i] << " ";
	  cout << endl;
        }

        //Greedy coloring of current level
        for( int i = 0; i < LevelIndices.size(); ++i )
        {
          BoundaryGraph.ExtractMyRowView( LevelIndices[i], NumIndices, Indices );

          set<int> usedColors;
          int color;
          for( int j = 0; j < NumIndices; ++j )
          {
            color = OverlapBoundaryColoring[ Indices[j] ];
            if( color > 0 ) usedColors.insert( color );
            color = 0;
            int testcolor = 1;
            while( !color )
            {
              if( !usedColors.count( testcolor ) ) color = testcolor;
              ++testcolor;
            }
          }
          OverlapBoundaryColoring[ LevelIndices[i] ] = color;
        }

	if( verbose_ ) cout << "Current Level Boundary Coloring:\n" << OverlapBoundaryColoring;

	//Update off processor coloring info
	BoundaryColoring.Import( OverlapBoundaryColoring, ReverseOverlapBoundaryImport, Insert );
	OverlapBoundaryColoring.Import( BoundaryColoring, OverlapBoundaryImport, Insert );
        Colored += LevelIndices.size();
	Level++;

	RowMap.Comm().SumAll( &Colored, &GlobalColored, 1 );
      }
    }

    if( verbose_ ) cout << "BoundaryColoring:\n " << BoundaryColoring;

    Epetra_MapColoring RowColorMap( RowMap );

    //Add Boundary Colors
    for( int i = 0; i < LocalBoundarySize; ++i )
    {
      int GID = BoundaryMap.GID(i);
      RowColorMap(GID) = BoundaryColoring(GID);
    }

    Epetra_MapColoring Adj2ColColorMap( Adj2.ColMap() );
    Epetra_Import Adj2Import( Adj2.ColMap(), RowMap );
    Adj2ColColorMap.Import( RowColorMap, Adj2Import, Insert );

    if( verbose_ ) cout << "RowColoringMap:\n " << RowColorMap;
    if( verbose_ ) cout << "Adj2ColColorMap:\n " << Adj2ColColorMap;

    vector<int> rowOrder( nRows );
    if( reordering_ == 0 || reordering_ == 1 ) 
    {
      multimap<int,int> adjMap;
      typedef multimap<int,int>::value_type adjMapValueType;
      for( int i = 0; i < nRows; ++i )
        adjMap.insert( adjMapValueType( Adj2.NumMyIndices(i), i ) );
      multimap<int,int>::iterator iter = adjMap.begin();
      multimap<int,int>::iterator end = adjMap.end();
      if( reordering_ == 0 ) //largest first (less colors)
      {
        for( int i = 1; iter != end; ++iter, ++i )
          rowOrder[nRows-i] = (*iter).second;
      }
      else                  //smallest first (better balance)
      {
        for( int i = 0; iter != end; ++iter, ++i )
          rowOrder[i] = (*iter).second;
      }
    }
    else if( reordering_ == 2 ) //random
    {
      for( int i = 0; i < nRows; ++i )
        rowOrder[i] = i;
#ifdef TFLOP
      random_shuffle( rowOrder.begin(), rowOrder.end() );
#endif
    }

    //Constrained greedy coloring of interior
    for( int row = 0; row < nRows; ++row )
    {
      if( !RowColorMap[ rowOrder[row] ] )
      {
        Adj2.ExtractMyRowView( rowOrder[row], NumIndices, Indices );

        set<int> usedColors;
        int color;
        for( int i = 0; i < NumIndices; ++i )
        {
          color = Adj2ColColorMap[ Indices[i] ];
          if( color > 0 ) usedColors.insert( color );
          color = 0;
          int testcolor = 1;
          while( !color )
          {
            if( !usedColors.count( testcolor ) ) color = testcolor;
            ++testcolor;
          }
        }
        Adj2ColColorMap[ rowOrder[row] ] = color;
      }
    }
    if( verbose_ ) cout << "RowColorMap after Greedy:\n " << RowColorMap;

    ColorMap = new Epetra_MapColoring( ColMap );
    Epetra_Import ColImport( ColMap, Adj2.ColMap() );
    ColorMap->Import( Adj2ColColorMap, ColImport, Insert );
  }

  if( verbose_ ) cout << "ColorMap!\n" << *ColorMap;

  newObj_ = ColorMap;

  return *ColorMap;
}

} // namespace EpetraExt
