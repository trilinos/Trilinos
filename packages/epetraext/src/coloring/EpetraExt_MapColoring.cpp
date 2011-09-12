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

#include <Epetra_Time.h>

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
  Epetra_Time timer( orig.Comm() );

  origObj_ = &orig;

  const Epetra_BlockMap & RowMap = orig.RowMap();
  int nRows = RowMap.NumMyElements();
  const Epetra_BlockMap & ColMap = orig.ColMap();
  int nCols = ColMap.NumMyElements();

  int MyPID = RowMap.Comm().MyPID();

  if( verbosity_ > 1 ) std::cout << "RowMap:\n" << RowMap;
  if( verbosity_ > 1 ) std::cout << "ColMap:\n" << ColMap;

  Epetra_MapColoring * ColorMap = 0;

  if( algo_ == GREEDY || algo_ == LUBY )
  {

    Epetra_CrsGraph * base = &( const_cast<Epetra_CrsGraph&>(orig) );

    int NumIndices;
    int * Indices;

    double wTime1 = timer.WallTime();

    // For parallel applications, add in boundaries to coloring
    bool distributedGraph = RowMap.DistributedGlobal();
    if( distributedGraph )
    {
      base = new Epetra_CrsGraph( Copy, ColMap, ColMap, 0 );

      for( int i = 0; i < nRows; ++i )
      {
        orig.ExtractMyRowView( i, NumIndices, Indices );
        base->InsertMyIndices( i, NumIndices, Indices );

        //Do this with a single insert
        //Is this the right thing?
        for( int j = 0; j < NumIndices; ++j )
          if( Indices[j] >= nRows )
            base->InsertMyIndices( Indices[j], 1, &i );
      } 
      base->FillComplete();
    }

    if( verbosity_ > 1 ) std::cout << "Base Graph:\n" << *base << std::endl;

    double wTime2 = timer.WallTime();
    if( verbosity_ > 0 )
      std::cout << "EpetraExt::MapColoring [INSERT BOUNDARIES] Time: " << wTime2-wTime1 << std::endl;

    //Generate Local Distance-1 Adjacency Graph
    //(Transpose of orig excluding non-local column indices)
    EpetraExt::CrsGraph_Transpose transposeTransform( true );
    Epetra_CrsGraph & Adj1 = transposeTransform( *base );
    if( verbosity_ > 1 ) std::cout << "Adjacency 1 Graph!\n" << Adj1;

    wTime1 = timer.WallTime();
    if( verbosity_ > 0 )
      std::cout << "EpetraExt::MapColoring [TRANSPOSE GRAPH] Time: " << wTime1-wTime2 << std::endl;

    int Delta = Adj1.MaxNumIndices();
    if( verbosity_ > 0 ) std::cout << std::endl << "Delta: " << Delta << std::endl;

    //Generation of Local Distance-2 Adjacency Graph
    Epetra_CrsGraph * Adj2 = &Adj1;
    if( !distance1_ )
    {
      Adj2 = new Epetra_CrsGraph( Copy, ColMap, ColMap, 0 );
      int NumAdj1Indices;
      int * Adj1Indices;
      for( int i = 0; i < nCols; ++i )
      {
        Adj1.ExtractMyRowView( i, NumAdj1Indices, Adj1Indices );

        set<int> Cols;
        for( int j = 0; j < NumAdj1Indices; ++j )
        {
          base->ExtractMyRowView( Adj1Indices[j], NumIndices, Indices );
          for( int k = 0; k < NumIndices; ++k )
            if( Indices[k] < nCols ) Cols.insert( Indices[k] );
        }
        int nCols2 = Cols.size();
        std::vector<int> ColVec( nCols2 );
        set<int>::iterator iterIS = Cols.begin();
        set<int>::iterator iendIS = Cols.end();
        for( int j = 0 ; iterIS != iendIS; ++iterIS, ++j ) ColVec[j] = *iterIS;
        Adj2->InsertMyIndices( i, nCols2, &ColVec[0] );
      }
      Adj2->FillComplete();

      if( verbosity_ > 1 ) std::cout << "Adjacency 2 Graph!\n" << *Adj2;

      wTime2 = timer.WallTime();
      if( verbosity_ > 0 )
        std::cout << "EpetraExt::MapColoring [GEN DIST-2 GRAPH] Time: " << wTime2-wTime1 << std::endl;
    }

    wTime2 = timer.WallTime();

    ColorMap = new Epetra_MapColoring( ColMap );

    std::vector<int> rowOrder( nCols );
    if( reordering_ == 0 || reordering_ == 1 ) 
    {
      std::multimap<int,int> adjMap;
      typedef std::multimap<int,int>::value_type adjMapValueType;
      for( int i = 0; i < nCols; ++i )
        adjMap.insert( adjMapValueType( Adj2->NumMyIndices(i), i ) );
      std::multimap<int,int>::iterator iter = adjMap.begin();
      std::multimap<int,int>::iterator end = adjMap.end();
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
      std::random_shuffle( rowOrder.begin(), rowOrder.end() );
#endif
    }

    wTime1 = timer.WallTime();
    if( verbosity_ > 0 )
      std::cout << "EpetraExt::MapColoring [REORDERING] Time: " << wTime1-wTime2 << std::endl;

    //Application of Greedy Algorithm to generate Color Map
    if( algo_ == GREEDY )
    {
      for( int col = 0; col < nCols; ++col )
      {
        Adj2->ExtractMyRowView( rowOrder[col], NumIndices, Indices );

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

      wTime2 = timer.WallTime();
      if( verbosity_ > 0 )
        std::cout << "EpetraExt::MapColoring [GREEDY COLORING] Time: " << wTime2-wTime1 << std::endl;
      if( verbosity_ > 0 )
        std::cout << "Num GREEDY Colors: " << ColorMap->NumColors() << std::endl;
    }
    else if( algo_ == LUBY )
    {
       //Assign Random Keys To Rows
      Epetra_Util util;
      std::vector<int> Keys(nCols);
      std::vector<int> State(nCols,-1);

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
              Adj2->ExtractMyRowView( col, NumIndices, Indices );
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
            if( State[col] < 0 )
              State[col] = CurrentColor;
          }

          //reinstate any zero not neighboring current color
          for( int col = 0; col < nCols; ++col )
            if( State[col] == 0 )
            {
              Adj2->ExtractMyRowView( col, NumIndices, Indices );
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

        if( verbosity_ > 2 )
        {
          std::cout << "Finished Color: " << CurrentColor << std::endl;
          std::cout << "NumRemaining: " << NumRemaining << std::endl;
        }

        //New color
        ++CurrentColor;
      }

      for( int col = 0; col < nCols; ++col )
        (*ColorMap)[col] = State[col]-1;

      wTime2 = timer.WallTime();
      if( verbosity_ > 0 )
        std::cout << "EpetraExt::MapColoring [LUBI COLORING] Time: " << wTime2-wTime1 << std::endl;
      if( verbosity_ > 0 )
        std::cout << "Num LUBI Colors: " << ColorMap->NumColors() << std::endl;

    }
    else
      abort(); //UNKNOWN ALGORITHM

    if( distributedGraph ) delete base;
    if( !distance1_ ) delete Adj2;
  }
  else
  {
    //Generate Parallel Adjacency-2 Graph
    EpetraExt::CrsGraph_Overlap OverlapTrans(1);
    Epetra_CrsGraph & OverlapGraph = OverlapTrans( orig );

    if( verbosity_ > 1 ) std::cout << "OverlapGraph:\n" << OverlapGraph;

    Epetra_CrsGraph * Adj2 = &orig;

    int NumIndices;
    int * Indices;
    if( !distance1_ )
    {
      Adj2 = new Epetra_CrsGraph( Copy, RowMap, OverlapGraph.ColMap(), 0 );
      int NumAdj1Indices;
      int * Adj1Indices;
      for( int i = 0; i < nRows; ++i )
      {
        OverlapGraph.ExtractMyRowView( i, NumAdj1Indices, Adj1Indices );

        set<int> Cols;
        for( int j = 0; j < NumAdj1Indices; ++j )
        {
          int GID = OverlapGraph.LRID( OverlapGraph.GCID( Adj1Indices[j] ) );
          OverlapGraph.ExtractMyRowView( GID, NumIndices, Indices );
          for( int k = 0; k < NumIndices; ++k ) Cols.insert( Indices[k] );
        }
        int nCols2 = Cols.size();
        std::vector<int> ColVec( nCols2 );
        set<int>::iterator iterIS = Cols.begin();
        set<int>::iterator iendIS = Cols.end();
        for( int j = 0 ; iterIS != iendIS; ++iterIS, ++j ) ColVec[j] = *iterIS;
        Adj2->InsertMyIndices( i, nCols2, &ColVec[0] );
      }
      int flag = Adj2->FillComplete();
      assert( flag == 0 );
      RowMap.Comm().Barrier();
      if( verbosity_ > 1 ) std::cout << "Adjacency 2 Graph!\n" << *Adj2;
    }

    //collect GIDs on boundary
    std::vector<int> boundaryGIDs;
    std::vector<int> interiorGIDs;
    for( int row = 0; row < nRows; ++row )
    {
      Adj2->ExtractMyRowView( row, NumIndices, Indices );
      bool testFlag = false;
      for( int i = 0; i < NumIndices; ++i )
        if( Indices[i] >= nRows ) testFlag = true;
      if( testFlag ) boundaryGIDs.push_back( Adj2->GRID(row) );
      else           interiorGIDs.push_back( Adj2->GRID(row) );
    }

    int LocalBoundarySize = boundaryGIDs.size();

    Epetra_Map BoundaryMap( -1, boundaryGIDs.size(),
      LocalBoundarySize ? &boundaryGIDs[0]: 0,
      0, RowMap.Comm() );
    if( verbosity_ > 1 ) std::cout << "BoundaryMap:\n" << BoundaryMap;
    
    int BoundarySize = BoundaryMap.NumGlobalElements();
    Epetra_MapColoring BoundaryColoring( BoundaryMap );

    if( algo_ == PSEUDO_PARALLEL )
    {
      Epetra_Map BoundaryIndexMap( BoundarySize, LocalBoundarySize, 0, RowMap.Comm() );
      if( verbosity_ > 1) std::cout << "BoundaryIndexMap:\n" << BoundaryIndexMap;

      Epetra_IntVector bGIDs( View, BoundaryIndexMap, &boundaryGIDs[0] );
      if( verbosity_ > 1) std::cout << "BoundaryGIDs:\n" << bGIDs;

      int NumLocalBs = 0;
      if( !RowMap.Comm().MyPID() ) NumLocalBs = BoundarySize;
     
      Epetra_Map LocalBoundaryIndexMap( BoundarySize, NumLocalBs, 0, RowMap.Comm() );
      if( verbosity_ > 1) std::cout << "LocalBoundaryIndexMap:\n" << LocalBoundaryIndexMap;

      Epetra_IntVector lbGIDs( LocalBoundaryIndexMap );
      Epetra_Import lbImport( LocalBoundaryIndexMap, BoundaryIndexMap );
      lbGIDs.Import( bGIDs, lbImport, Insert );
      if( verbosity_ > 1) std::cout << "LocalBoundaryGIDs:\n" << lbGIDs;

      Epetra_Map LocalBoundaryMap( BoundarySize, NumLocalBs, lbGIDs.Values(), 0, RowMap.Comm() );
      if( verbosity_ > 1) std::cout << "LocalBoundaryMap:\n" << LocalBoundaryMap;

      Epetra_CrsGraph LocalBoundaryGraph( Copy, LocalBoundaryMap, LocalBoundaryMap, 0 );
      Epetra_Import LocalBoundaryImport( LocalBoundaryMap, Adj2->RowMap() );
      LocalBoundaryGraph.Import( *Adj2, LocalBoundaryImport, Insert );
      LocalBoundaryGraph.FillComplete();
      if( verbosity_ > 1 ) std::cout << "LocalBoundaryGraph:\n " << LocalBoundaryGraph;

      EpetraExt::CrsGraph_MapColoring BoundaryTrans( GREEDY, reordering_, distance1_, verbosity_ );
      Epetra_MapColoring & LocalBoundaryColoring = BoundaryTrans( LocalBoundaryGraph );
      if( verbosity_ > 1 ) std::cout << "LocalBoundaryColoring:\n " << LocalBoundaryColoring;

      Epetra_Export BoundaryExport( LocalBoundaryMap, BoundaryMap );
      BoundaryColoring.Export( LocalBoundaryColoring, BoundaryExport, Insert );
    }
    else if( algo_ == JONES_PLASSMAN )
    {
    /* Alternative Distrib. Memory Coloring of Boundary based on JonesPlassman(sic) paper
     * 1.Random number assignment to all boundary nodes using GID as seed to function
     * (This allows any processor to compute adj. off proc values with a local computation)
     * 2.Find all nodes greater than any neighbor off processor, color them.
     * 3.Send colored node info to neighbors
     * 4.Constrained color all nodes with all off proc neighbors smaller or colored.
     * 5.Goto 3
     */

      std::vector<int> OverlapBoundaryGIDs( boundaryGIDs );
      for( int i = nRows; i < Adj2->ColMap().NumMyElements(); ++i )
        OverlapBoundaryGIDs.push_back( Adj2->ColMap().GID(i) );

      int OverlapBoundarySize = OverlapBoundaryGIDs.size();
      Epetra_Map BoundaryColMap( -1, OverlapBoundarySize,
        OverlapBoundarySize ? &OverlapBoundaryGIDs[0] : 0,
        0, RowMap.Comm() );

      Epetra_CrsGraph BoundaryGraph( Copy, BoundaryMap, BoundaryColMap, 0 );
      Epetra_Import BoundaryImport( BoundaryMap, Adj2->RowMap() );
      BoundaryGraph.Import( *Adj2, BoundaryImport, Insert );
      BoundaryGraph.FillComplete();
      if( verbosity_ > 1) std::cout << "BoundaryGraph:\n" << BoundaryGraph;

      Epetra_Import ReverseOverlapBoundaryImport( BoundaryMap, BoundaryColMap );
      Epetra_Import OverlapBoundaryImport( BoundaryColMap, BoundaryMap );

      int Colored = 0;
      int GlobalColored = 0;
      int Level = 0;
      Epetra_MapColoring OverlapBoundaryColoring( BoundaryColMap );

      //Setup random integers for boundary nodes
      Epetra_IntVector BoundaryValues( BoundaryMap );
      Epetra_Util Util;
      Util.SetSeed( 47954118 * (MyPID+1) );
      for( int i=0; i < LocalBoundarySize; ++i )
      {
        int val = Util.RandomInt();
        if( val < 0 ) val *= -1;
        BoundaryValues[i] = val;
      }

      //Get Random Values for External Boundary
      Epetra_IntVector OverlapBoundaryValues( BoundaryColMap );
      OverlapBoundaryValues.Import( BoundaryValues, OverlapBoundaryImport, Insert );

      while( GlobalColored < BoundarySize )
      {
	//Find current "Level" of boundary indices to color
	int NumIndices;
	int * Indices;
	std::vector<int> LevelIndices;
	for( int i = 0; i < LocalBoundarySize; ++i )
	{
          if( !OverlapBoundaryColoring[i] )
          {
            //int MyVal = PRAND(BoundaryColMap.GID(i));
            int MyVal = OverlapBoundaryValues[i];
            BoundaryGraph.ExtractMyRowView( i, NumIndices, Indices );
	    bool ColorFlag = true;
	    int Loc = 0;
	    while( Loc<NumIndices && Indices[Loc]<LocalBoundarySize ) ++Loc;
	    for( int j = Loc; j < NumIndices; ++j )
              if( (OverlapBoundaryValues[Indices[j]]>MyVal)
	          && !OverlapBoundaryColoring[Indices[j]] )
              {
                ColorFlag = false;
		break;
              }
            if( ColorFlag ) LevelIndices.push_back(i);
          }
        }

	if( verbosity_ > 1 )
        {
          std::cout << MyPID << " Level Indices: ";
	  int Lsize = (int) LevelIndices.size();
	  for( int i = 0; i < Lsize; ++i ) std::cout << LevelIndices[i] << " ";
	  std::cout << std::endl;
        }

        //Greedy coloring of current level
	set<int> levelColors;
	int Lsize = (int) LevelIndices.size();
        for( int i = 0; i < Lsize; ++i )
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
	  levelColors.insert( color );
        }

	if( verbosity_ > 2 ) std::cout << MyPID << " Level: " << Level << " Count: " << LevelIndices.size() << " NumColors: " << levelColors.size() << std::endl;

	if( verbosity_ > 2 ) std::cout << "Current Level Boundary Coloring:\n" << OverlapBoundaryColoring;

	//Update off processor coloring info
	BoundaryColoring.Import( OverlapBoundaryColoring, ReverseOverlapBoundaryImport, Insert );
	OverlapBoundaryColoring.Import( BoundaryColoring, OverlapBoundaryImport, Insert );
        Colored += LevelIndices.size();
	Level++;

	RowMap.Comm().SumAll( &Colored, &GlobalColored, 1 );
	if( verbosity_ > 2 ) std::cout << "Num Globally Colored: " << GlobalColored << " from Num Global Boundary Nodes: " << BoundarySize << std::endl;
      }
    }

    if( verbosity_ > 1 ) std::cout << "BoundaryColoring:\n " << BoundaryColoring;

    Epetra_MapColoring RowColorMap( RowMap );

    //Add Boundary Colors
    for( int i = 0; i < LocalBoundarySize; ++i )
    {
      int GID = BoundaryMap.GID(i);
      RowColorMap(GID) = BoundaryColoring(GID);
    }

    Epetra_MapColoring Adj2ColColorMap( Adj2->ColMap() );
    Epetra_Import Adj2Import( Adj2->ColMap(), RowMap );
    Adj2ColColorMap.Import( RowColorMap, Adj2Import, Insert );

    if( verbosity_ > 1 ) std::cout << "RowColoringMap:\n " << RowColorMap;
    if( verbosity_ > 1 ) std::cout << "Adj2ColColorMap:\n " << Adj2ColColorMap;

    std::vector<int> rowOrder( nRows );
    if( reordering_ == 0 || reordering_ == 1 ) 
    {
      std::multimap<int,int> adjMap;
      typedef std::multimap<int,int>::value_type adjMapValueType;
      for( int i = 0; i < nRows; ++i )
        adjMap.insert( adjMapValueType( Adj2->NumMyIndices(i), i ) );
      std::multimap<int,int>::iterator iter = adjMap.begin();
      std::multimap<int,int>::iterator end = adjMap.end();
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
    set<int> InteriorColors;
    for( int row = 0; row < nRows; ++row )
    {
      if( !RowColorMap[ rowOrder[row] ] )
      {
        Adj2->ExtractMyRowView( rowOrder[row], NumIndices, Indices );

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
	InteriorColors.insert( color );
      }
    }
    if( verbosity_ > 1 ) std::cout << MyPID << " Num Interior Colors: " << InteriorColors.size() << std::endl;
    if( verbosity_ > 1 ) std::cout << "RowColorMap after Greedy:\n " << RowColorMap;

    ColorMap = new Epetra_MapColoring( ColMap );
    Epetra_Import ColImport( ColMap, Adj2->ColMap() );
    ColorMap->Import( Adj2ColColorMap, ColImport, Insert );

    if( !distance1_ ) delete Adj2;
  }

  if( verbosity_ > 0 ) std::cout << MyPID << " ColorMap Color Count: " << ColorMap->NumColors() << std::endl;
  if( verbosity_ > 1 ) std::cout << "ColorMap!\n" << *ColorMap;

  newObj_ = ColorMap;

  return *ColorMap;
}

} // namespace EpetraExt
