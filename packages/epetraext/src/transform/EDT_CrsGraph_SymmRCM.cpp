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

#include <set>

#include <Epetra_Util.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

#include <EDT_CrsGraph_Transpose.h>
#include <EDT_CrsGraph_SymmRCM.h>

using std::set;
using std::vector;

namespace EpetraExt {

CrsGraph_SymmRCM::
~CrsGraph_SymmRCM()
{
  if( newObj_ ) delete newObj_;

  if( RCMMap_ ) delete RCMMap_;
}

CrsGraph_SymmRCM::NewTypeRef
CrsGraph_SymmRCM::
operator()( CrsGraph_SymmRCM::OriginalTypeRef orig )
{
  origObj_ = &orig;

  int err;

  //Generate Local Transpose Graph
  CrsGraph_Transpose transposeTransform;
  Epetra_CrsGraph & trans = transposeTransform( orig );

  //Generate Local Symmetric Adj. List
  //Find Min Degree Node While at it
  int NumNodes = orig.NumMyRows();
  int * LocalRow;
  int * LocalRowTrans;
  int RowSize, RowSizeTrans;
  vector< vector<int> > AdjList( NumNodes );
  int MinDegree = NumNodes;
  int MinDegreeNode;
  for( int i = 0; i < NumNodes; ++i )
  {
    orig.ExtractMyRowView( i, RowSize, LocalRow );
    trans.ExtractMyRowView( i, RowSizeTrans, LocalRowTrans );

    set<int> adjSet;
    for( int j = 0; j < RowSize; ++j ) adjSet.insert( LocalRow[j] );
    for( int j = 0; j < RowSizeTrans; ++j ) adjSet.insert( LocalRowTrans[j] );

    set<int>::iterator iterS = adjSet.begin();
    set<int>::iterator endS = adjSet.end();
    AdjList[i].resize( adjSet.size() );
    for( int j = 0; iterS != endS; ++iterS, ++j ) AdjList[i][j] = *iterS;
    
    if( AdjList[i].size() < MinDegree )
    {
      MinDegree = AdjList[i].size();
      MinDegreeNode = i;
    }
  }

  //Construct BFT for first
  bool TooWide;
  BFT * BestBFT = new BFT( AdjList, MinDegreeNode, NumNodes, TooWide );

  int MinWidth = BestBFT->Width();
  int BestWidth = MinWidth;
  int Diameter = BestBFT->Depth();
  vector<int> Leaves;
  BestBFT->NonNeighborLeaves( Leaves, testLeafWidth_ );

  bool DeeperFound;
  bool NarrowerFound;
  
  bool Finished = false;

  while( !Finished )
  {
    DeeperFound = false;
    NarrowerFound = false;

    for( int i = 0; i < Leaves.size(); ++i )
    {

      BFT * TestBFT = new BFT( AdjList, Leaves[i], MinWidth, TooWide );

      if( TooWide )
        delete TestBFT;
      else
      {
        if( TestBFT->Width() < MinWidth ) MinWidth = TestBFT->Width();

        if( TestBFT->Depth() > Diameter )
        {
          delete BestBFT;
          Diameter = TestBFT->Depth();
          BestWidth = TestBFT->Width();
          BestBFT = TestBFT;
          DeeperFound = true;
          NarrowerFound = false;
        }
        else if( (TestBFT->Depth()==Diameter) && (TestBFT->Width()<BestWidth) )
        {
          delete BestBFT;
          BestWidth = TestBFT->Width();
          BestBFT = TestBFT;
          NarrowerFound = true;
        }
        else delete TestBFT;
      }
    }

    if( DeeperFound )
      BestBFT->NonNeighborLeaves(Leaves,5);
    else if( NarrowerFound )
      Finished = true;
    else Finished = true;

  }

  vector<int> RCM;
  BestBFT->ReverseVector( RCM );

  //Generate New Row Map
  RCMMap_ = new Epetra_Map( orig.RowMap().NumGlobalElements(),
                                        NumNodes,
                                        &RCM[0],
                                        orig.RowMap().IndexBase(),
                                        orig.RowMap().Comm() );


  //Create New Graph
  Epetra_Import Importer( *RCMMap_, orig.RowMap() );
  Epetra_CrsGraph * RCMGraph( new Epetra_CrsGraph( Copy, *RCMMap_, 0 ) );
  RCMGraph->Import( orig, Importer, Insert );
  RCMGraph->TransformToLocal();

  newObj_ = RCMGraph;
  
  return *RCMGraph;
}

CrsGraph_SymmRCM::BFT::
BFT( const vector< vector<int> > & adjlist,
     int root,
     int max_width,
     bool & failed )
: width_(0),
  depth_(0),
  nodes_(adjlist.size()),
  adjList_(adjlist),
  failed_(false)
{
  set<int> touchedNodes;

  //setup level 0 of traversal
  levelSets_.push_back( vector<int>(1) );
  levelSets_[0][0] = root;
  ++depth_;

  //start set of touched nodes
  touchedNodes.insert( root );

  while( touchedNodes.size() < nodes_ )
  {
    //start new level set
    levelSets_.push_back( vector<int>() );
    ++depth_;

    for( int i = 0; i < levelSets_[depth_-2].size(); ++i )
    {
      int currNode = levelSets_[depth_-2][i];
      int adjSize  = adjlist[currNode].size();
      for( int j = 0; j < adjSize; ++j )
      {
        // add nodes to current level set when new
        int currAdj = adjlist[currNode][j];
        if( !touchedNodes.count( currAdj ) )
        {
          levelSets_[depth_-1].push_back( currAdj );
          touchedNodes.insert( currAdj );
        }
      }
    }

    int currWidth = levelSets_[depth_-1].size();

    if( currWidth ) //sort adj nodes by degree
    {
      if( currWidth > width_ ) width_ = currWidth;

      //fail if width is greater than allowed
      if( width_ > max_width )
      {
        failed_ = true;
        failed = failed_;
        return;
      }

      //Increasing Order By Degree
      vector<int> degrees( currWidth );
      for( int i = 0; i < currWidth; ++i )
        degrees[i] = adjList_[ levelSets_[depth_-1][i] ].size();
      int ** indices = new int*[1];
      indices[0] = &(levelSets_[depth_-1][0]);
      Epetra_Util().Sort( true, currWidth, &degrees[0], 0, 0, 1, indices );
    }
    else //it is a disconnected graph
    {
      //start again from minimum degree node of those remaining
      bool found = false;
      int minDegree = nodes_;
      int minDegreeNode;
      for( int i = 0; i < nodes_; ++i )
      {
        if( !touchedNodes.count( i ) && adjList_[i].size() < minDegree )
        {
          minDegree = adjList_[i].size();
          minDegreeNode = i;
          found = true;
        }
      }

      if( found )
      {
        touchedNodes.insert( minDegreeNode );
        levelSets_[depth_-1].push_back( minDegreeNode );
      }
      else
      {
        --depth_;
        failed_ = true;
        failed = failed_;
        return;
      }

    }

  }

  cout << "BFT<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
  cout << "Width: " << width_ << endl;
  cout << "Depth: " << depth_ << endl;
  cout << "Adj List: " << nodes_ << endl;
  for( int i = 0; i < nodes_; ++i )
  {
    cout << i << "\t";
    for( int j = 0; j < adjList_[i].size(); ++j )
      cout << adjList_[i][j] << " ";
    cout << endl;
  }
  cout << "Level Sets: " << depth_ << endl;
  for( int i = 0; i < depth_; ++i )
  {
    cout << i << "\t";
    for( int j = 0; j < levelSets_[i].size(); ++j )
      cout << levelSets_[i][j] << " ";
    cout << endl;
  }
  cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";

  failed = failed_;
}

void
CrsGraph_SymmRCM::BFT::
NonNeighborLeaves( vector<int> & leaves,
                   int count )
{
  assert( (depth_>0) && (failed_==false) );

  leaves.clear();
  int leafWidth = levelSets_[depth_-1].size();
  set<int> adjSeen;
  for( int i = 0; i < leafWidth; ++i )
  {
    int currLeaf = levelSets_[depth_-1][i];
    if( !adjSeen.count( currLeaf ) )
    {
      leaves.push_back( currLeaf );
      for( int j = 0; j < adjList_[currLeaf].size(); ++j )
        adjSeen.insert( adjList_[currLeaf][j] );
    }
    if( leaves.size() == count ) i = leafWidth;
  }
}

void
CrsGraph_SymmRCM::BFT::
ReverseVector( vector<int> & ordered )
{
  assert( (depth_>0) && (failed_==false) );

  ordered.resize( nodes_ );
  int loc = 0;
  for( int i = 0; i < (depth_-1); ++i )
  {
    int currLevel = depth_ - (i+1);
    int currWidth = levelSets_[currLevel].size();
    for( int j = 0; j < currWidth; ++j )
      ordered[loc++] = levelSets_[currLevel][j];
  }
}

} //namespace EpetraExt

