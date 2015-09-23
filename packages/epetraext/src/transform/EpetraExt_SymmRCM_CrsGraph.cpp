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

#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_EXPERIMENTAL
#ifdef HAVE_GRAPH_REORDERINGS

#include <EpetraExt_SymmRCM_CrsGraph.h>

#include <EpetraExt_Transpose_CrsGraph.h>

#include <set>

#include <Epetra_Util.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

namespace EpetraExt {

CrsGraph_SymmRCM::
~CrsGraph_SymmRCM()
{
  if( newObj_ ) delete newObj_;

  if( RCMColMap_ != RCMMap_ ) delete RCMColMap_;
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
  std::vector< std::vector<int> > AdjList( NumNodes );
  int MinDegree = NumNodes;
  int MinDegreeNode;
  for( int i = 0; i < NumNodes; ++i )
  {
    orig.ExtractMyRowView( i, RowSize, LocalRow );
    trans.ExtractMyRowView( i, RowSizeTrans, LocalRowTrans );

    std::set<int> adjSet;
    for( int j = 0; j < RowSize; ++j )
     if( LocalRow[j] < NumNodes ) adjSet.insert( LocalRow[j] );
    for( int j = 0; j < RowSizeTrans; ++j )
     if( LocalRowTrans[j] < NumNodes ) adjSet.insert( LocalRowTrans[j] );

    std::set<int>::iterator iterS = adjSet.begin();
    std::set<int>::iterator endS = adjSet.end();
    AdjList[i].resize( adjSet.size() );
    for( int j = 0; iterS != endS; ++iterS, ++j ) AdjList[i][j] = *iterS;
    
    if( AdjList[i].size() < MinDegree )
    {
      MinDegree = AdjList[i].size();
      MinDegreeNode = i;
    }
  }

  BFT * BestBFT;
  bool TooWide;

  //std::cout << "SymmRCM::bruteForce_ : " << bruteForce_ << std::endl;

  if( bruteForce_ )
  {
    int bestWidth = NumNodes;
    int bestDepth = 0;
    
    for( int i = 0; i < NumNodes; ++i )
    {
      BFT * TestBFT = new BFT( AdjList, i, NumNodes, TooWide );
      if( TestBFT->Depth() > bestDepth ||
          ( TestBFT->Depth() == bestDepth && TestBFT->Width() < bestWidth ) )
      {
        BestBFT = TestBFT;
        bestDepth = TestBFT->Depth();
        bestWidth = TestBFT->Width();
      }
      else
        delete TestBFT;
    }
  }
  else
  {
    //Construct BFT for first
    BestBFT = new BFT( AdjList, MinDegreeNode, NumNodes, TooWide );

    int MinWidth = BestBFT->Width();
    int BestWidth = MinWidth;
    int Diameter = BestBFT->Depth();
    std::vector<int> Leaves;
    BestBFT->NonNeighborLeaves( Leaves, AdjList, testLeafWidth_ );

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
        BestBFT->NonNeighborLeaves( Leaves, AdjList, testLeafWidth_ );
      else if( NarrowerFound )
        Finished = true;
      else Finished = true;
    }
  }

  //std::cout << "\nSymmRCM:\n";
  //std::cout << "----------------------------\n";
  //std::cout << " Depth: " << BestBFT->Depth() << std::endl;
  //std::cout << " Width: " << BestBFT->Width() << std::endl;
  //std::cout << "----------------------------\n\n";

  std::vector<int> RCM;
  BestBFT->ReverseVector( RCM );
  for( int i = 0; i < NumNodes; ++i )
    RCM[i] = orig.RowMap().GID( RCM[i] );

  //Generate New Row Map
  RCMMap_ = new Epetra_Map( orig.RowMap().NumGlobalElements(),
                                        NumNodes,
                                        &RCM[0],
                                        orig.RowMap().IndexBase(),
                                        orig.RowMap().Comm() );

  //Generate New Col Map
  if( RCMMap_->DistributedGlobal() )
  {
    std::vector<int> colIndices = RCM;
    const Epetra_BlockMap & origColMap = orig.ColMap();

    if( origColMap.NumMyElements() > RCMMap_->NumMyElements() )
    {
      for( int i = RCMMap_->NumMyElements(); i < origColMap.NumMyElements(); ++i )
        colIndices.push_back( origColMap.GID(i) );
    }

    RCMColMap_ = new Epetra_Map( orig.ColMap().NumGlobalElements(),
                                 colIndices.size(),
                                 &colIndices[0],
                                 orig.ColMap().IndexBase(),
                                 orig.ColMap().Comm() );
  }
  else
    RCMColMap_ = RCMMap_;

  //Create New Graph
  Epetra_Import Importer( *RCMMap_, orig.RowMap() );
  Epetra_CrsGraph * RCMGraph = new Epetra_CrsGraph( Copy, *RCMMap_, *RCMColMap_, 0 );
  RCMGraph->Import( orig, Importer, Insert );
  RCMGraph->FillComplete();

/*
  std::cout << "origGraph\n";
  std::cout << orig;
  std::cout << "RCMGraph\n";
  std::cout << *RCMGraph;
*/

  newObj_ = RCMGraph;
  
  return *RCMGraph;
}

CrsGraph_SymmRCM::BFT::
BFT( const std::vector< std::vector<int> > & adjlist,
     int root,
     int max_width,
     bool & failed )
: width_(0),
  depth_(0),
  nodes_(adjlist.size()),
  failed_(false)
{
  std::set<int> touchedNodes;

  //setup level 0 of traversal
  levelSets_.push_back( std::vector<int>(1) );
  levelSets_[0][0] = root;
  ++depth_;

  //start set of touched nodes
  touchedNodes.insert( root );

  while( touchedNodes.size() < nodes_ )
  {
    //start new level set
    levelSets_.push_back( std::vector<int>() );
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
      std::vector<int> degrees( currWidth );
      for( int i = 0; i < currWidth; ++i )
        degrees[i] = adjlist[ levelSets_[depth_-1][i] ].size();
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
        if( !touchedNodes.count( i ) && adjlist[i].size() < minDegree )
        {
          minDegree = adjlist[i].size();
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

/*
  std::cout << "BFT<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
  std::cout << "Width: " << width_ << std::endl;
  std::cout << "Depth: " << depth_ << std::endl;
  std::cout << "Adj List: " << nodes_ << std::endl;
  for( int i = 0; i < nodes_; ++i )
  {
    std::cout << i << "\t";
    for( int j = 0; j < adjlist[i].size(); ++j )
      std::cout << adjlist[i][j] << " ";
    std::cout << std::endl;
  }
  std::cout << "Level Sets: " << depth_ << std::endl;
  for( int i = 0; i < depth_; ++i )
  {
    std::cout << i << "\t";
    for( int j = 0; j < levelSets_[i].size(); ++j )
      std::cout << levelSets_[i][j] << " ";
    std::cout << std::endl;
  }
  std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
*/

  failed = failed_;
}

void
CrsGraph_SymmRCM::BFT::
NonNeighborLeaves( std::vector<int> & leaves,
                   const std::vector< std::vector<int> > & adjlist,
                   int count )
{
  assert( (depth_>0) && (failed_==false) );

  leaves.clear();
  int leafWidth = levelSets_[depth_-1].size();
  std::set<int> adjSeen;
  for( int i = 0; i < leafWidth; ++i )
  {
    int currLeaf = levelSets_[depth_-1][i];
    if( !adjSeen.count( currLeaf ) )
    {
      leaves.push_back( currLeaf );
      for( int j = 0; j < adjlist[currLeaf].size(); ++j )
        adjSeen.insert( adjlist[currLeaf][j] );
    }
    if( leaves.size() == count ) i = leafWidth;
  }
}

void
CrsGraph_SymmRCM::BFT::
ReverseVector( std::vector<int> & ordered )
{
  assert( (depth_>0) && (failed_==false) );

  ordered.resize( nodes_ );
  int loc = 0;
  for( int i = 0; i < depth_; ++i )
  {
    int currLevel = depth_ - (i+1);
    int currWidth = levelSets_[currLevel].size();
    for( int j = 0; j < currWidth; ++j )
      ordered[loc++] = levelSets_[currLevel][currWidth-j-1];
  }
}

} //namespace EpetraExt
#endif //HAVE_GRAPH_REORDERINGS
#endif //HAVE_EXPERIMENTAL
