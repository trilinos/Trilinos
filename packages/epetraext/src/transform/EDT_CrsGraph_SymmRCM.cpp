
#include <vector>
#include <set>

#include <Epetra_Util.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

#include <EDT_CrsGraph_Transpose.h>
#include <EDT_CrsGraph_SymmRCM.h>

namespace EpetraExt {

CrsGraph_SymmRCM::NewTypePtr CrsGraph_SymmRCM::operator()( CrsGraph_SymmRCM::OriginalTypeRef original )
{
  int err;

  //Generate Local Transpose Graph
  Epetra_CrsGraph * trans = CrsGraph_Transpose()( original );

  //Generate Local Symmetric Adj. List
  //Find Min Degree Node While at it
  int NumNodes = original.NumMyRows();
  int * LocalRow;
  int * LocalRowTrans;
  int RowSize, RowSizeTrans;
  vector< vector<int> > AdjList( NumNodes );
  int MinDegree = NumNodes;
  int MinDegreeNode;
  for( int i = 0; i < NumNodes; ++i )
  {
    original.ExtractMyRowView( i, RowSize, LocalRow );
    trans->ExtractMyRowView( i, RowSizeTrans, LocalRowTrans );
    int index = 0;
    int indexTrans = 0;
    while( ((index<RowSize)&&(LocalRow[index]<NumNodes)) ||
           ((indexTrans<RowSizeTrans)&&(LocalRowTrans[indexTrans]<NumNodes)) )
    {
      if( LocalRow[index] < LocalRowTrans[indexTrans] )
      {
        AdjList[i].push_back( LocalRow[index] );
        ++index;
      }
      else if( LocalRow[index] > LocalRowTrans[indexTrans] )
      {
        AdjList[i].push_back( LocalRowTrans[indexTrans] );
        ++indexTrans;
      }
      else
      {
        AdjList[i].push_back( LocalRow[index] );
        ++index;
        ++indexTrans;
      }
    }   
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
  BestBFT->NonNeighborLeaves(Leaves,5);

  BFT * DeeperBFT;
  bool DeeperFound;
  BFT * NarrowerBFT;
  bool NarrowerFound;
  
  bool Finished = false;

  while( !Finished )
  {
    DeeperBFT = 0;
    DeeperFound = false;

    NarrowerBFT = 0;
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
          if( DeeperBFT ) delete DeeperBFT;
          Diameter = TestBFT->Depth();
          BestWidth = TestBFT->Width();
          DeeperBFT = TestBFT;
          DeeperFound = true;
          if( NarrowerFound ) delete NarrowerBFT;
          NarrowerFound = false;
        }
        else if( (TestBFT->Depth()==Diameter) && (TestBFT->Width()<BestWidth) )
        {
          if( NarrowerBFT ) delete NarrowerBFT;
          NarrowerBFT = TestBFT;
          BestWidth = NarrowerBFT->Width();
          NarrowerFound = true;
        }
        else delete TestBFT;
      }
    }

    if( DeeperFound )
    {
      BestBFT = DeeperBFT;
      BestBFT->NonNeighborLeaves(Leaves,5);
    }
    else if( NarrowerFound )
    {
      BestBFT = NarrowerBFT;
      Finished = true;
    }
    else Finished = true;

  }

  vector<int> RCM;
  BestBFT->ReverseVector( RCM );

  //Generate New Row Map
  Epetra_Map * RCMMap = new Epetra_Map( original.RowMap().NumGlobalElements(),
                                        NumNodes,
                                        &RCM[0],
                                        original.RowMap().IndexBase(),
                                        original.RowMap().Comm() );


  //Create New Graph
  Epetra_Import Importer( *RCMMap, original.RowMap() );
  Epetra_CrsGraph * RCMGraph( new Epetra_CrsGraph( Copy, *RCMMap, 0 ) );
  RCMGraph->Import( original, Importer, Insert );
  RCMGraph->TransformToLocal();
  
  return RCMGraph;
}

CrsGraph_SymmRCM::BFT::BFT( const vector< vector<int> > & adjlist,
                            int root,
                            int max_width,
                            bool & failed )
: width_(0),
  depth_(0),
  nodes_(adjlist.size()),
  adjList_(adjlist)
{
  set<int> touchedNodes;
  levelSets_.push_back( vector<int>(1) );
  levelSets_[0][0] = root;
  touchedNodes.insert( root );

  while( touchedNodes.size() < (nodes_-1) )
  {
    levelSets_.push_back( vector<int>() );

    int currAdj;
    for( int i = 0; i < levelSets_[depth_].size(); ++i )
    {
      int currNode = levelSets_[depth_][i];

      for( int j = 0; j < adjlist[currNode].size(); ++j )
      {
        currAdj = adjlist[currNode][j];
        if( !touchedNodes.count( currAdj ) )
        {
          levelSets_[depth_+1].push_back( currAdj );
          touchedNodes.insert( currAdj );
        }
      }
    }

    int currWidth = levelSets_[depth_+1].size();

    if( currWidth )
    {
      if( currWidth > width_ ) width_ = currWidth;
      if( width_ > max_width ) { failed = true; return; }

      //Increasing Order By Degree
      vector<int> degrees( currWidth );
      for( int i = 0; i < currWidth; ++i )
        degrees[i] = adjList_[ levelSets_[depth_+1][i] ].size();
      int ** indices = new int*[1];
      indices[0] = &(levelSets_[depth_+1][0]);
      Epetra_Util().Sort( true, currWidth, &degrees[0],
                        0, 0, 1, indices );

      ++depth_;
    }
    else //it is a disconnected graph
    {
      bool found = false;
      for( int i = 0; i < nodes_; ++i )
        if( !touchedNodes.count( i ) )
        {
          levelSets_[depth_+1].push_back(i);
          ++depth_;
          found = true;
          i = nodes_;
        }
      if( !found ) { failed = true; return; }
    }

  }

}

void CrsGraph_SymmRCM::BFT::NonNeighborLeaves( vector<int> & leaves,
                                               int count )
{
  leaves.clear();
  int leafWidth = levelSets_[depth_].size();
  set<int> adjSeen;
  for( int i = 0; i < leafWidth; ++i )
  {
    int currLeaf = levelSets_[depth_][i];
    if( !adjSeen.count( currLeaf ) )
    {
      leaves.push_back( currLeaf );
      for( int j = 0; j < adjList_[currLeaf].size(); ++j )
        adjSeen.insert( adjList_[currLeaf][j] );
    }
    if( leaves.size() == count ) i = leafWidth;
  }
}

void CrsGraph_SymmRCM::BFT::ReverseVector( vector<int> & ordered )
{
  ordered.resize( nodes_ );
  int loc = 0;
  for( int i = 0; i < (depth_+1); ++i )
  {
    int currWidth = levelSets_[depth_-i].size();
    for( int j = 0; j < currWidth; ++j )
      ordered[loc++] = levelSets_[depth_-i][j];
  }
}

} //namespace EpetraExt

