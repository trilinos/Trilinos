// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <EpetraExt_AmesosBTFGlobal_LinearProblem.h>

#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Util.h>
	
#include <EpetraExt_AmesosBTF_CrsMatrix.h>
#include <EpetraExt_Reindex_CrsMatrix.h>
#include <EpetraExt_BlockAdjacencyGraph.h>

#include <Teuchos_ParameterList.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

using std::vector;

namespace EpetraExt {

AmesosBTFGlobal_LinearProblem::
~AmesosBTFGlobal_LinearProblem()
{
  if ( newObj_ ) delete newObj_;
}

AmesosBTFGlobal_LinearProblem::NewTypeRef
AmesosBTFGlobal_LinearProblem::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  // Extract the matrix and vectors from the linear problem
  OldRHS_ = Teuchos::rcp( orig.GetRHS(), false );
  OldLHS_ = Teuchos::rcp( orig.GetLHS(), false );
  OldMatrix_ = Teuchos::rcp( dynamic_cast<Epetra_CrsMatrix *>( orig.GetMatrix() ), false );
	
  int nGlobal = OldMatrix_->NumGlobalRows(); 
  int n = OldMatrix_->NumMyRows();

  // Check if the matrix is on one processor.
  int myMatProc = -1, matProc = -1;
  int myPID = OldMatrix_->Comm().MyPID();
  int numProcs = OldMatrix_->Comm().NumProc();

  const Epetra_BlockMap& oldRowMap = OldMatrix_->RowMap();

  // Get some information about the parallel distribution.
  int maxMyRows = 0;
  std::vector<int> numGlobalElem( numProcs );
  OldMatrix_->Comm().GatherAll(&n, &numGlobalElem[0], 1);
  OldMatrix_->Comm().MaxAll(&n, &maxMyRows, 1);

  for (int proc=0; proc<numProcs; proc++) 
  {
    if (OldMatrix_->NumGlobalNonzeros() == OldMatrix_->NumMyNonzeros())
      myMatProc = myPID;
  }

  OldMatrix_->Comm().MaxAll( &myMatProc, &matProc, 1 );

  Teuchos::RCP<Epetra_CrsMatrix> serialMatrix;
  Teuchos::RCP<Epetra_Map> serialMap;	
  if( oldRowMap.DistributedGlobal() && matProc == -1) 
  {
    // The matrix is distributed and needs to be moved to processor zero.
    // Set the zero processor as the master.
    matProc = 0;
    serialMap = Teuchos::rcp( new Epetra_Map( Epetra_Util::Create_Root_Map( OldMatrix_->RowMap(), matProc ) ) );
    
    Epetra_Import serialImporter( *serialMap, OldMatrix_->RowMap() );
    serialMatrix = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *serialMap, 0 ) );
    serialMatrix->Import( *OldMatrix_, serialImporter, Insert );
    serialMatrix->FillComplete();
  }
  else {
    // The old matrix has already been moved to one processor (matProc).
    serialMatrix = OldMatrix_;
  }

  if( debug_ )
  {
    std::cout << "Original (serial) Matrix:\n";
    std::cout << *serialMatrix << std::endl;
  }

  // Obtain the current row and column orderings
  std::vector<int> origGlobalRows(nGlobal), origGlobalCols(nGlobal);
  serialMatrix->RowMap().MyGlobalElements( &origGlobalRows[0] );
  serialMatrix->ColMap().MyGlobalElements( &origGlobalCols[0] );
  
  // Perform reindexing on the full serial matrix (needed for BTF).
  Epetra_Map reIdxMap( serialMatrix->RowMap().NumGlobalElements(), serialMatrix->RowMap().NumMyElements(), 0, serialMatrix->Comm() );
  Teuchos::RCP<EpetraExt::ViewTransform<Epetra_CrsMatrix> > reIdxTrans =
    Teuchos::rcp( new EpetraExt::CrsMatrix_Reindex( reIdxMap ) );
  Epetra_CrsMatrix newSerialMatrix = (*reIdxTrans)( *serialMatrix );
  reIdxTrans->fwd();
  
  // Compute and apply BTF to the serial CrsMatrix and has been filtered by the threshold
  EpetraExt::AmesosBTF_CrsMatrix BTFTrans( threshold_, upperTri_, verbose_, debug_ );
  Epetra_CrsMatrix newSerialMatrixBTF = BTFTrans( newSerialMatrix );
  
  rowPerm_ = BTFTrans.RowPerm();
  colPerm_ = BTFTrans.ColPerm();
  blockPtr_ = BTFTrans.BlockPtr();
  numBlocks_ = BTFTrans.NumBlocks();
 
  if ( (myPID == matProc) && verbose_) {
    bool isSym = true;
    for (int i=0; i<nGlobal; ++i) {
      if (rowPerm_[i] != colPerm_[i]) {
        isSym = false;
        break;
      }
    }
    std::cout << "The BTF permutation symmetry (0=false,1=true) is : " << isSym << std::endl;
  }
  
  // Compute the permutation w.r.t. the original row and column GIDs.
  std::vector<int> origGlobalRowsPerm(nGlobal), origGlobalColsPerm(nGlobal);
  if (myPID == matProc) {
    for (int i=0; i<nGlobal; ++i) {
      origGlobalRowsPerm[i] = origGlobalRows[ rowPerm_[i] ];
      origGlobalColsPerm[i] = origGlobalCols[ colPerm_[i] ];
    }
  }
  OldMatrix_->Comm().Broadcast( &origGlobalRowsPerm[0], nGlobal, matProc );
  OldMatrix_->Comm().Broadcast( &origGlobalColsPerm[0], nGlobal, matProc );

  // Generate the full serial matrix that imports according to the previously computed BTF.
  Epetra_CrsMatrix newSerialMatrixT( Copy, newSerialMatrixBTF.RowMap(), 0 );
  newSerialMatrixT.Import( newSerialMatrix, *(BTFTrans.Importer()), Insert );
  newSerialMatrixT.FillComplete();
  
  if( debug_ )
  {
    std::cout << "Original (serial) Matrix permuted via BTF:\n";
    std::cout << newSerialMatrixT << std::endl;
  }

  // Perform reindexing on the full serial matrix (needed for balancing).
  Epetra_Map reIdxMap2( newSerialMatrixT.RowMap().NumGlobalElements(), newSerialMatrixT.RowMap().NumMyElements(), 0, newSerialMatrixT.Comm() );
  Teuchos::RCP<EpetraExt::ViewTransform<Epetra_CrsMatrix> > reIdxTrans2 =
    Teuchos::rcp( new EpetraExt::CrsMatrix_Reindex( reIdxMap2 ) );
  Epetra_CrsMatrix tNewSerialMatrixT = (*reIdxTrans2)( newSerialMatrixT );
  reIdxTrans2->fwd();

  Teuchos::RCP<Epetra_Map> balancedMap;
  
  if (balance_ == "linear") 
  {
    // Distribute block somewhat evenly across processors
    std::vector<int> rowDist(numProcs+1,0);
    int balRows = nGlobal / numProcs + 1;
    int numRows = balRows, currProc = 1;
    for ( int i=0; i<numBlocks_ || currProc < numProcs; ++i ) {
      if (blockPtr_[i] > numRows) {
	rowDist[currProc++] = blockPtr_[i-1];
	numRows = blockPtr_[i-1] + balRows;
      }      
    }
    rowDist[numProcs] = nGlobal;
   
    // Create new Map based on this linear distribution.
    int numMyBalancedRows = rowDist[myPID+1]-rowDist[myPID];

    NewRowMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, numMyBalancedRows, &origGlobalRowsPerm[ rowDist[myPID] ], 0, OldMatrix_->Comm() ) );
    // Right now we do not explicitly build the column map and assume the BTF permutation is symmetric!
    //NewColMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, nGlobal, &colPerm_[0], 0, OldMatrix_->Comm() ) );
    
    if ( verbose_ ) 
    {
      std::cout << "Processor " << myPID << " has " << numMyBalancedRows << " rows." << std::endl;    
    }
  }
  else if (balance_ == "isorropia") 
  {
    // Compute block adjacency graph for partitioning.
    std::vector<double> weight;
    EpetraExt::BlockAdjacencyGraph adjGraph;
    Teuchos::RCP<const Epetra_CrsGraph> blkGraph =
      adjGraph.compute( const_cast<Epetra_CrsGraph&>(tNewSerialMatrixT.Graph()), 
                                                     numBlocks_, blockPtr_, weight, verbose_);

    // Get row weights for rebalancing the graph.
    Teuchos::RCP<const Epetra_Vector> rowWeights = Teuchos::rcp( new Epetra_Vector( View, blkGraph->Map(), &weight[0] ) );

    // Create the cost describer object for Isorropia, pass in row weights.
    Teuchos::RCP<Isorropia::Epetra::CostDescriber> costPtr = 
      Teuchos::rcp( new Isorropia::Epetra::CostDescriber() );
    costPtr->setVertexWeights(rowWeights);
   
    // Create the partitioner object, pass in original graph, cost describer, and empty parameter list.
    Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
      Teuchos::rcp(new Isorropia::Epetra::Partitioner(blkGraph, costPtr));

    // Create the redistributor object.
    Isorropia::Epetra::Redistributor rd(partitioner);

    // Call Isorropia to rebalance this graph.
    Teuchos::RCP<Epetra_CrsGraph> balancedGraph = rd.redistribute( *blkGraph );
 
    int myNumBlkRows = balancedGraph->NumMyRows();    
    
    //std::vector<int> myGlobalElements(nGlobal);
    std::vector<int> newRangeElements(nGlobal), newDomainElements(nGlobal);
    int grid = 0, myElements = 0;
    for (int i=0; i<myNumBlkRows; ++i) {
      grid = balancedGraph->GRID( i );
      for (int j=blockPtr_[grid]; j<blockPtr_[grid+1]; ++j) {
	newRangeElements[myElements++] = origGlobalRowsPerm[j];
	//myGlobalElements[myElements++] = j;
      }
    }

    NewRowMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, myElements, &newRangeElements[0], 0, OldMatrix_->Comm() ) );
    // Right now we do not explicitly build the column map and assume the BTF permutation is symmetric!
    //NewColMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, nGlobal, &colPerm_[0], 0, OldMatrix_->Comm() ) );
    //balancedMap = Teuchos::rcp( new Epetra_Map( nGlobal, myElements, &myGlobalElements[0], 0, serialMatrix->Comm() ) );

    if ( verbose_ ) 
    {
      std::cout << "Processor " << myPID << " has " << myElements << " rows." << std::endl;
    }
  }
  
  // Use New Domain and Range Maps to Generate Importer
  //for now, assume they start out as identical
  Epetra_Map OldRowMap = OldMatrix_->RowMap();
  Epetra_Map OldColMap = OldMatrix_->ColMap();
  
  if( debug_ )
  {
    std::cout << "New Row Map\n";
    std::cout << *NewRowMap_ << std::endl;
    //std::cout << "New Col Map\n";
    //std::cout << *NewColMap_ << std::endl;
  }

  // Generate New Graph
  // NOTE:  Right now we are creating the graph, assuming that the permutation is symmetric!
  // NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 ) );
  NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *NewRowMap_, 0 ) );
  Importer_ = Teuchos::rcp( new Epetra_Import( *NewRowMap_, OldRowMap ) );
  Importer2_ = Teuchos::rcp( new Epetra_Import( OldRowMap, *NewRowMap_ ) );
  NewGraph_->Import( OldMatrix_->Graph(), *Importer_, Insert );
  NewGraph_->FillComplete();

  if( debug_ )
  {
    std::cout << "NewGraph\n";
    std::cout << *NewGraph_;
  }

  // Create new linear problem and import information from old linear problem
  NewMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *NewGraph_ ) );
  NewMatrix_->Import( *OldMatrix_, *Importer_, Insert );
  NewMatrix_->FillComplete();

  NewLHS_ = Teuchos::rcp( new Epetra_MultiVector( *NewRowMap_, OldLHS_->NumVectors() ) );
  NewLHS_->Import( *OldLHS_, *Importer_, Insert );
  
  NewRHS_ = Teuchos::rcp( new Epetra_MultiVector( *NewRowMap_, OldRHS_->NumVectors() ) );
  NewRHS_->Import( *OldRHS_, *Importer_, Insert );

  if( debug_ )
  {
    std::cout << "New Matrix\n";
    std::cout << *NewMatrix_ << std::endl;
  }

  newObj_ = new Epetra_LinearProblem( &*NewMatrix_, &*NewLHS_, &*NewRHS_ );

  return *newObj_;
}

bool
AmesosBTFGlobal_LinearProblem::
fwd()
{
  NewLHS_->Import( *OldLHS_, *Importer_, Insert );
  NewRHS_->Import( *OldRHS_, *Importer_, Insert );
  NewMatrix_->Import( *OldMatrix_, *Importer_, Insert );

  return true;
}

bool
AmesosBTFGlobal_LinearProblem::
rvs()
{
  //  std::cout << "AmesosBTFGlobal_LinearProblem: NewLHS_" << std::endl;
  //  std::cout << *NewLHS_ << std::endl;
  OldLHS_->Import( *NewLHS_, *Importer2_, Insert );
  //std::cout << "AmesosBTFGlobal_LinearProblem: OldLHS_" << std::endl;
  //std::cout << *OldLHS_ << std::endl;

  return true;
}

} //namespace EpetraExt
