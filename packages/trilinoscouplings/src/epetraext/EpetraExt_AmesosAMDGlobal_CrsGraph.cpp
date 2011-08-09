// @HEADER
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
// @HEADER
#include <EpetraExt_AmesosAMDGlobal_CrsGraph.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_Util.h>

#include <amesos_amd.h>

using std::vector;

namespace EpetraExt {

AmesosAMDGlobal_CrsGraph::
~AmesosAMDGlobal_CrsGraph()
{
}

AmesosAMDGlobal_CrsGraph::NewTypeRef
AmesosAMDGlobal_CrsGraph::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;
  const Epetra_BlockMap & OldRowMap = orig.RowMap();
  const Epetra_BlockMap & OldColMap = orig.ColMap();
  
  // Check if the graph is on one processor.
  int myMatProc = -1, matProc = -1;
  int myPID = orig.Comm().MyPID();
  int numProcs = orig.Comm().NumProc();
  for (int proc=0; proc<orig.Comm().NumProc(); proc++) 
  {
    if (orig.NumGlobalNonzeros() == orig.NumMyNonzeros())
      myMatProc = myPID;
  }
  orig.Comm().MaxAll( &myMatProc, &matProc, 1 );

  // Get some information about the parallel distribution.
  int maxMyRows = 0;
  int n = orig.NumMyRows();
  std::vector<int> numGlobalElem( numProcs );
  orig.Comm().GatherAll(&n, &numGlobalElem[0], 1);
  orig.Comm().MaxAll(&n, &maxMyRows, 1);

  Teuchos::RCP<Epetra_CrsGraph> serialGraph;
  Teuchos::RCP<Epetra_Map> serialMap;

  if( orig.RowMap().DistributedGlobal() && matProc == -1)
  {
    // The matrix is distributed and needs to be moved to processor zero.
    // Set the zero processor as the master.
    matProc = 0;
    serialMap = Teuchos::rcp( new Epetra_Map( Epetra_Util::Create_Root_Map( dynamic_cast<const Epetra_Map&>(orig.RowMap()), matProc ) ) );

    Epetra_Import serialImporter( *serialMap, OldRowMap );
    serialGraph = Teuchos::rcp( new Epetra_CrsGraph( Copy, *serialMap, 0 ) );
    serialGraph->Import( *origObj_, serialImporter, Insert );
    serialGraph->FillComplete();
  }
  if( orig.IndicesAreGlobal() && matProc == -1)
    { std::cout << "FAIL for Global Indices!\n"; abort(); }
 
  int nGlobal = orig.NumGlobalRows(); 
  int nnz = orig.NumMyNonzeros();
  
  if( debug_ )
  {
    std::cout << "Orig Graph:\n";
    orig.Print(std::cout);
  }

  // Extract the CRS graph pattern to vector.
  vector<int> ia(n+1,0);
  vector<int> ja(nnz);
  int ptr = 0, cnt = 0;

  for( int i = 0; i < n; ++i )
  {
    orig.ExtractMyRowCopy( i, nnz, cnt, &ja[ptr] );
    ia[i+1] = ia[i] + cnt;
    ptr += cnt;
  }
  
  // Compute the AMD permutation only on the processor that has the graph.
  if ( matProc == myPID ) {
    
    if( debug_ )
      {
	std::cout << "-----------------------------------------\n";
	std::cout << "CRS Format Graph \n";
	std::cout << "-----------------------------------------\n";
	for( int i = 0; i < n; ++i )
	  {
	    std::cout << ia[i] << " - " << ia[i+1] << " : ";
	    for( int j = ia[i]; j<ia[i+1]; ++j )
	      std::cout << " " << ja[j];
	    std::cout << std::endl;
	  }
	std::cout << "-----------------------------------------\n";
      }
    
    // Control and information vector
   double control[AMD_CONTROL];
   double info[AMD_INFO]; 
    
    // Storage for the permutation.
    perm_.resize( n );
   
    // Call AMD permutation.
 
    int ret = amesos_amd_order( n, &ia[0], &ja[0], &perm_[0], control, info );
    
    if( debug_ ) {
      std::cout << "-----------------------------------------\n";
      std::cout << "AMD Output (ret = " << ret << ", n = " << n << ")\n";
      std::cout << "-----------------------------------------\n";
      std::cout << "Perm\n";
      for( int i = 0; i<n; ++i )
	std::cout << perm_[i] << std::endl;
      std::cout << "-----------------------------------------\n";
    }  
  }

  // Broadcast the AMD permutation information to all processors.
  perm_.resize( nGlobal );

  orig.Comm().Broadcast(&perm_[0], nGlobal, matProc);
  
  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  vector<int> myElements( n );
  OldRowMap.MyGlobalElements( &myElements[0] );
  
  vector<int> newDomainElements( n );
  vector<int> newRangeElements( n );
  for( int i = 0; i < n; ++i )
  {
    newRangeElements[ i ] = myElements[ perm_[i] ];
    newDomainElements[ i ] = myElements[ perm_[i] ];
  }

  NewRowMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, n, &newRangeElements[0], OldRowMap.IndexBase(), OldRowMap.Comm() ) );
  NewColMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, n, &newDomainElements[0], OldColMap.IndexBase(), OldColMap.Comm() ) );

  if( debug_ )
  {
    std::cout << "New Row Map\n";
    std::cout << *NewRowMap_ << std::endl;
    std::cout << "New Col Map\n";
    std::cout << *NewColMap_ << std::endl;
  }

  //Generate New Graph
  NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 ) );
  Importer_ = Teuchos::rcp( new Epetra_Import( *NewRowMap_, OldRowMap ) );
  NewGraph_->Import( orig, *Importer_, Insert );
  NewGraph_->FillComplete();

  if( debug_ )
  {
    std::cout << "NewGraph\n";
    std::cout << *NewGraph_;
  }

  newObj_ = &*NewGraph_;

  return *NewGraph_;
}

bool
AmesosAMDGlobal_CrsGraph::
fwd()
{
  NewGraph_->Import( *origObj_, *Importer_, Insert );
  return true;
}

bool
AmesosAMDGlobal_CrsGraph::
rvs()
{
  origObj_->Export( *NewGraph_, *Importer_, Insert );
  return true;
}

} //namespace EpetraExt
