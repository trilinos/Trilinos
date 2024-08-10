// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <EpetraExt_AmesosAMDGlobal_CrsGraph.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <trilinos_btf_decl.h>
#include <trilinos_amd.h>

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

  if( orig.RowMap().DistributedGlobal() )
  { std::cout << "FAIL for Global!\n"; abort(); }
  if( orig.IndicesAreGlobal() )
  { std::cout << "FAIL for Global Indices!\n"; abort(); }

  // Extract the CCS information
  int n = orig.NumMyRows();
  int nnz = orig.NumMyNonzeros();

  vector<int> Ap(n+1,0);  // column pointers
  vector<int> Ai(nnz);    // row indices
  int cnt;
  for( int i = 0; i < n; ++i )
  {
    int * tmpP = &Ai[Ap[i]];
    orig.ExtractMyRowCopy( i, nnz-Ap[i], cnt, tmpP );
    Ap[i+1] = Ap[i] + cnt;
  }

  if (verbose_) 
  {
    std::cout << "-----------------------------------------\n";
    std::cout << "CRS Format Graph\n";
    std::cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
      {
	std::cout << Ap[i] << " - " << Ap[i+1] << " : ";
	for( int j = Ap[i]; j<Ap[i+1]; ++j )
	  std::cout << " " << Ai[j];
	std::cout << std::endl;
      }
    std::cout << "-----------------------------------------\n";
  }

  // Control and information vector
  double control[TRILINOS_AMD_CONTROL];
  double info[TRILINOS_AMD_INFO];

  // Storage for the permutation.
  perm_.resize( n );

  // Call AMD permutation.

  int ret = trilinos_amd_order( n, &Ap[0], &Ai[0], &perm_[0], control, info );

  if( debug_ ) {
    std::cout << "-----------------------------------------\n";
    std::cout << "AMD Output (ret = " << ret << ", n = " << n << ")\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "Perm\n";
    for( int i = 0; i<n; ++i )
      std::cout << perm_[i] << std::endl;
    std::cout << "-----------------------------------------\n";
  }

  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  const Epetra_BlockMap & OldMap = orig.RowMap();

  NewRowMap_ = Teuchos::rcp( new Epetra_Map( n, n, &perm_[0], OldMap.IndexBase(), OldMap.Comm() ) );
  NewColMap_ = NewRowMap_;

  if (verbose_) 
  {
    std::cout << "New Row Map\n";
    std::cout << *NewRowMap_ << std::endl;
    std::cout << "New Col Map\n";
    std::cout << *NewColMap_ << std::endl;
  }

  //Generate New Graph
  NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 ) );
  Epetra_Import Importer( *NewRowMap_, OldMap );
  NewGraph_->Import( orig, Importer, Insert );
  NewGraph_->FillComplete();

  if (verbose_) 
  {
    std::cout << "New CrsGraph\n";
    std::cout << *NewGraph_ << std::endl;
  }

  newObj_ = &*NewGraph_;

  return *NewGraph_;
}

} //namespace EpetraExt
