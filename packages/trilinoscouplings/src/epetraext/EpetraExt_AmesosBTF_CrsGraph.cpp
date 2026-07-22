// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <EpetraExt_AmesosBTF_CrsGraph.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <trilinos_btf_decl.h>

using std::vector;

namespace EpetraExt {

AmesosBTF_CrsGraph::
~AmesosBTF_CrsGraph()
{
}

AmesosBTF_CrsGraph::NewTypeRef
AmesosBTF_CrsGraph::
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
  vector<int> Ai_tmp(nnz);

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

  // Transpose the graph, not the values
  int j=0, next=0;
  vector<int> Ap_tmp(n+1,0);

  // Compute row lengths
  for (int i = 0; i < n; i++)
      for (int k = Ap[i]; k < Ap[i+1]; k++)
          ++Ap_tmp[ Ai[k]+1 ];

  // Compute pointers from row lengths
  Ap_tmp[0] = 0;
  for (int i = 0; i < n; i++)
      Ap_tmp[i+1] += Ap_tmp[i];

  // Copy over indices
  for (int i = 0; i < n; i++) {
      for (int k = Ap[i]; k < Ap[i+1]; k++) {
          j = Ai[k];
          next = Ap_tmp[j];
          Ai_tmp[next] = i;
          Ap_tmp[j] = next + 1;
      }
  }

  // Reshift Ap_tmp
  for (int i=n-1; i >= 0; i--) Ap_tmp[i+1] = Ap_tmp[i];
  Ap_tmp[0] = 0;

  // Transformation information
  int numMatch = 0;       // number of nonzeros on diagonal after permutation.
  double maxWork =  0.0;  // no limit on how much work to perform in max-trans.
  double workPerf = 0.0;  // how much work was performed in max-trans.

  // Create a work vector for the BTF code.
  vector<int> work(5*n);

  // Storage for the row and column permutations.
  vector<int> rowperm(n);
  vector<int> colperm(n);
  vector<int> blkPtr(n+1);

  // NOTE:  The permutations are sent in backwards since the matrix is transposed.
  // On output, rowperm and colperm are the row and column permutations of A, where 
  // i = TRILINOS_BTF_UNFLIP(rowperm[k]) if row i of A is the kth row of P*A*Q, and j = colperm[k] 
  // if column j of A is the kth column of P*A*Q.  If rowperm[k] < 0, then the 
  // (k,k)th entry in P*A*Q is structurally zero.

  numBlocks_ = trilinos_btf_order( n, &Ap_tmp[0], &Ai_tmp[0], maxWork, &workPerf,
			  &rowperm[0], &colperm[0], &blkPtr[0], 
			  &numMatch, &work[0] );

  // Reverse ordering of permutation to get upper triangular form, if necessary
  rowPerm_.resize( n );
  colPerm_.resize( n );
  blkPtr_.resize( numBlocks_+1 );
  if (!upperTri_) {
    for( int i = 0; i < n; ++i )
    {
      rowPerm_[i] = TRILINOS_BTF_UNFLIP(rowperm[(n-1)-i]);
      colPerm_[i] = colperm[(n-1)-i];
    }
    for( int i = 0; i < numBlocks_+1; ++i )
      blkPtr_[i] = n - blkPtr[numBlocks_-i];
  }
  else {
    colPerm_ = colperm;
    blkPtr_ = blkPtr;
    for( int i = 0; i < n; ++i )
      rowPerm_[i] = TRILINOS_BTF_UNFLIP(rowperm[i]);
  }

  if (verbose_) 
  {
    std::cout << "-----------------------------------------\n";
    std::cout << "BTF Output (n = " << n << ")\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "Num Blocks: " << numBlocks_ << std::endl;
    std::cout << "Num NNZ Diags: " << numMatch << std::endl;
    std::cout << "RowPerm and ColPerm \n";
    for( int i = 0; i<n; ++i )
      std::cout << rowPerm_[i] << "\t" << colPerm_[i] << std::endl;
    std::cout << "-----------------------------------------\n";
  }

  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  const Epetra_BlockMap & OldMap = orig.RowMap();
  vector<int> newDomainElements( n );
  vector<int> newRangeElements( n );
  for( int i = 0; i < n; ++i )
  {
    newRangeElements[ i ] = rowPerm_[i];
    newDomainElements[ i ] = colPerm_[i];
  }

  NewRowMap_ = Teuchos::rcp( new Epetra_Map( n, n, &newRangeElements[0], OldMap.IndexBase(), OldMap.Comm() ) );
  NewColMap_ = Teuchos::rcp( new Epetra_Map( n, n, &newDomainElements[0], OldMap.IndexBase(), OldMap.Comm() ) );

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
