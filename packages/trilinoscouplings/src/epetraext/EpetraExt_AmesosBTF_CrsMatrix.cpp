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
#include <EpetraExt_AmesosBTF_CrsMatrix.h>

#include <Epetra_Import.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <amesos_btf_decl.h>

using std::vector;
using std::cout;
using std::endl;

namespace EpetraExt {

AmesosBTF_CrsMatrix::
~AmesosBTF_CrsMatrix()
{
}

AmesosBTF_CrsMatrix::NewTypeRef
AmesosBTF_CrsMatrix::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;
  const Epetra_BlockMap & OldRowMap = orig.RowMap();
  const Epetra_BlockMap & OldColMap = orig.ColMap();
  
  // Check if the matrix is on one processor.
  int myMatProc = -1, matProc = -1;
  int myPID = orig.Comm().MyPID();
  for (int proc=0; proc<orig.Comm().NumProc(); proc++) 
  {
    if (orig.NumGlobalNonzeros() == orig.NumMyNonzeros())
      myMatProc = myPID;
  }
  orig.Comm().MaxAll( &myMatProc, &matProc, 1 );
  
  if( orig.RowMap().DistributedGlobal() && matProc == -1)
    { cout << "FAIL for Global!\n"; abort(); }
  if( orig.IndicesAreGlobal() && matProc == -1)
    { cout << "FAIL for Global Indices!\n"; abort(); }
 
  int nGlobal = orig.NumGlobalRows(); 
  int n = orig.NumMyRows();
  int nnz = orig.NumMyNonzeros();
  
  if( debug_ )
  {
    cout << "Orig Matrix:\n";
    cout << orig << endl;
  }

  // Create std CRS format (without elements above the threshold)
  vector<int> ia(n+1,0);
  int maxEntries = orig.MaxNumEntries();
  vector<int> ja(nnz), ja_tmp(nnz);
  vector<double> jva_tmp(maxEntries);
  int cnt;

  Epetra_CrsGraph strippedGraph( Copy, OldRowMap, OldColMap, 0 );

  for( int i = 0; i < n; ++i )
  {
    orig.ExtractMyRowCopy( i, maxEntries, cnt, &jva_tmp[0], &ja_tmp[0] );
    ia[i+1] = ia[i];
    for( int j = 0; j < cnt; ++j )
      if( fabs(jva_tmp[j]) > threshold_ )
        ja[ ia[i+1]++ ] = ja_tmp[j];

    int new_cnt = ia[i+1] - ia[i];
    strippedGraph.InsertMyIndices( i, new_cnt, &ja[ ia[i] ] );
  }
  nnz = ia[n];
  strippedGraph.FillComplete();
  
  if( debug_ )
  {
    cout << "Stripped Graph\n";
    cout << strippedGraph;
  }

  // Compute the BTF permutation only on the processor that has the graph.
  if ( matProc == myPID ) {
    
    if( debug_ )
      {
	cout << "-----------------------------------------\n";
	cout << "CRS Format Graph (stripped) \n";
	cout << "-----------------------------------------\n";
	for( int i = 0; i < n; ++i )
	  {
	    cout << ia[i] << " - " << ia[i+1] << " : ";
	    for( int j = ia[i]; j<ia[i+1]; ++j )
	      cout << " " << ja[j];
	    cout << endl;
	  }
	cout << "-----------------------------------------\n";
      }
  
    // Transpose the graph, not the values
    int j=0, next=0;
    vector<int> ia_tmp(n+1,0);

    // Compute row lengths
    for (int i = 0; i < n; i++)
        for (int k = ia[i]; k < ia[i+1]; k++)
            ++ia_tmp[ ja[k]+1 ];

    // Compute pointers from row lengths
    ia_tmp[0] = 0;
    for (int i = 0; i < n; i++)
        ia_tmp[i+1] += ia_tmp[i];

    // Copy over indices
    for (int i = 0; i < n; i++) {
        for (int k = ia[i]; k < ia[i+1]; k++) {
            j = ja[k];
            next = ia_tmp[j];
            ja_tmp[next] = i;
            ia_tmp[j] = next + 1;
        }
    }

    // Reshift ia_tmp
    for (int i=n-1; i >= 0; i--) ia_tmp[i+1] = ia_tmp[i];
    ia_tmp[0] = 0;

    // Transformation information
    int numMatch = 0;       // number of nonzeros on diagonal after permutation.
    double maxWork =  0.0;  // no limit on how much work to perform in max-trans.
    double workPerf = 0.0;  // how much work was performed in max-trans.
    
    // Create a work vector for the BTF code.
    vector<int> work(5*n);
    
    // Storage for the row and column permutations.
    vector<int> rowperm(n);
    vector<int> colperm(n);
    vector<int> blockptr(n+1);
    
    // NOTE:  The permutations are sent in backwards since the matrix is transposed.
    // On output, rowperm and colperm are the row and column permutations of A, where 
    // i = BTF_UNFLIP(rowperm[k]) if row i of A is the kth row of P*A*Q, and j = colperm[k] 
    // if column j of A is the kth column of P*A*Q.  If rowperm[k] < 0, then the 
    // (k,k)th entry in P*A*Q is structurally zero.
    
    numBlocks_ = amesos_btf_order( n, &ia_tmp[0], &ja_tmp[0], maxWork, &workPerf,
			    &rowperm[0], &colperm[0], &blockptr[0], 
			    &numMatch, &work[0] );
    
    // Reverse ordering of permutation to get upper triangular form, if necessary.
    rowPerm_.resize( n );
    colPerm_.resize( n ); 
    blockptr.resize( numBlocks_+1 );
    blockPtr_.resize( numBlocks_+1 );
    if (!upperTri_) {
      for( int i = 0; i < n; ++i )
	{
	  rowPerm_[i] = BTF_UNFLIP(rowperm[(n-1)-i]);
	  colPerm_[i] = colperm[(n-1)-i];
	}
      for( int i = 0; i < numBlocks_+1; ++i ) 
	{
	  blockPtr_[i] = n - blockptr[numBlocks_-i];
	}
    }
    else {
      colPerm_ = colperm;
      blockPtr_ = blockptr;
      for( int i = 0; i < n; ++i )
	{
	  rowPerm_[i] = BTF_UNFLIP(rowperm[i]);
	}
    }
    
    if( debug_ ) {
      cout << "-----------------------------------------\n";
      cout << "BTF Output (n = " << n << ")\n";
      cout << "-----------------------------------------\n";
      cout << "Num Blocks: " << numBlocks_ << endl;
      cout << "Num NNZ Diags: " << numMatch << endl;
      cout << "RowPerm and ColPerm \n";
      for( int i = 0; i<n; ++i )
	cout << rowPerm_[i] << "\t" << colPerm_[i] << endl;
      cout << "-----------------------------------------\n";
    }  
  }

  // Broadcast the BTF permutation information to all processors.
  rowPerm_.resize( nGlobal );
  colPerm_.resize( nGlobal );

  orig.Comm().Broadcast(&rowPerm_[0], nGlobal, matProc);
  orig.Comm().Broadcast(&colPerm_[0], nGlobal, matProc);
  orig.Comm().Broadcast(&numBlocks_, 1, matProc);

  blockPtr_.resize( numBlocks_+1 );
  orig.Comm().Broadcast(&blockPtr_[0], numBlocks_+1, matProc);
  
  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  vector<int> myElements( n );
  OldRowMap.MyGlobalElements( &myElements[0] );
  
  vector<int> newDomainElements( n );
  vector<int> newRangeElements( n );
  for( int i = 0; i < n; ++i )
  {
    newRangeElements[ i ] = myElements[ rowPerm_[i] ];
    newDomainElements[ i ] = myElements[ colPerm_[i] ];
  }

  NewRowMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, n, &newRangeElements[0], OldRowMap.IndexBase(), OldRowMap.Comm() ) );
  NewColMap_ = Teuchos::rcp( new Epetra_Map( nGlobal, n, &newDomainElements[0], OldColMap.IndexBase(), OldColMap.Comm() ) );

  if( debug_ )
  {
    cout << "New Row Map\n";
    cout << *NewRowMap_ << endl;
    cout << "New Col Map\n";
    cout << *NewColMap_ << endl;
  }

  //Generate New Graph
  NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 ) );
  Importer_ = Teuchos::rcp( new Epetra_Import( *NewRowMap_, OldRowMap ) );
  NewGraph_->Import( strippedGraph, *Importer_, Insert );
  NewGraph_->FillComplete();

  if( debug_ )
  {
    cout << "NewGraph\n";
    cout << *NewGraph_;
  }

  NewMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *NewGraph_ ) );
  NewMatrix_->Import( orig, *Importer_, Insert );
  NewMatrix_->FillComplete();

  if( debug_ )
  {
    cout << "New CrsMatrix\n";
    cout << *NewMatrix_ << endl;
  }

  newObj_ = &*NewMatrix_;

  return *NewMatrix_;
}

bool
AmesosBTF_CrsMatrix::
fwd()
{
  NewMatrix_->Import( *origObj_, *Importer_, Insert );
  return true;
}

bool
AmesosBTF_CrsMatrix::
rvs()
{
  origObj_->Export( *NewMatrix_, *Importer_, Insert );
  return true;
}

} //namespace EpetraExt
