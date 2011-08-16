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

#include <amesos_btf_decl.h>
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

  if( orig.RowMap().DistributedGlobal() )
  { cout << "FAIL for Global!\n"; abort(); }
  if( orig.IndicesAreGlobal() )
  { cout << "FAIL for Global Indices!\n"; abort(); }

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
    cout << "-----------------------------------------\n";
    cout << "CRS Format Graph\n";
    cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
      {
	cout << Ap[i] << " - " << Ap[i+1] << " : ";
	for( int j = Ap[i]; j<Ap[i+1]; ++j )
	  cout << " " << Ai[j];
	cout << endl;
      }
    cout << "-----------------------------------------\n";
  }

  // Control and information vector
  double control[AMD_CONTROL];
  double info[AMD_INFO];

  // Storage for the permutation.
  perm_.resize( n );

  // Call AMD permutation.

  int ret = amesos_amd_order( n, &Ap[0], &Ai[0], &perm_[0], control, info );

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
    cout << "New Row Map\n";
    cout << *NewRowMap_ << endl;
    cout << "New Col Map\n";
    cout << *NewColMap_ << endl;
  }

  //Generate New Graph
  NewGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 ) );
  Epetra_Import Importer( *NewRowMap_, OldMap );
  NewGraph_->Import( orig, Importer, Insert );
  NewGraph_->FillComplete();

  if (verbose_) 
  {
    cout << "New CrsGraph\n";
    cout << *NewGraph_ << endl;
  }

  newObj_ = &*NewGraph_;

  return *NewGraph_;
}

} //namespace EpetraExt
