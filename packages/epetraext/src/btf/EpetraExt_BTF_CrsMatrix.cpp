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
#include <EpetraExt_BTF_CrsMatrix.h>

#include <EpetraExt_Transpose_CrsGraph.h>

#include <Epetra_Import.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <vector>

using std::vector;

#define MATTRANS_F77 F77_FUNC(mattrans,MATTRANS)
#define GENBTF_F77   F77_FUNC(genbtf,GENBTF)
                                                                                                  
extern "C" {
extern void MATTRANS_F77( int*, int*, int*, int*, int*, int* );
extern void GENBTF_F77( int*, int*, int*, int*, int*, int*, int*, int*, int*,
                     int*, int*, int*, int*, int*, int*, int*, int*, int*,
                     int*, int*, int* );
}

namespace EpetraExt {

CrsMatrix_BTF::
~CrsMatrix_BTF()
{
  if( NewRowMap_ ) delete NewRowMap_;
  if( NewColMap_ ) delete NewColMap_;

  if( Importer_ ) delete Importer_;

  if( NewMatrix_ ) delete NewMatrix_;
  if( NewGraph_ ) delete NewGraph_;
}

CrsMatrix_BTF::NewTypeRef
CrsMatrix_BTF::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  if( orig.RowMap().DistributedGlobal() )
  { cout << "FAIL for Global!\n"; abort(); }
  if( orig.IndicesAreGlobal() )
  { cout << "FAIL for Global Indices!\n"; abort(); }

  int n = orig.NumMyRows();
  int nnz = orig.NumMyNonzeros();

  if( verbose_ )
  {
    cout << "Orig Matrix:\n";
    cout << orig << endl;
  }

  //create std CRS format
  //also create graph without zero elements
  vector<int> ia(n+1,0);
  int maxEntries = orig.MaxNumEntries();
  vector<int> ja_tmp(maxEntries);
  vector<double> jva_tmp(maxEntries);
  vector<int> ja(nnz);
  int cnt;

  const Epetra_BlockMap & OldRowMap = orig.RowMap();
  const Epetra_BlockMap & OldColMap = orig.ColMap();
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

  if( verbose_ )
  {
    cout << "Stripped Graph\n";
    cout << strippedGraph;
  }

  vector<int> iat(n+1,0);
  vector<int> jat(nnz);
  for( int i = 0; i < n; ++i )
    for( int j = ia[i]; j < ia[i+1]; ++j )
      ++iat[ ja[j]+1 ];
  for( int i = 0; i < n; ++i )
    iat[i+1] += iat[i];
  for( int i = 0; i < n; ++i )
    for( int j = ia[i]; j < ia[i+1]; ++j )
      jat[ iat[ ja[j] ]++ ] = i;
  for( int i = 0; i < n; ++i )
    iat[n-i] = iat[n-i-1];
  iat[0] = 0;

  //convert to Fortran indexing
  for( int i = 0; i < n+1; ++i ) ++ia[i];
  for( int i = 0; i < nnz; ++i ) ++ja[i];
  for( int i = 0; i < n+1; ++i ) ++iat[i];
  for( int i = 0; i < nnz; ++i ) ++jat[i];

  if( verbose_ )
  {
    cout << "-----------------------------------------\n";
    cout << "CRS Format Graph\n";
    cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
    {
      cout << i+1 << ": " << ia[i+1] << ": ";
      for( int j = ia[i]-1; j<ia[i+1]-1; ++j )
        cout << " " << ja[j];
      cout << endl;
    }
    cout << "-----------------------------------------\n";
  }

/*
  vector<int> iat(n+1);
  vector<int> jat(nnz);
  int * jaf = &ja[0];
  int * iaf = &ia[0];
  int * jatf = &jat[0];
  int * iatf = &iat[0];
  MATTRANS_F77( &n, &n, &ja[0], &ia[0], &jat[0], &iat[0] );
*/
    
  if( verbose_ )
  {
    cout << "-----------------------------------------\n";
    cout << "CCS Format Graph\n";
    cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
    {
      cout << i+1 << ": " << iat[i+1] << ": ";
      for( int j = iat[i]-1; j<iat[i+1]-1; ++j )
        cout << " " << jat[j];
      cout << endl;
    }
    cout << "-----------------------------------------\n";
  }

  vector<int> w(10*n);

  vector<int> rowperm(n);
  vector<int> colperm(n);

  //horizontal block
  int nhrows, nhcols, hrzcmp;
  //square block
  int nsrows, sqcmpn;
  //vertial block
  int nvrows, nvcols, vrtcmp;

  vector<int> rcmstr(n+1);
  vector<int> ccmstr(n+1);

  int msglvl = 0;
  int output = 6;

  GENBTF_F77( &n, &n, &iat[0], &jat[0], &ia[0], &ja[0], &w[0],
          &rowperm[0], &colperm[0], &nhrows, &nhcols,
          &hrzcmp, &nsrows, &sqcmpn, &nvrows, &nvcols, &vrtcmp,
          &rcmstr[0], &ccmstr[0], &msglvl, &output );

  //convert back to C indexing
  for( int i = 0; i < n; ++i )
  {
    --rowperm[i];
    --colperm[i];
  }
  for( int i = 0; (i<n+1) && (rcmstr[i]!=n+1); ++i )
  {
    --rcmstr[i];
    --ccmstr[i];
  }

  if( verbose_ )
  {
    cout << "-----------------------------------------\n";
    cout << "BTF Output\n";
    cout << "-----------------------------------------\n";
    cout << "RowPerm and ColPerm\n";
    for( int i = 0; i<n; ++i )
      cout << rowperm[i] << "\t" << colperm[i] << endl;
    if( hrzcmp )
    {
      cout << "Num Horizontal: Rows, Cols, Comps\n";
      cout << nhrows << "\t" << nhcols << "\t" << hrzcmp << endl;
    }
    cout << "Num Square: Rows, Comps\n";
    cout << nsrows << "\t" << sqcmpn << endl;
    if( vrtcmp )
    {
      cout << "Num Vertical: Rows, Cols, Comps\n";
      cout << nvrows << "\t" << nvcols << "\t" << vrtcmp << endl;
    }
    cout << "Row, Col of upper left pt in blocks\n";
    for( int i = 0; (i<n+1) && (rcmstr[i]!=n+1); ++i )
      cout << i << " " << rcmstr[i] << " " << ccmstr[i] << endl;
    cout << "-----------------------------------------\n";
  }

  if( hrzcmp || vrtcmp )
  {
    cout << "FAILED! hrz cmp's:" << hrzcmp << " vrtcmp: " << vrtcmp << endl;
    exit(0);
  }

  //convert rowperm to OLD->NEW
  //reverse ordering of permutation to get upper triangular
  vector<int> rowperm_t( n );
  vector<int> colperm_t( n );
  for( int i = 0; i < n; ++i )
  {
//    rowperm_t[ rowperm[i] ] = i;
    rowperm_t[i] = rowperm[i];
    colperm_t[i] = colperm[i];
  }

  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  vector<int> myElements( n );
  OldRowMap.MyGlobalElements( &myElements[0] );

  vector<int> newDomainElements( n );
  vector<int> newRangeElements( n );
  for( int i = 0; i < n; ++i )
  {
    newDomainElements[ i ] = myElements[ rowperm_t[i] ];
    newRangeElements[ i ] = myElements[ colperm_t[i] ];
  }

  NewRowMap_ = new Epetra_Map( n, n, &newDomainElements[0], OldRowMap.IndexBase(), OldRowMap.Comm() );
  NewColMap_ = new Epetra_Map( n, n, &newRangeElements[0], OldColMap.IndexBase(), OldColMap.Comm() );

  if( verbose_ )
  {
    cout << "New Row Map\n";
    cout << *NewRowMap_ << endl;
    cout << "New ColMap\n";
    cout << *NewColMap_ << endl;
  }

  //Generate New Graph
  NewGraph_ = new Epetra_CrsGraph( Copy, *NewRowMap_, *NewColMap_, 0 );
  Importer_ = new Epetra_Import( *NewRowMap_, OldRowMap );
  NewGraph_->Import( strippedGraph, *Importer_, Insert );
  NewGraph_->FillComplete();
  if( verbose_ )
  {
    cout << "NewGraph\n";
    cout << *NewGraph_;
  }

  NewMatrix_ = new Epetra_CrsMatrix( Copy, *NewGraph_ );
  NewMatrix_->Import( orig, *Importer_, Insert );
  NewMatrix_->FillComplete();

  if( verbose_ )
  {
    cout << "New CrsMatrix\n";
    cout << *NewMatrix_ << endl;
  }

  newObj_ = NewMatrix_;

  return *NewMatrix_;
}

bool
CrsMatrix_BTF::
fwd()
{
  int ret = NewMatrix_->Import( *origObj_, *Importer_, Insert );
  if (ret<0) return false;
  return true;
}

bool
CrsMatrix_BTF::
rvs()
{
  int ret = origObj_->Export( *NewMatrix_, *Importer_, Insert );
  if (ret<0) return false;
  return true;
}

} //namespace EpetraExt
