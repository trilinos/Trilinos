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

#include <EpetraExt_AMD_CrsGraph.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <vector>

extern "C" {
#include <amd.h>
}

namespace EpetraExt {

CrsGraph_AMD::
~CrsGraph_AMD()
{
  if( NewMap_ ) delete NewMap_;

  if( NewGraph_ ) delete NewGraph_;
}

CrsGraph_AMD::NewTypeRef
CrsGraph_AMD::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int n = orig.NumMyRows();
  int nnz = orig.NumMyNonzeros();

  //create std CRS format
  vector<int> ia(n+1,0);
  vector<int> ja(nnz);
  int cnt;
  for( int i = 0; i < n; ++i )
  {
    int * tmpP = &ja[ia[i]];
    orig.ExtractMyRowCopy( i, nnz-ia[i], cnt, tmpP );
    ia[i+1] = ia[i] + cnt;
  }

  //trim down to local only
  vector<int> iat(n+1);
  vector<int> jat(nnz);
  int loc = 0;
  for( int i = 0; i < n; ++i )
  {
    iat[i] = loc;
    for( int j = ia[i]; j < ia[i+1]; ++j )
    {
      if( ja[j] < n )
        jat[loc++] = ja[j];
      else
	break;
    }
  }
  iat[n] = loc;


  if( verbose_ )
  {
    cout << "Orig Graph\n";
    cout << orig << endl;
    cout << "-----------------------------------------\n";
    cout << "CRS Format Graph\n";
    cout << "-----------------------------------------\n";
    for( int i = 0; i < n; ++i )
    {
      cout << i << ": " << iat[i+1] << ": ";
      for( int j = iat[i]; j<iat[i+1]; ++j )
        cout << " " << jat[j];
      cout << endl;
    }
    cout << "-----------------------------------------\n";
  }

  vector<int> perm(n);
  vector<double> info(AMD_INFO);

  amd_order( n, &iat[0], &jat[0], &perm[0], NULL, &info[0] ); 

  if( info[AMD_STATUS] == AMD_INVALID )
    cout << "AMD ORDERING: Invalid!!!!\n";

  if( verbose_ )
  {
    cout << "-----------------------------------------\n";
    cout << "AMD Output\n";
    cout << "-----------------------------------------\n";
    cout << "STATUS: " << info[AMD_STATUS] << endl;
    cout << "SYMM: " << info[AMD_SYMMETRY] << endl;
    cout << "N: " << info[AMD_N] << endl;
    cout << "NZ: " << info[AMD_NZ] << endl;
    cout << "SYMM: " << info[AMD_SYMMETRY] << endl;
    cout << "NZDIAG: " << info[AMD_NZDIAG] << endl;
    cout << "NZ A+At: " << info[AMD_NZ_A_PLUS_AT] << endl;
    cout << "NDENSE: " << info[AMD_SYMMETRY] << endl;
    cout << "Perm\n";
    for( int i = 0; i<n; ++i )
      cout << perm[i] << endl;
    cout << "-----------------------------------------\n";
  }

  //Generate New Domain and Range Maps
  //for now, assume they start out as identical
  const Epetra_BlockMap & OldMap = orig.RowMap();
  int nG = orig.NumGlobalRows();

  vector<int> newElements( n );
  for( int i = 0; i < n; ++i )
    newElements[i] = OldMap.GID( perm[i] );

  NewMap_ = new Epetra_Map( nG, n, &newElements[0], OldMap.IndexBase(), OldMap.Comm() );

  if( verbose_ )
  {
    cout << "Old Map\n";
    cout << OldMap << endl;
    cout << "New Map\n";
    cout << *NewMap_ << endl;
  }

  //Generate New Graph
  NewGraph_ = new Epetra_CrsGraph( Copy, *NewMap_, 0 );
  Epetra_Import Importer( *NewMap_, OldMap );
  NewGraph_->Import( orig, Importer, Insert );
  NewGraph_->TransformToLocal();

  if( verbose_ )
  {
    cout << "New CrsGraph\n";
    cout << *NewGraph_ << endl;
  }

  newObj_ = NewGraph_;

  return *NewGraph_;
}

} //namespace EpetraExt
