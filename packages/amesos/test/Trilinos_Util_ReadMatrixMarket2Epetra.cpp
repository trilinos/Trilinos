// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Trilinos_Util_CountMatrixMarket.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

int Trilinos_Util_ReadMatrixMarket2Epetra( char *data_file,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact ) {
  FILE *in_file ;
  int N_rows, nnz ; 

  bool diag = false;   // If set to false, this prevents any entries on the diagonal

  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  vector<int> non_zeros;   // Number of non-zeros in each row
#if 1
  Trilinos_Util_CountMatrixMarket( data_file, non_zeros, N_rows, nnz, comm ) ; 
#else
  Trilinos_Util_CountMatrixMarket( data_file, non_zeros, N_rows, nnz ) ; 
#endif

  vector<int> ptrs(N_rows+1) ; // Pointers into inds and vals for the start of each row
  vector<int> inds(nnz);     //  Column Indices
  vector<double> vals(nnz);  //  Matrix values
  vector<int> iptrs ;        //  Current pointers into inds and vals for each row

  if(comm.MyPID() == 0)  { 
    //  ptrs, inds and vals together constitute a compressed row storage of the matrix.


    in_file = fopen( data_file, "r");
    assert (in_file != NULL) ;  // Checked in T_U_CountMatrixMarket() 

    ptrs[0] = 0 ; 
    for( int i=0; i< N_rows; i++ ) { 
      ptrs[i+1] = ptrs[i] + non_zeros[i]; 
    }

    iptrs = ptrs ; //  Current pointers into inds and vals for each row

    fgets( buffer, BUFSIZE, in_file ) ;  // Pick symmetry info off of this string 
    bool symmetric = false ; 
    string headerline1 = buffer;
    if ( headerline1.find("symmetric") < headerline1.size() ) 
      symmetric = true; 

    //  http://www.yolinux.com/TUTORIALS/LinuxTutorialC++StringClass.html 
    //  says that Var.find() is the function that we want.
    fgets( buffer, BUFSIZE, in_file ) ;

    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int i, j; 
      double val ; 
      i = -13 ;   // Check for blank lines 
      sscanf( buffer, "%d %d %lg", &i, &j, &val ) ; 
      assert( i != -13) ; 
#if 0
      if ( true ) cout << " i,j = " << i <<  " " << j << " val = " << val << endl ; 
#endif
      if ( diag || i != j ) { 
	int iptr = iptrs[i-1] ; 
	iptrs[i-1]++ ;
	vals[iptr] = val ; 
	inds[iptr] = j-1 ; 
	//
	//  If this is a symmetric matrix, we need to enter the entry 
	//  for the other triangular half
	//
	if (symmetric && i != j ) {
	  iptr = iptrs[j-1] ; 
	  iptrs[j-1]++;
	  vals[iptr] = val ; 
	  inds[iptr] = i-1 ; 
	}
      }
    }
    fclose(in_file);


    if ( diag ) { 
      for (int i=0; i<N_rows; i++)
	assert( iptrs[i] == ptrs[i+1] ) ; 
    }
  }


  int nlocal = 0;
  if (comm.MyPID()==0) nlocal = N_rows;
  map = new Epetra_Map(N_rows, nlocal, 0, comm); // Create map with all elements on PE 0
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0
  if (comm.MyPID()==0)
    for (int i=0; i<N_rows; i++) {
      A->InsertGlobalValues(i, iptrs[i]-ptrs[i], &vals[ptrs[i]], &inds[ptrs[i]]);
    }
  A->TransformToLocal();

#if 0
  Epetra_Vector diagA(*map);
  Epetra_Vector Row1056(*map);

  A->ExtractDiagonalCopy( diagA ) ; 
  A->ExtractDiagonalCopy( Row1056 ) ; 

  comm.Barrier();
  comm.Barrier();

  vector <int> r1055ptrs[10000];
  vector <double> r1055vals[10000];

  diagA.PutScalar( 0.0 ) ; 
  diagA.ReplaceGlobalValue( 1055, 0, 1.0 ) ; //  Will generate an error code of 1 on processes which do not own entry 1055 

  A->Multiply( false, diagA, Row1056 ); 
  cout << " Row1056 = " << Row1056  << endl ; 

  //  cout << " A = " << *A  << endl ; 

#endif
  assert( diag || A->NoDiagonal() ) ;

  vector<double> hbx(N_rows);
  vector<double> hbb(N_rows);
  vector<double> hbxexact(N_rows);

  
  x = new Epetra_Vector(Copy, *map, &hbx[0]);
  b = new Epetra_Vector(Copy, *map, &hbb[0]);
  xexact = new Epetra_Vector(Copy, *map, &hbxexact[0]);

  EPETRA_CHK_ERR( x->PutScalar( 0.0 ) );
  EPETRA_CHK_ERR( xexact->Random( ) ) ; 
  EPETRA_CHK_ERR( A->Multiply( false, *xexact, *b ) ); 

  return 0;
}
