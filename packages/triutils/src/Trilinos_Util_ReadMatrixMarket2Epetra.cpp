// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
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

  bool diag = true;   // If set to false, this prevents any entries on the diagonal 
                      
  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  vector<int> non_zeros;   // Number of non-zeros in each row
  Trilinos_Util_CountMatrixMarket( data_file, non_zeros, N_rows, nnz, comm ) ; 

  vector<int> ptrs(N_rows+1) ; // Pointers into inds and vals for the start of each row
  vector<int> inds(nnz);     //  Column Indices
  vector<double> vals(nnz);  //  Matrix values
  vector<int> iptrs ;        //  Current pointers into inds and vals for each row

  if(comm.MyPID() == 0)  { 
    //  ptrs, inds and vals together constitute a compressed row storage of the matrix.


    in_file = fopen( data_file, "r");
    assert (in_file != NULL) ;  // Checked in Trilinos_Util_CountMatrixMarket() 

    ptrs[0] = 0 ; 
    for( int i=0; i< N_rows; i++ ) { 
      ptrs[i+1] = ptrs[i] + non_zeros[i]; 
    }

    iptrs = ptrs ; //  Current pointers into inds and vals for each row

    fgets( buffer, BUFSIZE, in_file ) ;  // Pick symmetry info off of this string 
    bool symmetric = false ; 
    string headerline1 = buffer;
#ifdef TFLOP
    // This is a fix for Janus since its string class does not recognize
    // string::npos.  The integer returned when "symmetric" is not found in 
    // the string makes no logical sense and must be a system dependent variable. 
    // However, since find returns the index of the first position in the 
    // string where "symmetric" is found, we will check to see if this 
    // number is larger than the buffer size.  This fix works for Janus, but is
    // not recommended for other platforms since find may also return a 
    // negative number. ( HKT 3/16/2004 )
    if ( headerline1.find("symmetric") < BUFSIZE ) symmetric = true;
#else
    if ( headerline1.find("symmetric") != string::npos) symmetric = true; 
#endif
    fgets( buffer, BUFSIZE, in_file ) ;

    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int i, j; 
      double val ; 
      i = -13 ;   // Check for blank lines 
      sscanf( buffer, "%d %d %lg", &i, &j, &val ) ; 
      assert( i != -13) ; 
      if ( diag || i != j ) { 
	//	if ( i == j && i == 1 ) val *= 1.0001 ;
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

  Epetra_Vector diagA(*map);

  A->ExtractDiagonalCopy( diagA ) ; 

  assert( diag || A->NoDiagonal() ) ;

  vector<double> xyz(N_rows);
  
  x = new Epetra_Vector(Copy, *map, &xyz[0]);
  b = new Epetra_Vector(Copy, *map, &xyz[0]);
  xexact = new Epetra_Vector(Copy, *map, &xyz[0]);

  EPETRA_CHK_ERR( x->PutScalar( 0.0 ) );
  EPETRA_CHK_ERR( xexact->Random( ) ) ; 
  EPETRA_CHK_ERR( A->Multiply( false, *xexact, *b ) ); 

  assert( map->SameAs(A->RowMap()) ) ;
 
  return 0;
}
