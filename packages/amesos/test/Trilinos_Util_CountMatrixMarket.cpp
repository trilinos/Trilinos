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
#include "Epetra_Object.h"
#include "Epetra_Comm.h"

//
//  This code reads a file which contains a sparse matrix
//  in triplet (i,j,val) form, counting the number of non-zero elements 
//  in each row.  
//
//  Returns:  N_rows and nnz replicated across all processes
//
void Trilinos_Util_CountMatrixMarket( const char *data_file, 
				      vector<int> &non_zeros,
				      int &N_rows, int &nnz, 
				      const Epetra_Comm  &comm) { 

  FILE *in_file ;
  
  N_rows = 0 ; 
  nnz = 0 ; 
  int vecsize = non_zeros.size(); 
  assert( vecsize == 0 ) ; 
  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  bool first_off_diag = true ; 
  bool upper ;

  if(comm.MyPID() == 0)  { 
    /* Get information about the array stored in the file specified in the  */
    /* argument list:                                                       */
    
    in_file = fopen( data_file, "r");
    if (in_file == NULL)
      {
	printf("Error: Cannot open file: %s\n",data_file);
	exit(1);
      }
    
    fgets( buffer, BUFSIZE, in_file ) ;
    bool symmetric = false ; 
    string headerline1 = buffer;
    if ( headerline1.find("symmetric") < headerline1.size() ) symmetric = true; 

    cout << " string::npos " << string::npos  << endl ; 
    cout << " headerline1.find(symmetric)" << headerline1.find("symmetric")  << endl ; 
    cout << "  headerline1.size() " <<  headerline1.size()  << endl ; 

    cout << ( symmetric?" symmetic ":" not symmetic " ) << endl ; 

    fgets( buffer, BUFSIZE, in_file ) ;
    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int i, j; 
      double val ; 
      sscanf( buffer, "%d %d %f", &i, &j, &val ) ; 
      int needvecsize = i;
      if (symmetric) needvecsize = EPETRA_MAX(i,j) ;
      if ( needvecsize >= vecsize ) {
	int oldvecsize = vecsize; 
	vecsize += EPETRA_MAX(1000,needvecsize-vecsize) ; 
	non_zeros.resize(vecsize) ; 
        for ( int i= oldvecsize; i < vecsize ; i++ ) non_zeros[i] = 0 ; 
      }
      N_rows = EPETRA_MAX( N_rows, i ) ; 
      if (symmetric) N_rows = EPETRA_MAX( N_rows, j ) ; 
      non_zeros[i-1]++ ; 
      nnz++; 
      if ( symmetric && i != j ) {
	if ( first_off_diag ) { 
	  upper = j > i ; 
	  first_off_diag = false ; 
	}
	if ( ( j > i && ! upper ) || ( i > j && upper ) ) { 
	  cout << "file not symmetric" << endl ; 
	  exit(1) ; 
	}
	non_zeros[j-1]++ ; 
	nnz++; 
      } 
    } 
    fclose(in_file);
  }
  comm.Broadcast( &N_rows, 1, 0 );
  comm.Broadcast( &nnz, 1, 0 );
  return;
}
