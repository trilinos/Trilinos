// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <stdlib.h>
#include <stdio.h>
#include <vector>
// #include "Trilinos_Util_Triples.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"

//
//  This code reads a file which contains a sparse matrix
//  in triplet (i,j,val) form, counting the number of non-zero elements 
//  in each row.  
//
//  Returns:  N_rows and nnz replicated across all processes
//
template<typename int_type>
void Trilinos_Util_CountTriples_internal( const char *data_file, 
				 bool symmetric, 
				 std::vector<int> &non_zeros,
				 int_type &N_rows, int_type &nnz, 
				 const Epetra_Comm  &comm, 
				 bool TimDavisHeader=false,
				 bool ZeroBased=false 
) { 

  FILE *in_file ;

  N_rows = 0 ; 
  nnz = 0 ; 
  int_type vecsize = non_zeros.size(); 
  assert( vecsize == 0 ) ; 
  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  bool first_off_diag = true ; 
  bool upper ;

  int_type num_rows = -1 ; 
  int_type num_cols = -1 ; 
  int_type num_nz = -1 ; 
  int hdr_type = -131313 ; 
    
  if(comm.MyPID() == 0)  { 
    /* Get information about the array stored in the file specified in the  */
    /* argument list:                                                       */
    
    in_file = fopen( data_file, "r");
    if (in_file == NULL)
      {
	printf("Error: Cannot open file: %s\n",data_file);
	exit(1);
      }
    
    if ( TimDavisHeader ) { 
      fgets( buffer, BUFSIZE, in_file );
      if(sizeof(int) == sizeof(int_type))
        sscanf( buffer, "%d %d %d %d", &num_rows, &num_cols, &num_nz, &hdr_type ) ; 
      else if(sizeof(long long) == sizeof(int_type))
        sscanf( buffer, "%lld %lld %lld %d", &num_rows, &num_cols, &num_nz, &hdr_type ) ; 
      else
        assert(false);
      if( hdr_type != 0 ) {
	if ( hdr_type == -131313 ) 
	  printf("Bad Tim Davis header line.  Should have four  values and the fourth must be zero.\n"); 
	else
	  printf("Bad Tim Davis header line.  Fourth value must be zero, but is %d.\n", hdr_type); 
	exit(1);
      }
      if ( num_rows != num_cols ) {
	printf( "Bad Tim Davis header line.  First two values, number of rows and columns must be equal.  We see %d and %d\n" ,
		num_rows, num_cols ) ; 
      }
    } 

    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int_type i, j; 
      float val ; 
      i = -13 ;   // Check for blank lines 
      if(sizeof(int) == sizeof(int_type))
        sscanf( buffer, "%d %d %f", &i, &j, &val ) ; 
      else if(sizeof(long long) == sizeof(int_type))
        sscanf( buffer, "%lld %lld %f", &i, &j, &val ) ; 
      else
        assert(false);
	  
      if ( ZeroBased ) { i++; j++ ; } 
      if ( i > 0 ) { 
	int_type needvecsize = i;
	if (symmetric) needvecsize = EPETRA_MAX(i,j) ;
	if ( needvecsize >= vecsize ) {
	  int_type oldvecsize = vecsize; 
	  vecsize += EPETRA_MAX((int_type) 1000,needvecsize-vecsize) ; 
	  non_zeros.resize(vecsize) ; 
	  for ( int_type i= oldvecsize; i < vecsize ; i++ ) non_zeros[i] = 0 ; 
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
    } 
    fclose(in_file);
  }

  if ( TimDavisHeader && comm.MyPID() == 0)  { 
    if ( num_rows != N_rows ) {
      if(sizeof(int) == sizeof(int_type))
        printf( " Bad Tim Davis Header Line.  The first value should be the number of rows.  We see %d, but the actual number of rows is: %d\n",
	      num_rows, N_rows );
      else if(sizeof(long long) == sizeof(int_type))
        printf( " Bad Tim Davis Header Line.  The first value should be the number of rows.  We see %lld, but the actual number of rows is: %lld\n",
	      num_rows, N_rows );
      else
        assert(false);
    }
    if ( num_nz != nnz ) {
      printf( " Bad Tim Davis Header Line.  The third value should be the number of non-zeros.  We see %d, but the actual number of non-zeros is: %d\n",
	      num_nz, nnz );
    }
  }

  comm.Broadcast( &N_rows, 1, 0 );
  comm.Broadcast( &nnz, 1, 0 );
  return;
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_CountTriples( const char *data_file, 
				 bool symmetric, 
				 std::vector<int> &non_zeros,
				 int &N_rows, int &nnz, 
				 const Epetra_Comm  &comm, 
				 bool TimDavisHeader=false, 
				 bool ZeroBased=false ) {
  Trilinos_Util_CountTriples_internal<int>(data_file, symmetric, non_zeros,
    N_rows, nnz, comm, TimDavisHeader, ZeroBased);
}

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_CountTriples( const char *data_file, 
				 bool symmetric, 
				 std::vector<int> &non_zeros,
				 long long &N_rows, long long &nnz, 
				 const Epetra_Comm  &comm, 
				 bool TimDavisHeader=false, 
				 bool ZeroBased=false ) {
  Trilinos_Util_CountTriples_internal<long long>(data_file, symmetric, non_zeros,
    N_rows, nnz, comm, TimDavisHeader, ZeroBased);
}

#endif
