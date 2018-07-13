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
void TTrilinos_Util_CountMatrixMarket(
    const char *data_file,
    std::vector<int> &non_zeros,
    int_type &N_rows,
    int_type &nnz,
    const Epetra_Comm  &comm
    )
{
  FILE *in_file ;

  N_rows = 0 ;
  nnz = 0 ;
  int_type vecsize = non_zeros.size();
  assert( vecsize == 0 ) ;
  const int BUFSIZE = 800 ;
  char buffer[BUFSIZE] ;
  bool first_off_diag = true ;
  bool upper = false ;

  if(comm.MyPID() == 0)  {
    /* Get information about the array stored in the file specified in the  */
    /* argument list:                                                       */

    in_file = fopen( data_file, "r");
    if (in_file == NULL)
      throw "TTrilinos_Util_CountMatrixMarket: I/O error.";

    if (NULL == fgets( buffer, BUFSIZE, in_file ))
      throw "TTrilinos_Util_CountMatrixMarket: I/O error.";
    std::string headerline1 = buffer;
    bool symmetric = (headerline1.find("symmetric") != std::string::npos);

    // Ignore optional comment lines before ignoring the matrix size information.
    bool found_count = false;
    while ( !found_count && fgets( buffer, BUFSIZE, in_file ) ) {
      if ( buffer[0] != '%' ) {
        found_count = true;
      }
    }

    while ( fgets( buffer, BUFSIZE, in_file ) ) {
      int_type i, j;
      float val ;

      // mfh 20 Mar 2015: We use temporaries of the type corresponding
      // to the sscanf format specifiers, in order to avoid compiler
      // warnings about the sscanf output arguments' types not
      // matching their corresponding format specifiers.  This was a
      // harmless warning, because the 'sizeof' branch prevented
      // incorrect execution, but the warning is easy to fix.
      if(sizeof(int) == sizeof(int_type)) {
        int i_int, j_int;
        sscanf( buffer, "%d %d %f", &i_int, &j_int, &val ) ;
        i = static_cast<int_type> (i_int);
        j = static_cast<int_type> (j_int);
      }
      else if(sizeof(long long) == sizeof(int_type)) {
        long long i_ll, j_ll;
        sscanf( buffer, "%lld %lld %f", &i_ll, &j_ll, &val ) ;
        i = static_cast<int_type> (i_ll);
        j = static_cast<int_type> (j_ll);
      }
      else {
        assert(false);
      }
      int_type needvecsize = i;
      if (symmetric) needvecsize = EPETRA_MAX(i,j) ;
      if ( needvecsize >= vecsize ) {
        int_type oldvecsize = vecsize;
        vecsize += EPETRA_MAX((int_type) 1000,needvecsize-vecsize) ;
        non_zeros.resize(vecsize) ;
        for ( int_type ii = oldvecsize; ii < vecsize ; ii++ ) {
          non_zeros[ii] = 0 ;
        }
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
          std::cout << "file not symmetric" << std::endl ;
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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_CountMatrixMarket(
    const char *data_file,
    std::vector<int> &non_zeros,
    int &N_rows, int &nnz,
    const Epetra_Comm &comm
    )
{
  TTrilinos_Util_CountMatrixMarket<int>(data_file, non_zeros, N_rows, nnz, comm);
}

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_CountMatrixMarket(
    const char *data_file,
    std::vector<int> &non_zeros,
    long long &N_rows, long long &nnz,
    const Epetra_Comm &comm
    )
{
  TTrilinos_Util_CountMatrixMarket<long long>(data_file, non_zeros, N_rows, nnz, comm);
}

#endif
