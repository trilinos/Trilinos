/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_MATRIX_DATA_H
#define EPETRA_MATRIX_DATA_H

#include "Epetra_ConfigDefs.h"

class Epetra_CrsMatrix;

namespace epetra_test {

/** matrix_data is a very simple data source to be used for filling test matrices.
    It is serial; the intent is that a test program declares this class to be
    of full (global) size on each processor, then fills the local portion of
    the test matrix from the appropriate section of the data in this class.
*/
class EPETRA_LIB_DLL_EXPORT matrix_data {
 public:
  matrix_data(int num_rows, int* rowlengths, int blocksize=1);
  matrix_data(int num_rows, int num_cols, int num_off_diagonals, int blocksize);
  matrix_data(int num_quad_elements, int num_dof_per_node,
	      bool make_numerically_nonsymmetric=false);
  virtual ~matrix_data();

  int numrows() { return(numrows_); }
  int numcols() { return(numcols_); }
  int blocksize() { return(blocksize_); }
  int* rows() { return(rows_); }
  int* rowlengths() { return(rowlengths_); }

  int** colindices() { return(colindices_); }
  double** coefs() { return(coefs_); }
  double* coefs(int row, int col);

  /** The portion of this matrix_data object's data that corresponds to
      the locally-owned rows of A, will be copied into A. A.FillComplete()
      will NOT be called.
  */
  //Not used.
  //void copy_local_data_to_matrix(Epetra_CrsMatrix& A);

  /** Compare the local rows of A to the corresponding rows of this
      matrix_data object's data.
  */ 
  bool compare_local_data(const Epetra_CrsMatrix& A);

 private:
  int numrows_;
  int numcols_;
  int* rows_;
  int* rowlengths_;
  int blocksize_;

  int** colindices_;
  double** coefs_;

  matrix_data(const matrix_data & data);
  matrix_data & operator=(const matrix_data & data);

};//class matrix_data

}//namespace epetra_test

#endif

