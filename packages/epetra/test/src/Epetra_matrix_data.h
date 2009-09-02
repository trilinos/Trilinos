/*
************************************************************************

              Epetra: Linear Algebra Services Package
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov)

************************************************************************
*/

#ifndef EPETRA_MATRIX_DATA_H
#define EPETRA_MATRIX_DATA_H

class Epetra_CrsMatrix;

namespace epetra_test {

/** matrix_data is a very simple data source to be used for filling test matrices.
    It is serial; the intent is that a test program declares this class to be
    of full (global) size on each processor, then fills the local portion of
    the test matrix from the appropriate section of the data in this class.
*/
class matrix_data {
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
  void copy_local_data_to_matrix(Epetra_CrsMatrix& A);

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

