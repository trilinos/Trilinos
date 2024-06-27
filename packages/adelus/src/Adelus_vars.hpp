/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef __ADELUS_VARS_HPP__
#define __ADELUS_VARS_HPP__


namespace Adelus {

template <class ScalarType,
          class ExecutionSpace,
          class MemorySpace>
class AdelusHandle {
 public:
  using value_type      = ScalarType;
  using execution_space = ExecutionSpace;
  using memory_space    = MemorySpace;

 private:

  int comm_id;       // communicator id
  MPI_Comm comm;     // communicator that I belong to
  MPI_Comm row_comm; // row sub-communicator that I belong to
  MPI_Comm col_comm; // column sub-communicator that I belong to

  int myrank;        // process id information

  int nrows_matrix;  // number of rows in the matrix
  int ncols_matrix;  // number of cols in the matrix

  int nprocs_cube;   // num of procs in the allocated cube
  int nprocs_row;    // num of procs to which a row is assigned
  int nprocs_col;    // num of procs to which a col is assigned

  int my_first_row;  // proc position in a row
  int my_first_col;  // proc position in a col

  int my_rows;       // num of rows I own
  int my_cols;       // num of cols I own

  int nrhs;          // number of right hand sides in the matrix
  int my_rhs;        // number of right hand sides that I own

  int blksz;         // block size for matrix update (matrix-matrix multiply)
                     // (e.g. blksz = 128 for GPU, or blksz = 96 for CPU)

  int myrow;         // process id in the col_comm
  int mycol;         // process id in the row_comm



 public:
  AdelusHandle( const int comm_id_,
                MPI_Comm comm_,
                const int matrix_size_,
                const int num_procsr_,
                const int num_rhs_,
                const int blksz_ = 128 )
      : comm_id(comm_id_),
        comm(comm_),
        nrows_matrix(matrix_size_),
        ncols_matrix(matrix_size_),
        nprocs_row(num_procsr_),
        nrhs(num_rhs_),
        blksz(blksz_) {
    // Determine who I am (myrank) and the total number of processes (nprocs_cube)
    MPI_Comm_size(comm, &nprocs_cube);
    MPI_Comm_rank(comm, &myrank);
    nprocs_col = nprocs_cube/nprocs_row;

    // Set up communicators for rows and columns
    mycol = myrank % nprocs_row;
    myrow = myrank / nprocs_row;

    MPI_Comm_split(comm, myrow, mycol, &row_comm);

    MPI_Comm_split(comm, mycol, myrow, &col_comm);

    // Distribution for the matrix on myrank
    my_first_col = myrank % nprocs_row;
    my_first_row = myrank / nprocs_row;

    my_rows = nrows_matrix / nprocs_col;
    if (my_first_row < nrows_matrix % nprocs_col) my_rows++;
    my_cols = ncols_matrix / nprocs_row;
    if (my_first_col < ncols_matrix % nprocs_row) my_cols++;

    // Distribution for the rhs on myrank
    my_rhs = nrhs / nprocs_row;
    if (my_first_col < nrhs % nprocs_row) my_rhs++;
  }

  ~AdelusHandle(){}

  KOKKOS_INLINE_FUNCTION
  int get_comm_id() const { return comm_id; }

  KOKKOS_INLINE_FUNCTION
  MPI_Comm get_comm() const { return comm; }

  KOKKOS_INLINE_FUNCTION
  MPI_Comm get_row_comm() const { return row_comm; }

  KOKKOS_INLINE_FUNCTION
  MPI_Comm get_col_comm() const { return col_comm; }

  KOKKOS_INLINE_FUNCTION
  int get_myrank() const { return myrank; }

  KOKKOS_INLINE_FUNCTION
  int get_myrow() const { return myrow; }

  KOKKOS_INLINE_FUNCTION
  int get_mycol() const { return mycol; }

  KOKKOS_INLINE_FUNCTION
  int get_nprocs_cube() const { return nprocs_cube; }

  KOKKOS_INLINE_FUNCTION
  int get_nprocs_row() const { return nprocs_row; }

  KOKKOS_INLINE_FUNCTION
  int get_nprocs_col() const { return nprocs_col; }

  KOKKOS_INLINE_FUNCTION
  int get_nrows_matrix() const { return nrows_matrix; }

  KOKKOS_INLINE_FUNCTION
  int get_ncols_matrix() const { return ncols_matrix; }

  KOKKOS_INLINE_FUNCTION
  int get_my_first_row() const { return my_first_row; }

  KOKKOS_INLINE_FUNCTION
  int get_my_first_col() const { return my_first_col; }

  KOKKOS_INLINE_FUNCTION
  int get_my_rows() const { return my_rows; }

  KOKKOS_INLINE_FUNCTION
  int get_my_cols() const { return my_cols; }

  KOKKOS_INLINE_FUNCTION
  int get_nrhs() const { return nrhs; }

  KOKKOS_INLINE_FUNCTION
  int get_my_rhs() const { return my_rhs; }

  KOKKOS_INLINE_FUNCTION
  int get_blksz() const { return blksz; }
};

}//namespace Adelus

#endif
