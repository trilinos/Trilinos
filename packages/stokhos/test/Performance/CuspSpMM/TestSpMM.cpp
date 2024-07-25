// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// CUSP
#include <cusp/array2d.h>
#include <cusp/csr_matrix.h>
#include <cusp/MVmultiply.h>
#include <cusp/gallery/poisson.h>
#include <cusp/detail/timer.h>

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[])
{
  typedef int                 IndexType;
  typedef double              ValueType;
  typedef cusp::device_memory MemorySpace;
  //typedef cusp::row_major     Orientation;

  bool success = true;
  bool verbose = false;
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString("This test performance of block multiply routines.\n");
    IndexType n = 32;
    CLP.setOption("n", &n, "Number of mesh points in the each direction");
    IndexType nrhs_begin = 32;
    CLP.setOption("begin", &nrhs_begin,
                  "Staring number of right-hand-sides");
    IndexType nrhs_end = 512;
    CLP.setOption("end", &nrhs_end,
                  "Ending number of right-hand-sides");
    IndexType nrhs_step = 32;
    CLP.setOption("step", &nrhs_step,
                  "Increment in number of right-hand-sides");
    IndexType nits = 10;
    CLP.setOption("nits", &nits,
                  "Number of multiply iterations");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
    CLP.parse( argc, argv );

    // Set CUDA device
    cudaSetDevice(device_id);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    // create 3D Poisson problem
    cusp::csr_matrix<IndexType, ValueType, MemorySpace> A;
    cusp::gallery::poisson27pt(A, n, n, n);

    std::cout << "nrhs , num_rows , num_entries , row_time , row_gflops , "
              << "col_time , col_gflops" << std::endl;

    for (IndexType nrhs = nrhs_begin; nrhs <= nrhs_end; nrhs += nrhs_step) {

      double flops =
        2.0 * static_cast<double>(A.num_entries) * static_cast<double>(nrhs);

      // test row-major storage
      cusp::array2d<ValueType, MemorySpace, cusp::row_major> x_row(
        A.num_rows, nrhs, 1);
      cusp::array2d<ValueType, MemorySpace, cusp::row_major> y_row(
        A.num_rows, nrhs, 0);

      cusp::detail::timer row_timer;
      row_timer.start();
      for (IndexType iter=0; iter<nits; ++iter) {
        cusp::MVmultiply(A, x_row, y_row);
      }
      cudaDeviceSynchronize();
      double row_time = row_timer.seconds_elapsed() / nits;
      double row_gflops = 1.0e-9 * flops / row_time;

      // test column-major storage
      cusp::array2d<ValueType, MemorySpace, cusp::column_major> x_col(
        A.num_rows, nrhs, 1);
      cusp::array2d<ValueType, MemorySpace, cusp::column_major> y_col(
        A.num_rows, nrhs, 0);

      cusp::detail::timer col_timer;
      col_timer.start();
      for (IndexType iter=0; iter<nits; ++iter) {
        cusp::MVmultiply(A, x_col, y_col);
      }
      cudaDeviceSynchronize();
      double col_time = col_timer.seconds_elapsed() / nits;
      double col_gflops = 1.0e-9 * flops / col_time;

      std::cout << nrhs << " , "
                << A.num_rows << " , "  << A.num_entries << " , "
                << row_time << " , " << row_gflops << " , "
                << col_time << " , " << col_gflops
                << std::endl;

    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
