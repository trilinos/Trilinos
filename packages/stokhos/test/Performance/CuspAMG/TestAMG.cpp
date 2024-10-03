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
#include <cusp/precond/block_smoothed_aggregation.h>
#include <cusp/gallery/poisson.h>
#include <cusp/csr_matrix.h>
#include <cusp/krylov/blockcg.h>
#include <cusp/block_monitor.h>

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

template <typename Orientation,
          typename IndexType, typename ValueType, typename MemorySpace>
void cusp_sa_block_cg(
  const cusp::csr_matrix<IndexType, ValueType, MemorySpace>& A,
  IndexType nrhs, IndexType max_its, ValueType tol) {

  // allocate storage for solution (x) and right hand side (b)
  cusp::array2d<ValueType, MemorySpace, Orientation> x(A.num_rows, nrhs, 0);
  cusp::array2d<ValueType, MemorySpace, Orientation> b(A.num_rows, nrhs, 1);

  // set stopping criteria
  cusp::default_block_monitor<ValueType> monitor(b, max_its, tol, 0.0, false);

  // setup preconditioner
  typedef cusp::relaxation::block_polynomial<ValueType,MemorySpace,Orientation> Smoother;
  cusp::precond::aggregation::block_smoothed_aggregation<IndexType, ValueType, MemorySpace, Smoother> M(A, nrhs);

  // solve
  {
    TEUCHOS_FUNC_TIME_MONITOR("Total Block-CG Solve Time");
    cusp::krylov::blockcg(A, x, b, monitor, M);
  }
}

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
    IndexType max_its = 100;
    CLP.setOption("max_iterations", &max_its,
                  "Maximum number of CG iterations");
    double tol = 1e-6; // has to be double
    CLP.setOption("tolerance", &tol, "Convergence tolerance");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
    CLP.parse( argc, argv );

    // Set CUDA device
    cudaSetDevice(device_id);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    // create 3D Poisson problem
    cusp::csr_matrix<IndexType, ValueType, MemorySpace> A;
    cusp::gallery::poisson27pt(A, n, n, n);

    // Create timers
    Teuchos::RCP<Teuchos::Time> tm_cg =
      Teuchos::TimeMonitor::getNewTimer("Total Block-CG Solve Time");
    Teuchos::RCP<Teuchos::Time> tm_prec =
      Teuchos::TimeMonitor::getNewTimer("CUSP Block Multilevel Solve");
    Teuchos::RCP<Teuchos::Time> tm_coarse =
      Teuchos::TimeMonitor::getNewTimer("CUSP Coarse-grid Solve");
    Teuchos::RCP<Teuchos::Time> tm_op =
      Teuchos::TimeMonitor::getNewTimer("CUSP Operator block-apply");
    Teuchos::RCP<Teuchos::Time> tm_prec_op =
      Teuchos::TimeMonitor::getNewTimer("CUSP Matrix block-apply");

    std::cout << "nrhs , num_rows , num_entries , "
              << "row_cg , row_op , row_prec , row_prec_op , row_coarse , "
              << "col_cg , col_op , col_prec , col_prec_op , col_coarse"
              << std::endl;

    for (IndexType nrhs = nrhs_begin; nrhs <= nrhs_end; nrhs += nrhs_step) {

       std::cout << nrhs << " , "
                 << A.num_rows << " , "  << A.num_entries << " , ";

      // test row-major storage
      Teuchos::TimeMonitor::zeroOutTimers();
      cusp_sa_block_cg<cusp::row_major>(A, nrhs, max_its, tol);

      std::cout << tm_cg->totalElapsedTime() << " , "
                << tm_op->totalElapsedTime() << " , "
                << tm_prec->totalElapsedTime() << " , "
                << tm_prec_op->totalElapsedTime() << " , "
                << tm_coarse->totalElapsedTime() << " , ";

      // test column-major storage
      Teuchos::TimeMonitor::zeroOutTimers();
      cusp_sa_block_cg<cusp::column_major>(A, nrhs, max_its, tol);

      std::cout << tm_cg->totalElapsedTime() << " , "
                << tm_op->totalElapsedTime() << " , "
                << tm_prec->totalElapsedTime() << " , "
                << tm_prec_op->totalElapsedTime() << " , "
                << tm_coarse->totalElapsedTime() << std::endl;

    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
