// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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

template <typename IndexType, typename ValueType, typename MemorySpace,
          typename Orientation>
void cusp_sa_block_cg(IndexType nx, IndexType ny, IndexType nz, IndexType nrhs,
                      IndexType max_its, ValueType tol, bool verbose) {
  // create an empty sparse matrix structure
  cusp::csr_matrix<IndexType, ValueType, MemorySpace> A;

  // create 2D Poisson problem
  cusp::gallery::poisson27pt(A, nx, ny, nz);
  std::cout << "N =" << A.num_rows<< std::endl;
  std::cout << "nnz of A = " << A.num_entries << std::endl;

  // solve with multiple RHS
  std::cout << "\nSolving Ax = b with multiple RHS..." << std::endl;

  // allocate storage for solution (x) and right hand side (b)
  cusp::array2d<ValueType, MemorySpace, Orientation> x(A.num_rows, nrhs, 0);
  cusp::array2d<ValueType, MemorySpace, Orientation> b(A.num_rows, nrhs, 1);

  std::cout << "numRHS = " << nrhs << std::endl;

  // set stopping criteria
  cusp::default_block_monitor<ValueType> monitor(b, max_its, tol, 0);

  // setup preconditioner
  typedef cusp::relaxation::block_polynomial<ValueType,MemorySpace,Orientation> Smoother;
  cusp::precond::aggregation::block_smoothed_aggregation<IndexType, ValueType, MemorySpace, Smoother> M(A, nrhs);
  if (verbose) {
    // print hierarchy information
    std::cout << "\nPreconditioner statistics" << std::endl;
    M.print();
  }

  // solve
  {
    TEUCHOS_FUNC_TIME_MONITOR("Total Block-CG Solve Time");
    cusp::krylov::blockcg(A, x, b, monitor, M);
  }
}

// Orientation types
enum Orient { ROW, COL };
const int num_orient = 2;
const Orient orient_values[] = { ROW, COL };
const char *orient_names[] = { "row", "col" };

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
    CLP.setDocString(
      "This example runs AMG-preconditioned block-CG with CUSP.\n");

    IndexType nx = 32;
    CLP.setOption("nx", &nx, "Number of mesh points in the x-direction");
    IndexType ny = 32;
    CLP.setOption("ny", &ny, "Number of mesh points in the y-direction");
    IndexType nz = 32;
    CLP.setOption("nz", &nz, "Number of mesh points in the z-direction");
    IndexType nrhs = 32;
    CLP.setOption("nrhs", &nrhs, "Number of right-hand-sides");
    Orient orient = ROW;
    CLP.setOption("orient", &orient, num_orient, orient_values, orient_names,
                  "Orientation of block RHS");
    IndexType max_its = 100;
    CLP.setOption("max_iterations", &max_its,
                  "Maximum number of CG iterations");
    double tol = 1e-6; // has to be double
    CLP.setOption("tolerance", &tol, "Convergence tolerance");
    int device_id = 0;
    CLP.setOption("device", &device_id, "CUDA device ID");
    CLP.setOption("verbose", "quiet", &verbose, "Verbose output");
    CLP.parse( argc, argv );

    // Set CUDA device
    cudaSetDevice(device_id);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    if (orient == ROW)
      cusp_sa_block_cg<IndexType,ValueType,MemorySpace,cusp::row_major>(
        nx, ny, nz, nrhs, max_its, tol, verbose);
    else if (orient == COL)
      cusp_sa_block_cg<IndexType,ValueType,MemorySpace,cusp::column_major>(
        nx, ny, nz, nrhs, max_its, tol, verbose);

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if (success)
    return 0;
  return -1;
}
