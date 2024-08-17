//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#include "Kokkos_Core.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_gauss_seidel.hpp"
#include "KokkosBlas1_nrm2.hpp"

// Parallel Gauss-Seidel Preconditioner/Smoother
//  -Uses graph coloring to find independent row sets,
//   and applies GS to each set in parallel
//  -Here, use to solve a diagonally dominant linear system directly.

// Helper to print out colors in the shape of the grid
int main() {
  using Scalar    = default_scalar;
  using Mag       = Kokkos::ArithTraits<Scalar>::mag_type;
  using Ordinal   = default_lno_t;
  using Offset    = default_size_type;
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using MemSpace  = typename ExecSpace::memory_space;
  using Device    = Kokkos::Device<ExecSpace, MemSpace>;
  using Handle =
      KokkosKernels::Experimental::KokkosKernelsHandle<Offset, Ordinal, default_scalar, ExecSpace, MemSpace, MemSpace>;
  using Matrix              = KokkosSparse::CrsMatrix<Scalar, Ordinal, Device, void, Offset>;
  using Vector              = typename Matrix::values_type;
  constexpr Ordinal numRows = 10000;
  const Scalar one          = Kokkos::ArithTraits<Scalar>::one();
  const Mag magOne          = Kokkos::ArithTraits<Mag>::one();
  // Solve tolerance
  const Mag tolerance = 1e-6 * magOne;
  Kokkos::initialize();
  {
    // Generate a square, strictly diagonally dominant, but nonsymmetric matrix
    // on which Gauss-Seidel should converge. Get approx. 20 entries per row
    // Diagonals are 2x the absolute sum of all other entries.
    Offset nnz = numRows * 20;
    Matrix A = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<Matrix>(numRows, numRows, nnz, 2, 100,
                                                                                         1.05 * one);
    std::cout << "Generated a matrix with " << numRows << " rows/cols, and " << nnz << " entries.\n";
    // Create a kernel handle, then a Gauss-Seidel handle with the default
    // algorithm
    Handle handle;
    handle.create_gs_handle(KokkosSparse::GS_DEFAULT);
    // Set up Gauss-Seidel for the graph (matrix sparsity pattern)
    KokkosSparse::Experimental::gauss_seidel_symbolic(&handle, numRows, numRows, A.graph.row_map, A.graph.entries,
                                                      false);
    // Set up Gauss-Seidel for the matrix values (numeric)
    // Another matrix with the same sparsity pattern could re-use the handle and
    // symbolic phase, and only call numeric.
    KokkosSparse::Experimental::gauss_seidel_numeric(&handle, numRows, numRows, A.graph.row_map, A.graph.entries,
                                                     A.values, false);
    // Now, preconditioner is ready to use. Set up an unknown vector
    // (uninitialized) and randomized right-hand-side vector.
    Vector x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), numRows);
    Vector b(Kokkos::view_alloc(Kokkos::WithoutInitializing, "b"), numRows);
    Vector res(Kokkos::view_alloc(Kokkos::WithoutInitializing, "res"), numRows);
    auto bHost = Kokkos::create_mirror_view(b);
    for (Ordinal i = 0; i < numRows; i++) bHost(i) = 3 * ((one * rand()) / RAND_MAX);
    Kokkos::deep_copy(b, bHost);
    // Measure initial residual norm ||Ax - b||, where x is 0
    Mag initialRes    = KokkosBlas::nrm2(b);
    Mag scaledResNorm = magOne;
    bool firstIter    = true;
    // Iterate until reaching the tolerance
    int numIters = 0;
    while (scaledResNorm > tolerance) {
      // Run one sweep of forward Gauss-Seidel (SOR with omega = 1.0)
      // If this is the first iteration, tell apply:
      //  * to zero out x (it was uninitialized)
      //  * that b has changed since the previous apply (since there was no
      //  previous apply)
      KokkosSparse::Experimental::forward_sweep_gauss_seidel_apply(
          &handle, numRows, numRows, A.graph.row_map, A.graph.entries, A.values, x, b, firstIter, firstIter, one, 1);
      firstIter = false;
      // Now, compute the new residual norm using SPMV
      Kokkos::deep_copy(res, b);
      // Compute res := Ax - res (since res is now equal to b, this is Ax - b)
      KokkosSparse::spmv("N", one, A, x, -one, res);
      // Recompute the scaled norm
      scaledResNorm = KokkosBlas::nrm2(res) / initialRes;
      numIters++;
      std::cout << "Iteration " << numIters << " scaled residual norm: " << scaledResNorm << '\n';
    }
    std::cout << "SUCCESS: converged in " << numIters << " iterations.\n";
  }
  Kokkos::finalize();
  return 0;
}
