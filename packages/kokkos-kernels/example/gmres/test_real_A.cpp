/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Jennifer Loe (jloe@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <math.h>
#include "KokkosKernels_IOUtils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas3_trsm.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_IOUtils.hpp>

#include "gmres.hpp"

int main(int /*argc*/, char** /*argv[]*/) {
  typedef double ST;
  typedef int OT;
  typedef Kokkos::DefaultExecutionSpace EXSP;

  //  std::cout << "ST: " << typeid(ST).name() << std::endl;
  //  std::cout << "EXSP: " << typeid(EXSP).name() << std::endl;

  using sp_matrix_type = KokkosSparse::CrsMatrix<ST, OT, EXSP>;
  using ViewVectorType = Kokkos::View<ST*, Kokkos::LayoutLeft, EXSP>;
  typedef sp_matrix_type::non_const_ordinal_type ncOT;
  typedef sp_matrix_type::const_ordinal_type cOT;
  typedef sp_matrix_type::non_const_size_type ncST;

  GmresOpts<ST> solverOpts;
  solverOpts.ortho      = "CGS2";  // orthog type
  solverOpts.m          = 15;      // Max subspace size before restarting.
  solverOpts.tol        = 1e-10;   // Relative residual convergence tolerance.
  solverOpts.maxRestart = 50;
  bool pass1            = false;
  bool pass2            = false;

  std::cout << "Convergence tolerance is: " << solverOpts.tol << std::endl;

  // Initialize Kokkos AFTER parsing parameters:
  Kokkos::initialize();
  {
    // Create a diagonally dominant sparse matrix to test:
    ncST nnz;
    cOT n             = 5000;
    cOT numRows       = n;
    cOT numCols       = n;
    cOT diagDominance = 1;
    nnz               = 10 * numRows;
    sp_matrix_type A =
        KokkosKernels::Impl::kk_generate_diagonally_dominant_sparse_matrix<
            sp_matrix_type>(numRows, numCols, nnz, 0, ncOT(0.01 * numRows),
                            diagDominance);

    // Set initial vectors:
    ViewVectorType X("X", n);    // Solution and initial guess
    ViewVectorType Wj("Wj", n);  // For checking residuals at end.
    ViewVectorType B(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"),
                     n);  // right-hand side vec

    // Make rhs ones so that results are repeatable:
    Kokkos::deep_copy(B, 1.0);

    std::cout << "Testing GMRES with CGS2 ortho:" << std::endl;
    GmresStats solveStats =
        gmres<ST, Kokkos::LayoutLeft, EXSP>(A, B, X, solverOpts);

    // Double check residuals at end of solve:
    double nrmB = static_cast<double>(KokkosBlas::nrm2(B));
    KokkosSparse::spmv("N", ST(1.0), A, X, ST(0.0), Wj);  // wj = Ax
    KokkosBlas::axpy(ST(-1.0), Wj, B);                    // b = b-Ax.
    double endRes = KokkosBlas::nrm2(B) / nrmB;
    std::cout << "=======================================" << std::endl;
    std::cout << "Verify from main: Ending residual is " << endRes << std::endl;
    std::cout << "Number of iterations is: " << solveStats.numIters
              << std::endl;
    std::cout << "Diff of residual from main - residual from solver: "
              << solveStats.endRelRes - endRes << std::endl;
    std::cout << "Convergence flag is : " << solveStats.convFlag() << std::endl;

    if (solveStats.numIters < 40 && solveStats.numIters > 20 &&
        endRes < static_cast<double>(solverOpts.tol)) {
      std::cout << "Test CGS2 Passed!" << std::endl;
      pass1 = true;
    } else {
      std::cout << "Solver did not converge within the expected number of "
                   "iterations. "
                << std::endl
                << "CGS2 Test Failed." << std::endl;
    }
    std::cout << "=======================================" << std::endl
              << std::endl
              << std::endl;

    solverOpts.ortho = "MGS";
    Kokkos::deep_copy(X, 0.0);
    Kokkos::deep_copy(B, 1.0);

    std::cout << "Testing GMRES with MGS ortho:" << std::endl;
    solveStats = gmres<ST, Kokkos::LayoutLeft, EXSP>(A, B, X, solverOpts);

    // Double check residuals at end of solve:
    nrmB = static_cast<double>(KokkosBlas::nrm2(B));
    KokkosSparse::spmv("N", ST(1.0), A, X, ST(0.0), Wj);  // wj = Ax
    KokkosBlas::axpy(ST(-1.0), Wj, B);                    // b = b-Ax.
    endRes = KokkosBlas::nrm2(B) / nrmB;
    std::cout << "=======================================" << std::endl;
    std::cout << "Verify from main: Ending residual is " << endRes << std::endl;
    std::cout << "Number of iterations is: " << solveStats.numIters
              << std::endl;
    std::cout << "Diff of residual from main - residual from solver: "
              << solveStats.endRelRes - endRes << std::endl;
    std::cout << "Convergence flag is : " << solveStats.convFlag() << std::endl;

    if (solveStats.numIters < 40 && solveStats.numIters > 20 &&
        endRes < static_cast<double>(solverOpts.tol)) {
      std::cout << "Test MGS Passed!" << std::endl;
      if (pass1) {
        pass2 = true;
      };
    } else {
      std::cout << "Solver did not converge within the expected number of "
                   "iterations. "
                << std::endl
                << "MGS Test Failed." << std::endl;
    }
    std::cout << "=======================================" << std::endl
              << std::endl
              << std::endl;
  }
  Kokkos::finalize();

  if (pass2) {
    std::cout << "Both tests have passed!!" << std::endl;
  } else {
    std::cout << "One or more tests has failed." << std::endl;
  }

  return (pass2 ? EXIT_SUCCESS : EXIT_FAILURE);
}
