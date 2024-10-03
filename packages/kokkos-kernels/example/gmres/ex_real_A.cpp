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

#include <math.h>
#include "KokkosSparse_IOUtils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas3_trsm.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_gmres.hpp>

int main(int argc, char* argv[]) {
  using ST   = double;
  using OT   = int;
  using EXSP = Kokkos::DefaultExecutionSpace;
  using MESP = typename EXSP::memory_space;
  using CRS  = KokkosSparse::CrsMatrix<ST, OT, EXSP, void, OT>;

  using ViewVectorType = Kokkos::View<ST*, Kokkos::LayoutLeft, EXSP>;
  using KernelHandle   = KokkosKernels::Experimental::KokkosKernelsHandle<OT, OT, ST, EXSP, MESP, MESP>;

  std::string filename("bcsstk09.mtx");  // example matrix
  std::string ortho("CGS2");             // orthog type
  int m          = 50;                   // Max subspace size before restarting.
  double convTol = 1e-10;                // Relative residual convergence tolerance.
  int cycLim     = 50;                   // Maximum number of times to restart the solver.
  bool rand_rhs  = false;                // Generate random right-hand side.

  for (int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (token == std::string("--filename")) filename = argv[++i];
    if (token == std::string("--max-subsp")) m = std::atoi(argv[++i]);
    if (token == std::string("--max-restarts")) cycLim = std::atoi(argv[++i]);
    if (token == std::string("--tol")) convTol = std::stod(argv[++i]);
    if (token == std::string("--ortho")) ortho = argv[++i];
    if (token == std::string("--rand_rhs")) rand_rhs = true;
    if (token == std::string("--help") || token == std::string("-h")) {
      std::cout << "Kokkos GMRES solver options:" << std::endl
                << "--filename    :  The name of a matrix market (.mtx) file for "
                   "matrix A (Default bcsstk09.mtx)."
                << std::endl
                << "--max-subsp   :  The maximum size of the Kyrlov subspace before "
                   "restarting (Default 50)."
                << std::endl
                << "--max-restarts:  Maximum number of GMRES restarts (Default 50)." << std::endl
                << "--tol         :  Convergence tolerance.  (Default 1e-10)." << std::endl
                << "--ortho       :  Type of orthogonalization. Use 'CGS2' or 'MGS'. "
                   "(Default 'CGS2')"
                << std::endl
                << "--rand_rhs    :  Generate a random right-hand side b.  (Else, "
                   "default uses b = vector of ones.)"
                << std::endl
                << "--help  -h    :  Display this help message." << std::endl
                << "Example Call  :  ./Gmres.exe --filename Laplace3D100.mtx --tol "
                   "1e-5 --max-subsp 100 "
                << std::endl
                << std::endl;
      return 0;
    }
  }
  std::cout << "File to process is: " << filename << std::endl;
  std::cout << "Convergence tolerance is: " << convTol << std::endl;

  // Initialize Kokkos AFTER parsing parameters:
  Kokkos::initialize();
  {
    // Read in a matrix Market file and use it to test the Kokkos Operator.
    CRS A = KokkosSparse::Impl::read_kokkos_crst_matrix<CRS>(filename.c_str());

    int n = A.numRows();
    ViewVectorType X("X", n);    // Solution and initial guess
    ViewVectorType Wj("Wj", n);  // For checking residuals at end.
    ViewVectorType B(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"),
                     n);  // right-hand side vec

    // Make kernel handles and set options
    KernelHandle kh;
    kh.create_gmres_handle(m, convTol, cycLim);
    auto gmres_handle = kh.get_gmres_handle();
    // Get full gmres handle type using decltype. Deferencing a pointer gives a
    // reference, so we need to strip that too.
    using GMRESHandle = typename std::remove_reference<decltype(*gmres_handle)>::type;
    gmres_handle->set_ortho(ortho == "CGS2" ? GMRESHandle::Ortho::CGS2 : GMRESHandle::Ortho::MGS);

    if (rand_rhs) {
      // Make rhs random.
      int rand_seed = 123;
      Kokkos::Random_XorShift64_Pool<> pool(rand_seed);
      Kokkos::fill_random(B, pool, -1, 1);
    } else {
      // Make rhs ones so that results are repeatable:
      Kokkos::deep_copy(B, 1.0);
    }

    // Run GMRS solve:
    KokkosSparse::Experimental::gmres(&kh, A, B, X);

    const auto numIters  = gmres_handle->get_num_iters();
    const auto convFlag  = gmres_handle->get_conv_flag_val();
    const auto endRelRes = gmres_handle->get_end_rel_res();

    // Double check residuals at end of solve:
    ST nrmB = KokkosBlas::nrm2(B);
    KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
    KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
    ST endRes = KokkosBlas::nrm2(B) / nrmB;
    std::cout << "=========================================" << std::endl;
    std::cout << "Verify from main: Ending residual is " << endRes << std::endl;
    std::cout << "Number of iterations is: " << numIters << std::endl;
    std::cout << "Diff of residual from main - residual from solver: " << endRelRes - endRes << std::endl;
    std::cout << "Convergence flag is : " << convFlag << std::endl;
  }
  Kokkos::finalize();
}
