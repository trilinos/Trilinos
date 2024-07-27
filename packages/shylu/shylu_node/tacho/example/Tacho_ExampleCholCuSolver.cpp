// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include "Tacho_CommandLineParser.hpp"
#include "Tacho_Internal.hpp"

#if defined(KOKKOS_ENABLE_CUDA)
#include "Tacho_CuSolver.hpp"
#endif

using namespace Tacho;

int main(int argc, char *argv[]) {
  CommandLineParser opts("This example program measure the performance of cuSolver on Kokkos::Cuda");

  bool verbose = true;
  bool sanitize = false;
  std::string file = "test.mtx";
  int nrhs = 1;

  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("sanitize", "Flag to sanitize input matrix (remove zeros)", &sanitize);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);

  const bool detail = false;

  typedef double value_type;

  typedef UseThisDevice<Kokkos::Cuda>::type device_type;
  typedef UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;

  Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace", detail);

  Kokkos::Timer timer;
  int r_val = 0;
#if defined(KOKKOS_ENABLE_CUDA)
  {
    ///
    /// read from crs matrix
    ///
    typedef Tacho::CrsMatrixBase<value_type, host_device_type> CrsMatrixBaseTypeHost;
    typedef Tacho::CrsMatrixBase<value_type, device_type> CrsMatrixBaseType;
    typedef Kokkos::View<value_type **, Kokkos::LayoutLeft, device_type> DenseMultiVectorType;

    /// read a spd matrix of matrix market format
    CrsMatrixBaseTypeHost h_A;
    {
      std::ifstream in;
      in.open(file);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file << std::endl;
        return -1;
      }
      Tacho::MatrixMarket<value_type>::read(file, h_A, sanitize, verbose);
    }

    ///
    /// cuSolver
    ///
    CuSolver cusolver;
    cusolver.setVerbose(verbose);

    ///
    /// reorder matrix
    ///
#if defined(TACHO_HAVE_METIS)
    typedef GraphTools_Metis graph_tools_type;
#else
    typedef GraphTools graph_tools_type;
#endif
    Graph graph(h_A.NumRows(), h_A.NumNonZeros(), h_A.RowPtr(), h_A.Cols());
    graph_tools_type G(graph);
    G.reorder(verbose);

    const auto h_perm = G.PermVector();
    const auto h_peri = G.InvPermVector();

    const auto perm = Kokkos::create_mirror_view(typename device_type::memory_space(), h_perm);
    Kokkos::deep_copy(perm, h_perm);
    const auto peri = Kokkos::create_mirror_view(typename device_type::memory_space(), h_peri);
    Kokkos::deep_copy(peri, h_peri);

    CrsMatrixBaseType A;
    A.createConfTo(h_A);
    A.copy(h_A);

    /// permute ondevice
    CrsMatrixBaseType Ap;
    {
      timer.reset();
      Ap.createConfTo(A);
      Tacho::applyPermutationToCrsMatrixLower(Ap, A, perm, peri);
      Kokkos::fence();
      const double t_permute_A = timer.seconds();

      if (verbose) {
        printf("ExampleCuSolver: Construction of permuted matrix A\n");
        printf("  Time\n");
        printf("             time for permutation of A:                       %10.6f s\n", t_permute_A);
        printf("\n");
      }
    }

    ///
    /// analyze
    ///
    { cusolver.analyze(Ap.NumRows(), Ap.RowPtr(), Ap.Cols()); }

    ///
    /// factorize
    ///
    { cusolver.factorize(Ap.Values()); }

    ///
    /// random right hand side
    ///
    DenseMultiVectorType b("b", A.NumRows(), nrhs), // rhs multivector
        x("x", A.NumRows(), nrhs),                  // solution multivector
        bb("bb", A.NumRows(), nrhs),                // temp workspace (store permuted rhs)
        xx("t", A.NumRows(), nrhs);                 // temp workspace (store permuted rhs)

    {
      Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);
      Kokkos::fill_random(b, random, value_type(1));
    }

    ///
    /// solve
    ///
    {
      timer.reset();
      const auto exec_instance = typename device_type::execution_space();
      ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, perm, bb);
      cusolver.solve(xx, bb);
      ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, xx, peri, x);
      Kokkos::fence();
      const double t_solve = timer.seconds();
      if (verbose) {
        printf("ExampleCuSolver: P b, solve, and P^{-1} x\n");
        printf("  Time\n");
        printf("             time for permute and solve:                      %10.6f s\n", t_solve);
        printf("\n");
      }
    }

    ///
    /// compute residual to check solutions
    ///
    const double res = computeRelativeResidual(A, x, b);

    std::cout << "cuSolver: residual = " << res << "\n\n";
  }
#else
  r_val = -1;
  std::cout << "CUDA is NOT configured in Trilinos" << std::endl;
#endif

  Kokkos::finalize();

  return r_val;
}
