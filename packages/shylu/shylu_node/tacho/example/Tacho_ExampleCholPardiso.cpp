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

#if defined(TACHO_HAVE_MKL)
#include "mkl.h"
#include "Tacho_Pardiso.hpp"
#endif

using namespace Tacho;

int main(int argc, char *argv[]) {
  CommandLineParser opts("This example program measure the performance of Pardiso on Kokkos::OpenMP");

  int nthreads = 1;
  bool verbose = true;
  bool sanitize = false;
  std::string file_input = "test.mtx";
  int nrhs = 1;

  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<bool>("sanitize", "Flag to sanitize input matrix (remove zeros)", &sanitize);
  opts.set_option<std::string>("file", "Input file (MatrixMarket SPD matrix)", &file_input);
  opts.set_option<int>("nrhs", "Number of RHS vectors", &nrhs);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  const bool skip_factorize = false, skip_solve = false;

  Kokkos::initialize(argc, argv);

  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;
  printExecSpaceConfiguration<typename host_device_type::execution_space>("HostDevice", false);

  int r_val = 0;
#if defined(__INTEL_MKL__)
  {
    typedef double value_type;
    typedef CrsMatrixBase<value_type, host_device_type> CrsMatrixBaseType;
    typedef Kokkos::View<value_type **, Kokkos::LayoutLeft, host_device_type> DenseMatrixBaseType;

    // mkl nthreads setting
    mkl_set_dynamic(0);
    mkl_set_num_threads(nthreads);

    Kokkos::Timer timer;
    double t = 0.0;
    Pardiso pardiso;

    constexpr int AlgoChol = 2;
    std::cout << "PardisoChol:: init" << std::endl;
    {
      timer.reset();
      r_val = pardiso.init<value_type, AlgoChol>();
      t = timer.seconds();

      if (r_val) {
        std::cout << "PardisoChol:: Pardiso init error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      }
    }
    std::cout << "PardisoChol:: init ::time = " << t << std::endl;

    std::cout << "PardisoChol:: import input file = " << file_input << std::endl;
    CrsMatrixBaseType A, Asym;
    timer.reset();
    {
      {
        std::ifstream in;
        in.open(file_input);
        if (!in.good()) {
          std::cout << "Failed in open the file: " << file_input << std::endl;
          return -1;
        }
      }
      MatrixMarket<value_type>::read(file_input, A, sanitize, verbose);

      // somehow pardiso does not like symmetric full matrix (store only half)
      Asym.createConfTo(A);
      {
        size_type nnz = 0;
        for (ordinal_type i = 0; i < A.NumRows(); ++i) {
          Asym.RowPtrBegin(i) = nnz;
          for (ordinal_type idx = A.RowPtrBegin(i); idx < A.RowPtrEnd(i); ++idx) {
            if (i <= A.Col(idx)) {
              Asym.Col(nnz) = A.Col(idx);
              Asym.Value(nnz) = A.Value(idx);
              ++nnz;
            }
          }
          Asym.RowPtrEnd(i) = nnz;
        }
      }
    }
    t = timer.seconds();

    // 32bit vs 64bit integers; A uses size_t for size array
    Kokkos::View<ordinal_type *, host_device_type> rowptr("rowptr", Asym.NumRows() + 1);
    {
      for (ordinal_type i = 0; i <= Asym.NumRows(); ++i)
        rowptr(i) = Asym.RowPtrBegin(i);
    }
    std::cout << "PardisoChol:: import input file::time = " << t << std::endl;

    DenseMatrixBaseType B("B", Asym.NumRows(), nrhs), X("X", Asym.NumRows(), nrhs), P("P", Asym.NumRows(), 1);

    {
      const auto m = Asym.NumRows();
      Random<value_type> random;
      for (ordinal_type rhs = 0; rhs < nrhs; ++rhs)
        for (ordinal_type i = 0; i < m; ++i)
          B(i, rhs) = random.value();
      Kokkos::deep_copy(X, B);
    }

    pardiso.setProblem(Asym.NumRows(), (double *)Asym.Values().data(),
                       (int *)rowptr.data(), // (int*)Asym.RowPtr().data(),
                       (int *)Asym.Cols().data(), (int *)P.data(), nrhs, (double *)B.data(), (double *)X.data());

    std::cout << "PardisoChol:: analyze matrix" << std::endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::Analyze);
      t = timer.seconds();

      if (r_val) {
        std::cout << "PardisoChol:: Pardiso analyze error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::Analyze) << std::endl;
      }
    }
    std::cout << "PardisoChol:: analyze matrix::time = " << t << std::endl;

    if (!skip_factorize) {
      std::cout << "PardisoChol:: factorize matrix" << std::endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Factorize);
        t = timer.seconds();

        if (r_val) {
          std::cout << "PardisoChol:: Pardiso factorize error = " << r_val << std::endl;
          pardiso.showErrorCode(std::cout) << std::endl;
        } else {
          pardiso.showStat(std::cout, Pardiso::Factorize) << std::endl;
        }
      }
      std::cout << "PardisoChol:: factorize matrix::time = " << t << std::endl;
    }

    if (!skip_factorize && !skip_solve) {
      std::cout << "PardisoChol:: solve matrix" << std::endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Solve);
        t = timer.seconds();

        if (r_val) {
          std::cout << "PardisoChol:: Pardiso solve error = " << r_val << std::endl;
          pardiso.showErrorCode(std::cout) << std::endl;
        } else {
          pardiso.showStat(std::cout, Pardiso::Solve) << std::endl;
        }
      }
      std::cout << "PardisoChol:: solve matrix::time = " << t << std::endl;
    }

    {
      const double res = computeRelativeResidual(A, X, B);
      std::cout << "PardisoChol:: residual = " << res << std::endl;
    }

    std::cout << "PardisoChol:: release all" << std::endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::ReleaseAll);
      t = timer.seconds();

      if (r_val) {
        std::cout << "PardisoChol:: release error = " << r_val << std::endl;
        pardiso.showErrorCode(std::cout) << std::endl;
      } else {
        pardiso.showStat(std::cout, Pardiso::ReleaseAll) << std::endl;
      }
    }
    std::cout << "PardisoChol:: release all::time = " << t << std::endl;
  }
#else
  r_val = -1;
  std::cout << "MKL is NOT configured in Trilinos" << std::endl;
#endif

  Kokkos::finalize();

  return r_val;
}
