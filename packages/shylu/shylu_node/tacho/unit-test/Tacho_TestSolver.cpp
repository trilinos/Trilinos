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
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "Tacho_CommandLineParser.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_Driver.hpp"
#include "Tacho_MatrixMarket.hpp"


/// matrix generator
template <typename CrsMatrixBaseTypeHost>
int generate2dLaplace(const int nx, CrsMatrixBaseTypeHost &A) {
  // generate Laplace matrix on 5pt stencil
  int n = nx*nx;
  int nnz = n + 4*nx*(nx-1);

  using value_type = typename CrsMatrixBaseTypeHost::value_type;
  typename CrsMatrixBaseTypeHost::size_type_array    ap("ap", n+1);
  typename CrsMatrixBaseTypeHost::ordinal_type_array aj("aj", nnz);
  typename CrsMatrixBaseTypeHost::value_type_array   ax("ax", nnz);

  nnz = 0;
  ap(0) = 0;
  for (int i=0; i<nx; i++) {
    for (int j=0; j<nx; j++) {
      int k = i*nx + j;
      if (i > 0) {
        aj(nnz) = k-nx;
        ax(nnz) = value_type(1.0);
        nnz ++;
      }
      if (j > 0) {
        aj(nnz) = k-1;
        ax(nnz) = value_type(1.0);
        nnz ++;
      }
      aj(nnz) = k;
      ax(nnz) = value_type(4.0);
      nnz++;
      if (j < nx-1) {
        aj(nnz) = k+1;
        ax(nnz) = value_type(1.0);
        nnz ++;
      }
      if (i < nx-1) {
        aj(nnz) = k+nx;
        ax(nnz) = value_type(1.0);
        nnz ++;
      }
      ap(k+1) = nnz;
    }
  }
  A.clear();
  A.setExternalMatrix(n, n, nnz, ap, aj, ax);
  return nnz;
}


// main test driver
template <typename value_type>
int driver(const std::string file, const std::string method_name, const int variant, const int nrhs,
           const bool store_transpose = false, const bool single_solve = true, const bool single_setup = true) {
  int nx = 100;
  int method = 1; // 1 - Chol, 2 - LDL, 3 - SymLU
  if (method_name == "ldl-nopiv")
    method = 0;
  else if (method_name == "chol")
    method = 1;
  else if (method_name == "ldl")
    method = 2;
  else if (method_name == "lu")
    method = 3;
  else {
    std::cout << "Error: not supported solution method\n";
    return -1;
  }

  using device_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;
  using host_device_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;

  const bool detail = false;
  const bool verbose = !single_solve;
  Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace", detail);
  std::cout << std::endl << "    --------------------- " << std::endl;
  std::cout << "     Method Name:: " << method_name << std::endl;
  std::cout << "     Solver Type:: " << variant << std::endl;
  std::cout << "          # RHSs:: " << nrhs << std::endl;
  std::cout << "          Matrix:: " << file << std::endl;
  if (!single_setup) std::cout << "     Multiple Init calls" << std::endl;
  if (!single_solve) std::cout << "     Mixed LU/Chol solve" << std::endl;
  std::cout << "    --------------------- " << std::endl << std::endl;

  int r_val = -1;
  try {
    /// crs matrix format and dense multi vector
    using CrsMatrixBaseTypeHost = Tacho::CrsMatrixBase<value_type, host_device_type>;
    using DenseMultiVectorType = Kokkos::View<value_type **, Kokkos::LayoutLeft, device_type>;

    /// read a spd matrix of matrix market format
    CrsMatrixBaseTypeHost A;
    if (file != "") {
      std::cout << " Read matrix from " << file << std::endl;
      std::ifstream in;
      in.open(file);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file << std::endl;
        return -1;
      }
      bool sanitize = false;
      Tacho::MatrixMarket<value_type>::read(file, A, sanitize, verbose);
    } else {
      std::cout << " Generate 2D Laplace matrix(nx = " << nx << ")" << std::endl;
      generate2dLaplace(nx, A);
    }

    /// create tacho solver
    Tacho::Driver<value_type, device_type> solver;

    /// set common options
    solver.setVerbose(verbose);
    solver.setSolutionMethod(method);
    solver.setLevelSetOptionAlgorithmVariant(variant);
    if (store_transpose) {
      std::cout << "  > Using explicit transpose " << std::endl;
      solver.storeExplicitTranspose(true);
    }

    /// levelset options
    ///  forcing to have a few device factor/solve tasks
    int device_factor_thres = 10;
    int device_solve_thres = 10;
    solver.setLevelSetOptionDeviceFunctionThreshold(device_factor_thres, device_solve_thres);

    auto values_on_device = Kokkos::create_mirror_view(typename device_type::memory_space(), A.Values());

    int num_setups = (single_setup ? 1 : 2); // number of symbolic calls
    for (int s = 0; s < num_setups; s++) {
      /// initialize
      r_val = solver.analyze(A.NumRows(), A.RowPtr(), A.Cols());
      if(r_val == 0) {
        r_val = solver.initialize();
      }
      int num_solves = (single_solve ? 1 : 5); // number of numeric + solve calls
      for (int step = 0; step < num_solves && r_val == 0; step++) {
        if (step > 0) {
          // perturb the first element (diagonal if Laplace), on host
          A.Values()[0]+=value_type(step);
        }
        // copy A to device
        Kokkos::deep_copy(values_on_device, A.Values());
        /// do numerical factorization
        if (single_solve) {
          r_val = solver.factorize(values_on_device);
        } else {
          if (step%2 == 1) {
            // User-specified method
            r_val = solver.factorize(values_on_device);
          } else {
            // Chol
            r_val = solver.factorize(values_on_device, 1);
          }
        }

        /// solve
        DenseMultiVectorType b("b", A.NumRows(), nrhs), // rhs multivector
          x("x", A.NumRows(), nrhs),                  // solution multivector
          t("t", A.NumRows(), nrhs);                  // temp workspace (store permuted rhs)
        if(r_val == 0) {
          const value_type zero(0.0);
          const value_type one (1.0);
          if (true) {
            Kokkos::deep_copy (b, one);
          } else {
            Kokkos::deep_copy (x, one);
            solver.computeSpMV(values_on_device, x, b);
            Kokkos::deep_copy (x, zero);
          }
          r_val = solver.solve(x, b, t);
        }
        if(r_val == 0) {
          const double tol = 0.0000000001;
          const double res = solver.computeRelativeResidual(values_on_device, x, b);
          r_val = (res <= tol ? 0 : -1);
          std::cout << " res = " << res << " tol = " << tol << "(rval = " << r_val << ")"
                    << std::endl << std::endl;
        }
      } // end of for steps
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: exception is caught: \n" << e.what() << "\n";
  } catch (...) {
    std::cerr << "Error: unknown exception is caught\n";
  }
  return r_val;
}


// tests
TEST( Solver, Chol ) {
  // Chol
  std::string file = "";
  // > single RHS
  EXPECT_EQ(driver<double>(file, "chol", 0, 1), 0);
  EXPECT_EQ(driver<double>(file, "chol", 1, 1), 0);
  EXPECT_EQ(driver<double>(file, "chol", 2, 1), 0);
  EXPECT_EQ(driver<double>(file, "chol", 3, 1), 0);
  EXPECT_EQ(driver<double>(file, "chol", 3, 1, true), 0);
  // > multiple RHSs
  EXPECT_EQ(driver<double>(file, "chol", 0, 5), 0);
  EXPECT_EQ(driver<double>(file, "chol", 1, 5), 0);
  EXPECT_EQ(driver<double>(file, "chol", 2, 5), 0);
  EXPECT_EQ(driver<double>(file, "chol", 3, 5), 0);
  EXPECT_EQ(driver<double>(file, "chol", 3, 5, true), 0);
  #if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_SYCL)
  // > sequential path
  EXPECT_EQ(driver<double>(file, "chol", -1, 1), 0);
  EXPECT_EQ(driver<double>(file, "chol", -1, 5), 0);
  #endif
}

TEST( Solver, LU ) {
  // LU
  std::string file = "";
  // > single RHS
  EXPECT_EQ(driver<double>(file, "lu", 0, 1), 0);
  EXPECT_EQ(driver<double>(file, "lu", 1, 1), 0);
  EXPECT_EQ(driver<double>(file, "lu", 2, 1), 0);
  EXPECT_EQ(driver<double>(file, "lu", 3, 1), 0);
  EXPECT_EQ(driver<double>(file, "lu", 0, 1, false, false), 0); // with mixed LU or Chol
  EXPECT_EQ(driver<double>(file, "lu", 1, 1, false, false), 0);
  EXPECT_EQ(driver<double>(file, "lu", 2, 1, false, false), 0);
  EXPECT_EQ(driver<double>(file, "lu", 3, 1, false, false), 0);
  // > multiple RHSs
  EXPECT_EQ(driver<double>(file, "lu", 0, 5), 0);
  EXPECT_EQ(driver<double>(file, "lu", 1, 5), 0);
  EXPECT_EQ(driver<double>(file, "lu", 2, 5), 0);
  EXPECT_EQ(driver<double>(file, "lu", 3, 5), 0);
  EXPECT_EQ(driver<double>(file, "lu", 0, 5, false, false), 0); // with mixed LU or Chol
  EXPECT_EQ(driver<double>(file, "lu", 1, 5, false, false), 0);
  EXPECT_EQ(driver<double>(file, "lu", 2, 5, false, false), 0);
  EXPECT_EQ(driver<double>(file, "lu", 3, 5, false, false), 0);
  EXPECT_EQ(driver<double>(file, "lu", 3, 5, false, false, false), 0); // multiple symbolic calls
  #if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_SYCL)
  // > sequential path
  EXPECT_EQ(driver<double>(file, "lu", -1, 1), 0);
  EXPECT_EQ(driver<double>(file, "lu", -1, 5), 0);
  #endif
}

TEST( Solver, LDL ) {
  // LDL
  std::string file = "";
  // > single RHS
  EXPECT_EQ(driver<double>(file, "ldl", 0, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 1, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 2, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 0, 1, false, false), 0); // with mixed LDL or Chol
  EXPECT_EQ(driver<double>(file, "ldl", 1, 1, false, false), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 2, 1, false, false), 0);
  // > multiple RHSs
  EXPECT_EQ(driver<double>(file, "ldl", 0, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 1, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 2, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 0, 5, false, false), 0); // with mixed LDL or Chol
  EXPECT_EQ(driver<double>(file, "ldl", 1, 5, false, false), 0);
  EXPECT_EQ(driver<double>(file, "ldl", 2, 5, false, false), 0);
  #if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_SYCL)
  // > sequential path
  EXPECT_EQ(driver<double>(file, "ldl", -1, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl", -1, 5), 0);
  #endif
}

TEST( Solver, NonPivLDL ) {
  // Non-piv LDL
  std::string file = "";
  // > single RHS
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 0, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 1, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 2, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 1), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 0, 1, false, false), 0); // with mixed nopiv-LDL or Chol
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 1, 1, false, false), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 2, 1, false, false), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 1, false, false), 0);
  //EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 1), 0);
  //EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 1, true), 0);
  // > multiple RHSs
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 0, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 1, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 2, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 5), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 0, 5, false, false), 0); // with mixed nopiv-LDL or Chol
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 1, 5, false, false), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 2, 5, false, false), 0);
  EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 5, false, false), 0);
  //EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 5), 0);
  //EXPECT_EQ(driver<double>(file, "ldl-nopiv", 3, 5, true), 0);
  //#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_SYCL)
  //// > sequential path
  //EXPECT_EQ(driver<double>(file, "ldl-nopiv", -1, 1), 0);
  //EXPECT_EQ(driver<double>(file, "ldl-nopiv", -1, 5), 0);
  //#endif
}

int main(int argc, char *argv[]) {
  int info = 0;
  Kokkos::initialize(argc, argv);
  {
     ::testing::InitGoogleTest(&argc, argv);
     info = RUN_ALL_TESTS();
  }
  Kokkos::finalize();
  return info;
}

