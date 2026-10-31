// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   PAMatrixImplTests.cpp
    \brief  Tests against Intrepid2::PAMatrixImpl.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_PAMatrix.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <chrono>

namespace
{
  using namespace Intrepid2;

  TEUCHOS_UNIT_TEST( PAMatrixImpl, GEMM )
{
    double relTol = 1e-13;
    double absTol = 1e-13;

    // Define column-major matrices A, B, C for each cell.
    // Compute C = A * B.
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;

    const int cellCount = 100;
    const int M = 100; // row count for A, C
    const int N = 9; // new output dimension (column count for B, C)
    const int K = 11; // contraction dimension (A's column count, B's row count)

    // pick a formula to generate the matrix A
    auto formula_A = [&] (int m, int k) -> double
    {
      return double((m+1) % 3) + double(k * M);
    };

    // pick a formula to generate the matrices B
    auto formula_B = [&] (int cellOrdinal, int k, int n) -> double
    {
      return 1.0 + double(k * cellCount) + double(n * N);
    };

    auto formula_C = [&] (int cellOrdinal, int m, int n) -> Scalar
    {
      double value = 0;
      for (int k=0; k<K; k++)
      {
        const double A_k = formula_A(m,k);
        const double B_k = formula_B(cellOrdinal,k,n);
        value += A_k * B_k;
      }
      return value;
    };


    using View1D = Kokkos::View<Scalar*,DeviceType>;
    View1D AView("PAMatrixImpl.GEMM test: AView", M * K);
    View1D BView("PAMatrixImpl.GEMM test: BView", cellCount * K * N);
    View1D CView("PAMatrixImpl.GEMM test: CView", cellCount * M * N);
    View1D CViewExpected("PAMatrixImpl.GEMM test: CViewExpected", cellCount * M * N);

    auto AViewHost = Kokkos::create_mirror(AView);
    auto BViewHost = Kokkos::create_mirror(BView);
    auto CViewExpectedHost = Kokkos::create_mirror(CViewExpected);

    std::vector<bool> BTransposedChoices {false, true};
    for (const auto storeBTransposed : BTransposedChoices) // if true, we store B as a N x K matrix
    {
      Scalar* Aptr = AViewHost.data();
      for (int k=0; k<K; k++)
      {
        for (int m=0; m<M; m++)
        {
          *Aptr++ = formula_A(m, k);
        }
      }
      Scalar* Bptr = BViewHost.data();
      if (!storeBTransposed)
      {
        for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
        {
          for (int n=0; n<N; n++)
          {
            for (int k=0; k<K; k++)
            {
              *Bptr++ = formula_B(cellOrdinal, k, n);
            }
          }
        }
      }
      else
      {
        for (int k=0; k<K; k++)
        {
          for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
          {
            for (int n=0; n<N; n++)
            {
              *Bptr++ = formula_B(cellOrdinal, k, n);
            }
          }
        }
      }
      Scalar* Cptr = CViewExpectedHost.data();
      for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
      {
        for (int n=0; n<N; n++)
        {
          for (int m=0; m<M; m++)
          {
            *Cptr++ = formula_C(cellOrdinal, m, n);
          }
        }
      }
      Kokkos::deep_copy(AView, AViewHost);
      Kokkos::deep_copy(BView, BViewHost);
      Kokkos::deep_copy(CViewExpected, CViewExpectedHost);

      const Scalar alpha = 1.0;
      const Scalar beta  = 0.0;
      const ordinal_type LDA = M;
      const char trB = storeBTransposed ? 'T' : 'N';
      Intrepid2::Impl::gemm<typename DeviceType::execution_space>('N', trB, M, N*cellCount, K, alpha, AView.data(), LDA, BView.data(), beta, CView.data());

      testFloatingEquality1(CView, CViewExpected, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( PAMatrixImpl, GEMMPerformance )
  {
    double relTol = 1e-13;
    double absTol = 1e-13;

    // Define column-major matrices A, B, C for each cell.
    // Compute C = A * B.
    // Here, we are emulating a 1D contraction such as happens in matrix-free apply
    // on an 8x8x8 hexahedral mesh, and estimating/measuring GEMM() throughput.
    // For polynomial order p, we construct a (p+1)^2 matrix representing basis values
    // at quadrature points.  The contraction is this matrix against a (p+1) x (C * (p+1)^2)
    // matrix.
    using DeviceType = DefaultTestDeviceType;
    using Scalar = double;

    const int cellCount = 512;

    std::vector<int> p_orders {1,2,4,6,8,10,12,20,50};

    for (auto p : p_orders)
    {
      const int M = p+1; // row count for A, C
      const int N = p+1; // new output dimension (column count for B, C)
      const int K = (p+1) * (p+1); // contraction dimension (A's column count, B's row count)

      // pick a formula to generate the matrix A
      auto formula_A = [&] (int m, int k) -> double
      {
        return double((m+1) % 3) + double(k * M);
      };

      // pick a formula to generate the matrices B
      auto formula_B = [&] (int cellOrdinal, int k, int n) -> double
      {
        return 1.0 + double(k * cellCount) + double(n * N);
      };

      auto formula_C = [&] (int cellOrdinal, int m, int n) -> Scalar
      {
        double value = 0;
        for (int k=0; k<K; k++)
        {
          const double A_k = formula_A(m,k);
          const double B_k = formula_B(cellOrdinal,k,n);
          value += A_k * B_k;
        }
        return value;
      };

      using View1D = Kokkos::View<Scalar*,DeviceType>;
      View1D AView("PAMatrixImpl.GEMM test: AView", M * K);
      View1D BView("PAMatrixImpl.GEMM test: BView", cellCount * K * N);
      View1D CView("PAMatrixImpl.GEMM test: CView", cellCount * M * N);
      View1D CViewExpected("PAMatrixImpl.GEMM test: CViewExpected", cellCount * M * N);

      auto AViewHost = Kokkos::create_mirror(AView);
      auto BViewHost = Kokkos::create_mirror(BView);
      auto CViewExpectedHost = Kokkos::create_mirror(CViewExpected);

      std::vector<bool> BTransposedChoices {false, true};
      for (const auto storeBTransposed : BTransposedChoices) // if true, we store B as a N x K matrix
      {
        Scalar* Aptr = AViewHost.data();
        for (int k=0; k<K; k++)
        {
          for (int m=0; m<M; m++)
          {
            *Aptr++ = formula_A(m, k);
          }
        }
        Scalar* Bptr = BViewHost.data();
        if (!storeBTransposed)
        {
          for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
          {
            for (int n=0; n<N; n++)
            {
              for (int k=0; k<K; k++)
              {
                *Bptr++ = formula_B(cellOrdinal, k, n);
              }
            }
          }
        }
        else
        {
          for (int k=0; k<K; k++)
          {
            for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
            {
              for (int n=0; n<N; n++)
              {
                *Bptr++ = formula_B(cellOrdinal, k, n);
              }
            }
          }
        }
        Scalar* Cptr = CViewExpectedHost.data();
        for (int cellOrdinal=0; cellOrdinal<cellCount; cellOrdinal++)
        {
          for (int n=0; n<N; n++)
          {
            for (int m=0; m<M; m++)
            {
              *Cptr++ = formula_C(cellOrdinal, m, n);
            }
          }
        }
        Kokkos::deep_copy(AView, AViewHost);
        Kokkos::deep_copy(BView, BViewHost);
        Kokkos::deep_copy(CViewExpected, CViewExpectedHost);

        const Scalar alpha = 1.0;
        const Scalar beta  = 0.0;
        const ordinal_type LDA = M;
        const char trB = storeBTransposed ? 'T' : 'N';
        const int numRuns = 10;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i=0; i<numRuns; i++)
        {
          Intrepid2::Impl::gemm<typename DeviceType::execution_space>('N', trB, M, N*cellCount, K, alpha, AView.data(), LDA, BView.data(), beta, CView.data());
        }
        auto end = std::chrono::high_resolution_clock::now();

        const double flops = double(M) * double(N) * cellCount * (2 * K - 1) * numRuns;
        std::chrono::duration<double> duration = end - start;
        const double throughputGFlops = flops / duration.count() / 1.0e9;

        out << "p=" << p << " (trB: " << trB << "): " << throughputGFlops << " GFlops.\n";

        testFloatingEquality1(CView, CViewExpected, relTol, absTol, out, success);
      }
    }
  }

} // anonymous namespace
