// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Syrk.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched syrk
/// ArgTrans == NoTranspose/ConjNoTranspose
/// Perform a symmetric rank-k update of matrix C by matrix A, where
///   A: [[1,  -3, -2,  0],
///       [-3,  3, -1, -2],
///       [-2, -1,  9,  5],
///       [0,  -2,  5, 27]]
///   C == A
///   alpha: 1.5, beta: 1.2
///
/// After, C = alpha * U * U**T + beta * C, it will give
///   C: [[22.2, -18.6, -27.9,  -6.0],
///       [-3,    38.1, -25.2, -99.9],
///       [-2,      -1, 177.3, 279.0],
///       [0,       -2,   5, 1169.4]]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Matrix A, C
    View3DType A("A", Nb, n, n), C("C", Nb, n, n), Ref("Ref", Nb, n, n);

    // Initialize A, C and Ref
    auto h_A   = Kokkos::create_mirror_view(A);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      h_A(ib, 0, 0) = 1;
      h_A(ib, 0, 1) = -3;
      h_A(ib, 0, 2) = -2;
      h_A(ib, 0, 3) = 0;
      h_A(ib, 1, 0) = -3;
      h_A(ib, 1, 1) = 3;
      h_A(ib, 1, 2) = -1;
      h_A(ib, 1, 3) = -2;
      h_A(ib, 2, 0) = -2;
      h_A(ib, 2, 1) = -1;
      h_A(ib, 2, 2) = 9;
      h_A(ib, 2, 3) = 5;
      h_A(ib, 3, 0) = 0;
      h_A(ib, 3, 1) = -2;
      h_A(ib, 3, 2) = 5;
      h_A(ib, 3, 3) = 27;

      h_Ref(ib, 0, 0) = 22.2;
      h_Ref(ib, 0, 1) = -18.6;
      h_Ref(ib, 0, 2) = -27.9;
      h_Ref(ib, 0, 3) = -6.0;
      h_Ref(ib, 1, 0) = -3;
      h_Ref(ib, 1, 1) = 38.1;
      h_Ref(ib, 1, 2) = -25.2;
      h_Ref(ib, 1, 3) = -99.9;
      h_Ref(ib, 2, 0) = -2;
      h_Ref(ib, 2, 1) = -1;
      h_Ref(ib, 2, 2) = 177.3;
      h_Ref(ib, 2, 3) = 279.0;
      h_Ref(ib, 3, 0) = 0;
      h_Ref(ib, 3, 1) = -2;
      h_Ref(ib, 3, 2) = 5;
      h_Ref(ib, 3, 3) = 1169.4;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(C, A);

    // Compute C = alpha * A * A**T + beta * C
    const double alpha = 1.5, beta = 1.2;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "syrk", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_C = Kokkos::subview(C, ib, Kokkos::ALL, Kokkos::ALL);

          // C = alpha * A * A**T + beta * C
          KokkosBatched::SerialSyrk<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose>::invoke(alpha, sub_A,
                                                                                                           beta, sub_C);
        });

    // Confirm that the results are correct
    auto h_C     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (Kokkos::abs(h_C(ib, i, j) - h_Ref(ib, i, j)) > eps) correct = false;
        }
      }
    }

    if (correct) {
      std::cout << "syrk works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
