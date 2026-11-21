// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Syr2k.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched syr2k
/// Perform a symmetric rank-2 update of matrix A by vectors x and y, where
///   alpha = 1.5, beta = 1.2
///   4x4 matrix (upper)
///   A: [[1,  -3, -2,  0],
///       [3,  3, -1, -2],
///       [2,  1,  9, 5],
///       [0,  2, -5, 27]]
///   B: [[2, 1, 0, -1],
///       [-1, 4, 2, 0],
///       [3, -2, 1, 5],
///       [1, 0, -3, 2]]
///   C = A
///   ArgTrans == NoTranspose
///   ArgUplo == Upper
///   C: = alpha*A*B^T + alpha*B*A^T + beta * C
///   Ref: [[ -1.8 -12.6   8.1 -27. ]
///         [  3.   24.6  16.8  -2.4]
///         [  2.    1.  124.8 172.5]
///         [  0.    2.   -5.  239.4]]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Matrix A, B, and C
    View3DType A("A", Nb, n, n), B("B", Nb, n, n), C("C", Nb, n, n), Ref("Ref", Nb, n, n);

    // Initialize A, B and C
    auto h_A   = Kokkos::create_mirror_view(A);
    auto h_B   = Kokkos::create_mirror_view(B);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      h_A(ib, 0, 0) = 1;
      h_A(ib, 0, 1) = -3;
      h_A(ib, 0, 2) = -2;
      h_A(ib, 0, 3) = 0;
      h_A(ib, 1, 0) = 3;
      h_A(ib, 1, 1) = 3;
      h_A(ib, 1, 2) = -1;
      h_A(ib, 1, 3) = -2;
      h_A(ib, 2, 0) = 2;
      h_A(ib, 2, 1) = 1;
      h_A(ib, 2, 2) = 9;
      h_A(ib, 2, 3) = 5;
      h_A(ib, 3, 0) = 0;
      h_A(ib, 3, 1) = 2;
      h_A(ib, 3, 2) = -5;
      h_A(ib, 3, 3) = 27;

      h_B(ib, 0, 0) = 2;
      h_B(ib, 0, 1) = 1;
      h_B(ib, 0, 2) = 0;
      h_B(ib, 0, 3) = -1;
      h_B(ib, 1, 0) = -1;
      h_B(ib, 1, 1) = 4;
      h_B(ib, 1, 2) = 2;
      h_B(ib, 1, 3) = 0;
      h_B(ib, 2, 0) = 3;
      h_B(ib, 2, 1) = -2;
      h_B(ib, 2, 2) = 1;
      h_B(ib, 2, 3) = 5;
      h_B(ib, 3, 0) = 1;
      h_B(ib, 3, 1) = 0;
      h_B(ib, 3, 2) = -3;
      h_B(ib, 3, 3) = 2;

      // Fill the reference matrix
      h_Ref(ib, 0, 0) = -1.8;
      h_Ref(ib, 0, 1) = -12.6;
      h_Ref(ib, 0, 2) = 8.1;
      h_Ref(ib, 0, 3) = -27.0;
      h_Ref(ib, 1, 0) = 3;
      h_Ref(ib, 1, 1) = 24.6;
      h_Ref(ib, 1, 2) = 16.8;
      h_Ref(ib, 1, 3) = -2.4;
      h_Ref(ib, 2, 0) = 2;
      h_Ref(ib, 2, 1) = 1;
      h_Ref(ib, 2, 2) = 124.8;
      h_Ref(ib, 2, 3) = 172.5;
      h_Ref(ib, 3, 0) = 0;
      h_Ref(ib, 3, 1) = 2;
      h_Ref(ib, 3, 2) = -5;
      h_Ref(ib, 3, 3) = 239.4;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(B, h_B);
    Kokkos::deep_copy(C, A);  // Initialize C with A

    // Compute C = alpha * A * B^T + alpha * B * A^T + beta * C
    const double alpha = 1.5, beta = 1.2;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "syr2k", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_B = Kokkos::subview(B, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_C = Kokkos::subview(C, ib, Kokkos::ALL, Kokkos::ALL);

          // C = alpha * A * B^T + alpha * B * A^T + beta * C
          KokkosBatched::SerialSyr2k<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose>::invoke(
              alpha, sub_A, sub_B, beta, sub_C);
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
      std::cout << "syr2k works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
