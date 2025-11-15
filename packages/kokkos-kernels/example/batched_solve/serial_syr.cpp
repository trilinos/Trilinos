// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Syr.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched syr
/// Perform a symmetric rank-1 update of matrix A by vector x, where
///   A: [[1, -3, -2,  0],
///       [0,  1, -3, -2],
///       [0,  0,  1, -3],
///       [0,  0,  0,  1]]
///   x: [1, 2, 3, 4]
///   alpha: 1.5
///
/// After, A = A + alpha * (x * x^T), it will give
///   A: [[ 2.5,  0.,   2.5,  6.],
///       [ 0.,   7.,   6.,  10.],
///       [ 0.,   0.,  14.5, 15.],
///       [ 0.,   0.,   0.,  25.]]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Matrix A
    View3DType A("A", Nb, n, n), Ref("Ref", Nb, n, n);

    // Vector x
    View2DType x("x", Nb, n);

    // Initialize A and x
    Kokkos::deep_copy(x, 1.0);
    auto h_A   = Kokkos::create_mirror_view(A);
    auto h_x   = Kokkos::create_mirror_view(x);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vector x
      for (int i = 0; i < n; i++) {
        h_x(ib, i) = i + 1;
      }

      // Fill the matrix A
      h_A(ib, 0, 0) = 1.0;
      h_A(ib, 0, 1) = -3.0;
      h_A(ib, 0, 2) = -2.0;
      h_A(ib, 0, 3) = 0.0;

      h_A(ib, 1, 0) = 0.0;
      h_A(ib, 1, 1) = 1.0;
      h_A(ib, 1, 2) = -3.0;
      h_A(ib, 1, 3) = -2.0;

      h_A(ib, 2, 0) = 0.0;
      h_A(ib, 2, 1) = 0.0;
      h_A(ib, 2, 2) = 1.0;
      h_A(ib, 2, 3) = -3.0;

      h_A(ib, 3, 0) = 0.0;
      h_A(ib, 3, 1) = 0.0;
      h_A(ib, 3, 2) = 0.0;
      h_A(ib, 3, 3) = 1.0;

      // Fill the reference matrix
      h_Ref(ib, 0, 0) = 2.5;
      h_Ref(ib, 0, 1) = 0.0;
      h_Ref(ib, 0, 2) = 2.5;
      h_Ref(ib, 0, 3) = 6.0;

      h_Ref(ib, 1, 0) = 0.0;
      h_Ref(ib, 1, 1) = 7.0;
      h_Ref(ib, 1, 2) = 6.0;
      h_Ref(ib, 1, 3) = 10.0;

      h_Ref(ib, 2, 0) = 0.0;
      h_Ref(ib, 2, 1) = 0.0;
      h_Ref(ib, 2, 2) = 14.5;
      h_Ref(ib, 2, 3) = 15.0;

      h_Ref(ib, 3, 0) = 0.0;
      h_Ref(ib, 3, 1) = 0.0;
      h_Ref(ib, 3, 2) = 0.0;
      h_Ref(ib, 3, 3) = 25.0;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(x, h_x);

    // Compute A = A + alpha * (x * x^T)
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "syr", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);

          // A = A + alpha * (x * x^T)
          KokkosBatched::SerialSyr<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose>::invoke(alpha, sub_x,
                                                                                                        sub_A);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_A, A);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (Kokkos::abs(h_A(ib, i, j) - h_Ref(ib, i, j)) > eps) correct = false;
        }
      }
    }

    if (correct) {
      std::cout << "syr works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
