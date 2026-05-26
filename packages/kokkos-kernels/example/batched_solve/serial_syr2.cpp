// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Syr2.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched syr2
/// Perform a symmetric rank-2 update of matrix A by vectors x and y, where
///   A: [[1, -3, -2,  0],
///       [0,  1, -3, -2],
///       [0,  0,  1, -3],
///       [0,  0,  0,  1]]
///   x: [1, 2, 3, 4]
///   y: [4, 3, 2, 1]
///   alpha: 1.5
///
/// After, A = A + alpha * (x * y^T) + alpha * (y * x^T), it will give
///   A: [[13.,  13.5,  19.,  25.5],
///       [ 0.,  19.,  16.5,  19. ],
///       [ 0.,   0.,  19.,  13.5 ],
///       [ 0.,   0.,   0.,  13.  ]]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, n = 4;

    // Matrix A
    View3DType A("A", Nb, n, n), Ref("Ref", Nb, n, n);

    // Vectors x and y
    View2DType x("x", Nb, n), y("y", Nb, n);

    // Initialize A, x, and y
    auto h_A   = Kokkos::create_mirror_view(A);
    auto h_x   = Kokkos::create_mirror_view(x);
    auto h_y   = Kokkos::create_mirror_view(y);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vectors x and y
      for (int i = 0; i < n; i++) {
        h_x(ib, i) = i + 1;
        h_y(ib, i) = n - i;
      }

      // Fill the matrix A (upper triangular)
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
      h_Ref(ib, 0, 0) = 13.0;
      h_Ref(ib, 0, 1) = 13.5;
      h_Ref(ib, 0, 2) = 19.0;
      h_Ref(ib, 0, 3) = 25.5;

      h_Ref(ib, 1, 0) = 0.0;
      h_Ref(ib, 1, 1) = 19.0;
      h_Ref(ib, 1, 2) = 16.5;
      h_Ref(ib, 1, 3) = 19.0;

      h_Ref(ib, 2, 0) = 0.0;
      h_Ref(ib, 2, 1) = 0.0;
      h_Ref(ib, 2, 2) = 19.0;
      h_Ref(ib, 2, 3) = 13.5;

      h_Ref(ib, 3, 0) = 0.0;
      h_Ref(ib, 3, 1) = 0.0;
      h_Ref(ib, 3, 2) = 0.0;
      h_Ref(ib, 3, 3) = 13.0;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    // Compute A = A + alpha * (x * y^T) + alpha * (y * x^T)
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "syr2", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);
          auto sub_y = Kokkos::subview(y, ib, Kokkos::ALL);

          // A = A + alpha * (x * y^T) + alpha * (y * x^T)
          KokkosBatched::SerialSyr2<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose>::invoke(alpha, sub_x,
                                                                                                         sub_y, sub_A);
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
      std::cout << "syr2 works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
