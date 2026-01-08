// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Ger.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched ger
/// Perform a rank-1 update of matrix A by vectors x and y, where
///   A: [[1, -3, -2,  0],
///       [-1, 1, -3, -2],
///       [2, -1,  1, -3]]
///   x: [1, 2, 3]
///   y: [0, 1, 2, 3]
///   alpha: 1.5
///
/// After, A = A + alpha * (x * y^T), it will give
///   A: [[ 1.,  -1.5,  1.,   4.5],
///       [-1.,   4.,   3.,   7. ],
///       [ 2.,   3.5, 10.,  10.5],]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    using View3DType = Kokkos::View<double***, ExecutionSpace>;
    const int Nb = 10, m = 3, n = 4;

    // Matrix A
    View3DType A("A", Nb, m, n), Ref("Ref", Nb, m, n);

    // vectros x and y
    View2DType x("x", Nb, m), y("y", Nb, n);

    // Initialize A and x
    Kokkos::deep_copy(x, 1.0);
    auto h_A   = Kokkos::create_mirror_view(A);
    auto h_x   = Kokkos::create_mirror_view(x);
    auto h_y   = Kokkos::create_mirror_view(y);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vectors x and y
      for (int i = 0; i < m; i++) {
        h_x(ib, i) = i + 1;
      }
      for (int i = 0; i < n; i++) {
        h_y(ib, i) = i;
      }

      // Fill the matrix A
      h_A(ib, 0, 0) = 1.0;
      h_A(ib, 0, 1) = -3.0;
      h_A(ib, 0, 2) = -2.0;
      h_A(ib, 0, 3) = 0.0;
      h_A(ib, 1, 0) = -1.0;
      h_A(ib, 1, 1) = 1.0;
      h_A(ib, 1, 2) = -3.0;
      h_A(ib, 1, 3) = -2.0;
      h_A(ib, 2, 0) = 2.0;
      h_A(ib, 2, 1) = -1.0;
      h_A(ib, 2, 2) = 1.0;
      h_A(ib, 2, 3) = -3.0;

      // Fill the reference matrix
      h_Ref(ib, 0, 0) = 1.0;
      h_Ref(ib, 0, 1) = -1.5;
      h_Ref(ib, 0, 2) = 1.0;
      h_Ref(ib, 0, 3) = 4.5;
      h_Ref(ib, 1, 0) = -1.0;
      h_Ref(ib, 1, 1) = 4.0;
      h_Ref(ib, 1, 2) = 3.0;
      h_Ref(ib, 1, 3) = 7.0;
      h_Ref(ib, 2, 0) = 2.0;
      h_Ref(ib, 2, 1) = 3.5;
      h_Ref(ib, 2, 2) = 10.0;
      h_Ref(ib, 2, 3) = 10.5;
    }
    Kokkos::deep_copy(A, h_A);
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    // Compute A = A + alpha * (x * y^T)
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "ger", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_A = Kokkos::subview(A, ib, Kokkos::ALL, Kokkos::ALL);
          auto sub_x = Kokkos::subview(x, ib, Kokkos::ALL);
          auto sub_y = Kokkos::subview(y, ib, Kokkos::ALL);

          // A = A + alpha * (x * y^T)
          KokkosBatched::SerialGer<KokkosBatched::Trans::Transpose>::invoke(alpha, sub_x, sub_y, sub_A);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_A, A);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          if (Kokkos::abs(h_A(ib, i, j) - h_Ref(ib, i, j)) > eps) correct = false;
        }
      }
    }

    if (correct) {
      std::cout << "ger works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
