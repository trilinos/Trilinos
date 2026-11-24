// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Spr.hpp>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;

/// \brief Example of batched spr
/// Perform a symmetric rank-1 update of matrix A by vector x, where
///   A: [[1, -3, -2,  0],
///       [0,  3, -1, -2],
///       [0,  0,  9,  5],
///       [0,  0,  0, 27]]
///   A_packed: [1, -3, -2, 0, 3, -1, -2, 9, 5, 27]
///   x: [1, 2, 3, 4]
///   alpha: 1.5
///
/// After, A = A + alpha * (x * x^T), it will give
///   A: [2.5,  0.,  9.,  2.5,  8.,  22.5,  6.,  10.,  23.,  51.]
///
int main(int /*argc*/, char** /*argv*/) {
  Kokkos::initialize();
  {
    using View2DType = Kokkos::View<double**, ExecutionSpace>;
    const int Nb = 10, n = 4;
    const int packed_size = n * (n + 1) / 2;

    // Matrix A in packed storage
    View2DType ap("ap", Nb, packed_size), Ref("Ref", Nb, packed_size);

    // Vector x
    View2DType x("x", Nb, n);

    // Initialize A and x
    auto h_ap  = Kokkos::create_mirror_view(ap);
    auto h_x   = Kokkos::create_mirror_view(x);
    auto h_Ref = Kokkos::create_mirror_view(Ref);

    // Upper triangular matrix
    for (int ib = 0; ib < Nb; ib++) {
      // Fill vector x
      for (int i = 0; i < n; i++) {
        h_x(ib, i) = i + 1;
      }

      // Fill the matrix A
      h_ap(ib, 0) = 1;
      h_ap(ib, 1) = -3;
      h_ap(ib, 2) = 3;
      h_ap(ib, 3) = -2;
      h_ap(ib, 4) = -1;
      h_ap(ib, 5) = 9;
      h_ap(ib, 6) = 0;
      h_ap(ib, 7) = -2;
      h_ap(ib, 8) = 5;
      h_ap(ib, 9) = 27;

      h_Ref(ib, 0) = 2.5;
      h_Ref(ib, 1) = 0;
      h_Ref(ib, 2) = 9;
      h_Ref(ib, 3) = 2.5;
      h_Ref(ib, 4) = 8;
      h_Ref(ib, 5) = 22.5;
      h_Ref(ib, 6) = 6;
      h_Ref(ib, 7) = 10;
      h_Ref(ib, 8) = 23;
      h_Ref(ib, 9) = 51;
    }
    Kokkos::deep_copy(ap, h_ap);
    Kokkos::deep_copy(x, h_x);

    // Compute A = A + alpha * (x * x^T)
    const double alpha = 1.5;
    ExecutionSpace exec;
    using policy_type = Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<int>>;
    policy_type policy{exec, 0, Nb};
    Kokkos::parallel_for(
        "syr", policy, KOKKOS_LAMBDA(int ib) {
          auto sub_ap = Kokkos::subview(ap, ib, Kokkos::ALL);
          auto sub_x  = Kokkos::subview(x, ib, Kokkos::ALL);

          // A = A + alpha * (x * x^T)
          KokkosBatched::SerialSpr<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose>::invoke(alpha, sub_x,
                                                                                                        sub_ap);
        });

    // Confirm that the results are correct
    Kokkos::deep_copy(h_ap, ap);
    bool correct = true;
    double eps   = 1.0e-12;
    for (int ib = 0; ib < Nb; ib++) {
      for (int i = 0; i < packed_size; i++) {
        if (Kokkos::abs(h_ap(ib, i) - h_Ref(ib, i)) > eps) correct = false;
      }
    }

    if (correct) {
      std::cout << "spr works correctly!" << std::endl;
    }
  }
  Kokkos::finalize();
}
