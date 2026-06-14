// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosLapack_potrf.hpp>
#include <KokkosLapack_potrs.hpp>
#include <KokkosKernels_TestMatrixUtils.hpp>

#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    // Solve A*X = B for a 4x4 symmetric positive definite (SPD) matrix A
    // and two right-hand sides.
    //
    // Strategy:
    //   1. Choose a known X_exact and compute B = A * X_exact.
    //   2. Factor A with potrf (Cholesky): A = L * L^H.
    //   3. Solve with potrs:              L * L^H * X = B  ->  X overwrites B.
    //   4. Verify X ≈ X_exact.
    //
    // A (diagonally dominant, SPD):        X_exact:    B = A * X_exact:
    //   [10  1  2  1]                        [1  0]      [22   8]
    //   [ 1 10  3  2]                        [2  1]      [38  22]
    //   [ 2  3 10  1]                        [3  2]      [42  26]
    //   [ 1  2  1 10]                        [4  3]      [48  34]

    const int N = 4, nrhs = 2;
    using ViewType = Kokkos::View<double**, Kokkos::LayoutLeft>;

    ViewType A("A", N, N);
    ViewType B("B", N, nrhs);
    // clang-format off
    std::vector<std::vector<double>> A_data = {
        {10,  1,  2,  1},
        { 1, 10,  3,  2},
        { 2,  3, 10,  1},
        { 1,  2,  1, 10}};
    std::vector<std::vector<double>> B_data = {
        {22,  8},
        {38, 22},
        {42, 26},
        {48, 34}};
    // clang-format on

    Test::fill_view_from_fixture(A, A_data);
    Test::fill_view_from_fixture(B, B_data);

    // Step 1: Cholesky factorization — A is overwritten with the lower factor L
    KokkosLapack::potrf("L", A);

    // Step 2: Triangular solve — B is overwritten with the solution X
    Kokkos::View<const double**, Kokkos::LayoutLeft> Aconst(A);
    KokkosLapack::potrs("L", Aconst, B);

    auto h_B = Kokkos::create_mirror_view(B);
    Kokkos::deep_copy(h_B, B);

    // X_exact col-0 = [1, 2, 3, 4],  col-1 = [0, 1, 2, 3]
    std::cout << "Solution X (col 0 expected [1,2,3,4], col 1 expected [0,1,2,3]):\n";
    for (int i = 0; i < N; ++i)
      std::cout << "  X(" << i << ",0) = " << h_B(i, 0) << "  X(" << i << ",1) = " << h_B(i, 1) << "\n";
  }
  Kokkos::finalize();
}
