// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosLapack_getrf.hpp>
#include <iostream>

int main(void) {
  bool correct = true;
  Kokkos::initialize();
  {
    using execution_space = Kokkos::DefaultExecutionSpace;
    using KAT             = KokkosKernels::ArithTraits<double>;

    Kokkos::View<double**, Kokkos::LayoutLeft> A("A", 3, 3);
    Kokkos::View<int*, Kokkos::LayoutLeft> Ipiv("tau", 3);
    Kokkos::View<int*, Kokkos::LayoutLeft> Info("Info", 1);

    auto h_A  = Kokkos::create_mirror_view(A);
    h_A(0, 0) = -1;
    h_A(0, 1) = 2;
    h_A(0, 2) = -1;
    h_A(1, 0) = 2;
    h_A(1, 1) = -1;
    h_A(1, 2) = 0;
    h_A(2, 0) = 0;
    h_A(2, 1) = -1;
    h_A(2, 2) = 2;
    Kokkos::deep_copy(A, h_A);

    execution_space space{};
    KokkosLapack::getrf(space, A, Ipiv, Info);
    Kokkos::fence();
    Kokkos::deep_copy(h_A, A);

    auto h_Ipiv = Kokkos::create_mirror_view(Ipiv);
    Kokkos::deep_copy(h_Ipiv, Ipiv);

    Kokkos::View<double**, Kokkos::DefaultHostExecutionSpace> analytic("result", 3, 3);
    analytic(0, 0) = 2.;
    analytic(0, 1) = -1.;
    analytic(0, 2) = 0.;
    analytic(1, 0) = -1. / 2.;
    analytic(1, 1) = 3. / 2.;
    analytic(1, 2) = -1;
    analytic(2, 0) = 0.;
    analytic(2, 1) = -2. / 3.;
    analytic(2, 2) = 4. / 3.;

    for (int rowIdx = 0; rowIdx < A.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < A.extent_int(1); ++colIdx) {
        if (KAT::abs(h_A(rowIdx, colIdx) - analytic(rowIdx, colIdx)) >
            KAT::abs(analytic(rowIdx, colIdx)) * 10 * KAT::epsilon()) {
          std::cout << "Error: h_A(" << rowIdx << ", " << colIdx << ")=" << h_A(rowIdx, colIdx) << ", analytic is "
                    << analytic(rowIdx, colIdx) << std::endl;
          correct = false;
        }
      }
    }
    if (correct) std::cout << "KokkosLapack::getrf() returned correct results!" << std::endl;
  }
  Kokkos::finalize();
}
