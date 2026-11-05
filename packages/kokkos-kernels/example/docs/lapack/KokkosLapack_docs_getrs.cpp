// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosLapack_getrs.hpp>
#include <iostream>

int main(void) {
  bool correct = true;
  Kokkos::initialize();
  {
    using execution_space = Kokkos::DefaultExecutionSpace;
    using KAT             = KokkosKernels::ArithTraits<double>;

    Kokkos::View<double**, Kokkos::LayoutLeft> A("A", 3, 3);
    Kokkos::View<double**, Kokkos::LayoutLeft> B("B", 3, 2);
    Kokkos::View<int*, Kokkos::LayoutLeft> Ipiv("Ipiv", 3);
    Kokkos::View<int*, Kokkos::LayoutLeft> Info("Info", 1);

    {
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
    }

    auto h_B  = Kokkos::create_mirror_view(B);
    h_B(0, 0) = 0;
    h_B(0, 1) = 0;
    h_B(1, 0) = 4;
    h_B(1, 1) = -2;
    h_B(2, 0) = 4;
    h_B(2, 1) = 4;

    execution_space space{};
    KokkosLapack::getrs(space, "N", A, Ipiv, B, Info);
    Kokkos::fence();
    Kokkos::deep_copy(h_B, B);

    auto h_Ipiv = Kokkos::create_mirror_view(Ipiv);
    Kokkos::deep_copy(h_Ipiv, Ipiv);

    if (KAT::abs(h_B(0, 0) - 1) > 10 * KAT::epsilon()) {
      correct = false;
    }
    if (KAT::abs(h_B(1, 0) - 2) / 2 > 10 * KAT::epsilon()) {
      correct = false;
    }
    if (KAT::abs(h_B(2, 0) - 3) / 3 > 10 * KAT::epsilon()) {
      correct = false;
    }
    if (KAT::abs(h_B(0, 1) - 3) / 3 > 10 * KAT::epsilon()) {
      correct = false;
    }
    if (KAT::abs(h_B(1, 1) - 2) / 2 > 10 * KAT::epsilon()) {
      correct = false;
    }
    if (KAT::abs(h_B(2, 1) - 3) / 3 > 10 * KAT::epsilon()) {
      correct = false;
    }
    if (correct) std::cout << "KokkosLapack::getrs() returned correct results!" << std::endl;
  }
  Kokkos::finalize();
}
