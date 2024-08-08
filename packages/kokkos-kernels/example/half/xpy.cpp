//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosKernels_default_types.hpp"
#include <iostream>

template <class ViewType>
struct Functor_xpy {
  ViewType x, y;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const { x(i) += y(i); }
};

template <class ScalarType, class DeviceType, class LayoutType>
void do_xpy(size_t n, bool time_only = false) {
  using namespace Kokkos;
  using ExecutionSpace      = typename DeviceType::execution_space;
  using ViewType            = View<ScalarType *, LayoutType, DeviceType>;
  using ReferenceScalarType = double;

  ViewType x("x", n);
  ViewType y("y", n);
  View<ReferenceScalarType *, LayoutType, DeviceType> x_rand("x_rand", n);
  View<ReferenceScalarType *, LayoutType, DeviceType> y_rand("y_rand", n);

  View<ReferenceScalarType *, LayoutType, HostSpace> expected("expected", n);
  View<ReferenceScalarType *, LayoutType, HostSpace> relative_error("relative_error", n);
  typename ViewType::HostMirror x_host = create_mirror_view(x);
  typename ViewType::HostMirror y_host = create_mirror_view(y);
  // TODO: Report segfault in random_pool creation with:
  // typename ViewType::HostMirror y_host = create_mirror_view(y_host);

  Random_XorShift64_Pool<ExecutionSpace> random_pool(12345);
  fill_random(x_rand, random_pool, ReferenceScalarType(1.0), ReferenceScalarType(2.0));
  fill_random(y_rand, random_pool, ReferenceScalarType(1.0), ReferenceScalarType(2.0));
  ExecutionSpace().fence();

  deep_copy(x, x_rand);
  deep_copy(y, y_rand);
  ExecutionSpace().fence();

  deep_copy(x_host, x);
  deep_copy(y_host, y);
  ExecutionSpace().fence();

  Functor_xpy<ViewType> xpy;
  xpy.x = x;
  xpy.y = y;
  Timer timer;
  parallel_for("xpy", n, xpy);
  ExecutionSpace().fence();
  double s = timer.seconds();

  if (!time_only) {
    for (size_t i = 0; i < n; i++)
      expected(i) = static_cast<ReferenceScalarType>(y_host(i)) + static_cast<ReferenceScalarType>(x_host(i));
  }

  deep_copy(x_host, x);
  ExecutionSpace().fence();

  std::cout << "n: " << n << ", " << typeid(ScalarType).name() << " Runtime(s): " << s << std::endl;

  if (!time_only) {
    std::cout << "n: " << n << ", " << typeid(ScalarType).name() << " Relative Errors:" << std::endl;
    for (size_t i = 0; i < n; i++) {
      std::cout << ", " << std::abs(expected(i) - x_host(i)) / expected(i) << std::endl;
    }
    std::cout << std::endl << std::endl;
  }
}

int main(int argc, char **argv) {
  Kokkos::initialize();
  if (argc < 2) {
    std::cout << "./" << argv[0] << " N:Z TIME_ONLY:{0,1}" << std::endl;
    Kokkos::finalize();
    return 1;
  }
  using LayoutType = Kokkos::LayoutLeft;
  using DeviceType = default_device;
  size_t n         = atoi(argv[1]);
  bool time_only   = static_cast<bool>(atoi(argv[2]));
  do_xpy<float, DeviceType, LayoutType>(n, time_only);
  do_xpy<Kokkos::Experimental::half_t, DeviceType, LayoutType>(n, time_only);
  do_xpy<Kokkos::Experimental::bhalf_t, DeviceType, LayoutType>(n, time_only);
  Kokkos::finalize();
  return 0;
}