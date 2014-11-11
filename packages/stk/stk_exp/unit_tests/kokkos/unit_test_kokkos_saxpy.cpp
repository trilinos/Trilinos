#include <gtest/gtest.h>

#include <stk_util/stk_config.h>

// restrict this file to only build if KokkosCore is enabled
#ifdef HAVE_STK_KokkosCore

#include <Kokkos_Core.hpp>

#include <iostream>

#if defined(KOKKOS_HAVE_CUDA)
#define KOKKOS_DEVICE Kokkos::Cuda
#define CALL_KOKKOS_SAXPY_FUNCTION call_kokkos_saxpy_cuda
#elif defined(KOKKOS_HAVE_OPENMP)
#define KOKKOS_DEVICE Kokkos::OpenMP
#define CALL_KOKKOS_SAXPY_FUNCTION call_kokkos_saxpy_openmp
#else
#define KOKKOS_DEVICE Kokkos::Serial
#define CALL_KOKKOS_SAXPY_FUNCTION call_kokkos_saxpy_serial
#endif

#include <copy_kokkos_memory.hpp>
#include <create_kokkos_view.hpp>
#include <saxpy_kokkos.hpp>

TEST(stk_exp_kokkos, kokkos_saxpy)
{
  KOKKOS_DEVICE::initialize();
#if defined(KOKKOS_HAVE_CUDA)
  KOKKOS_DEVICE::print_configuration(std::cout);
#endif
  const size_t N = 1000000;
  const float alpha = 0.5;

  Kokkos::View<float*,KOKKOS_DEVICE> x = create_kokkos_view_float1D<KOKKOS_DEVICE>("x", N, 2.0);
  Kokkos::View<float*,KOKKOS_DEVICE> y = create_kokkos_view_float1D<KOKKOS_DEVICE>("y", N, 1.0);

  CALL_KOKKOS_SAXPY_FUNCTION(N, y, alpha, x);

  std::vector<float> host_y;
  copy_kokkos_device_memory_to_host(y, host_y);

  double expected_value = 2.0;
  EXPECT_EQ(expected_value, host_y[0]);

  KOKKOS_DEVICE::finalize();
}

#endif

