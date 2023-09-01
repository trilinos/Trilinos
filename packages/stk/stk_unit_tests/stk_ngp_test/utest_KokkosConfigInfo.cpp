#include <Kokkos_Core.hpp>
#include <iostream>
#include <stk_ngp_test/ngp_test.hpp>

NGP_TEST(KokkosConfigInfo, print_device_information) {
  std::cout << "Kokkos::DefaultExecutionSpace::device_type::execution_space::name() returns "
            <<  Kokkos::DefaultExecutionSpace::device_type::execution_space::name() << std::endl;
  std::cout << "Kokkos::DefaultExecutionSpace::device_type::execution_space::memory_space::name() returns "
            <<  Kokkos::DefaultExecutionSpace::device_type::execution_space::memory_space::name() << std::endl;
#ifdef KOKKOS_ENABLE_CUDA
#ifdef KOKKOS_ENABLE_CUDA_UVM
  std::cout << "KOKKOS_ENABLE_CUDA_UVM IS defined" << std::endl;
#else
  std::cout << "KOKKOS_ENABLE_CUDA_UVM IS NOT defined" << std::endl;
#endif
#endif

  std::cout << "Kokkos::DefaultHostExecutionSpace().concurrency() returns "
            <<  Kokkos::DefaultHostExecutionSpace().concurrency() << std::endl;
}

