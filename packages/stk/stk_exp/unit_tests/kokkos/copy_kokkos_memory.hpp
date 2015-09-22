#ifndef _kokkos_copy_memory_hpp_
#define _kokkos_copy_memory_hpp_

#include <Kokkos_Core.hpp>

void copy_kokkos_device_memory_to_host(Kokkos::View<float*,KOKKOS_DEVICE>& y, std::vector<float>& host_y);

#endif

