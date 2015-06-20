
#ifndef _create_kokkos_view_hpp_
#define _create_kokkos_view_hpp_

#include <Kokkos_Core.hpp>
#include <string>

template<class KokkosDevice>
Kokkos::View<float*,KokkosDevice> create_kokkos_view_float1D(const std::string& name, size_t N, float init_value);

#endif

