#include <create_kokkos_view.cpp>

template Kokkos::View<float*,Kokkos::Cuda> create_kokkos_view_float1D(const std::string& name, size_t N, float init_value);

