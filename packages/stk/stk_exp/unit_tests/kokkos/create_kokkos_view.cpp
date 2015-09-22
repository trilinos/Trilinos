
#include <Kokkos_Core.hpp>
#include <string>

//Init Functor
template<class KokkosDevice>
struct Init {
  typedef KokkosDevice device_type;
  Kokkos::View<float*,KokkosDevice> x;
  float init_value;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    x(i) = init_value;
  }
};

template<class KokkosDevice>
Kokkos::View<float*,KokkosDevice> create_kokkos_view_float1D(const std::string& name, size_t N, float init_value)
{
  Kokkos::View<float*,KokkosDevice> x(name, N);

  Init<KokkosDevice> init = {x, init_value};

  Kokkos::parallel_for(N, init);

  return x;
}

//explicit template instantiations:
//
#ifdef KOKKOS_HAVE_OPENMP
template Kokkos::View<float*,Kokkos::OpenMP> create_kokkos_view_float1D(const std::string& name, size_t N, float init_value);
#endif

template Kokkos::View<float*,Kokkos::Serial> create_kokkos_view_float1D(const std::string& name, size_t N, float init_value);

