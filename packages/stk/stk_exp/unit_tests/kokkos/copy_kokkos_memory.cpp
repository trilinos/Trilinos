
#include <Kokkos_Core.hpp>
#include <vector>

#if defined(KOKKOS_HAVE_CUDA) && defined(__CUDACC__)
#define KOKKOS_DEVICE Kokkos::Cuda
#elif defined(KOKKOS_HAVE_OPENMP)
#define KOKKOS_DEVICE Kokkos::OpenMP
#else
#define KOKKOS_DEVICE Kokkos::Serial
#endif

void copy_kokkos_device_memory_to_host(Kokkos::View<float*,KOKKOS_DEVICE>& y, std::vector<float>& host_y)
{
  Kokkos::View<float*,KOKKOS_DEVICE>::HostMirror host_view = Kokkos::create_mirror_view(y);

  Kokkos::deep_copy(host_view, y);

  size_t N = y.dimension(0);
  host_y.resize(N);

  for(size_t i=0; i<N; ++i) {
    host_y[i] = host_view(i);
  }
}

