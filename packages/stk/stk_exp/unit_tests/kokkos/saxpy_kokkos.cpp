
#include <Kokkos_Core.hpp>

#if defined(KOKKOS_HAVE_CUDA) && defined(__CUDACC__)
#define KOKKOS_DEVICE Kokkos::Cuda
#define CALL_KOKKOS_SAXPY_FUNCTION call_kokkos_saxpy_cuda
#elif defined(KOKKOS_HAVE_OPENMP)
#define KOKKOS_DEVICE Kokkos::OpenMP
#define CALL_KOKKOS_SAXPY_FUNCTION call_kokkos_saxpy_openmp
#else
#define KOKKOS_DEVICE Kokkos::Serial
#define CALL_KOKKOS_SAXPY_FUNCTION call_kokkos_saxpy_serial
#endif

//Axpy Functor
template<class DEVICE>
struct Axpy {
  typedef DEVICE device_type;
  Kokkos::View<float*,DEVICE> x;
  Kokkos::View<float*,DEVICE> y;
  float alpha;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    y(i) = y(i) + alpha*x(i);
  }
};

void CALL_KOKKOS_SAXPY_FUNCTION(size_t N, Kokkos::View<float*,KOKKOS_DEVICE> y, float alpha, Kokkos::View<float*,KOKKOS_DEVICE> x)
{
  Axpy<KOKKOS_DEVICE> axpy = {x, y, alpha};

  Kokkos::parallel_for(N, axpy);
}

