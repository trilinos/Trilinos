
#ifndef _kokkos_saxpy_hpp_
#define _kokkos_saxpy_hpp_

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_HAVE_CUDA
void call_kokkos_saxpy_cuda(size_t N, Kokkos::View<float*,Kokkos::Cuda> y, float alpha, Kokkos::View<float*,Kokkos::Cuda> x);
#endif

#ifdef KOKKOS_HAVE_OPENMP
void call_kokkos_saxpy_openmp(size_t N, Kokkos::View<float*,Kokkos::OpenMP> y, float alpha, Kokkos::View<float*,Kokkos::OpenMP> x);
#endif

#endif

