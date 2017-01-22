#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSPACES_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSPACES_H_

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_HAVE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif

namespace ngp {

#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP   ExecSpace ;
#elif KOKKOS_HAVE_CUDA
  typedef Kokkos::Cuda     ExecSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
#endif

#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP   HostExecSpace ;
#elif KOKKOS_HAVE_CUDA
  typedef Kokkos::Serial   HostExecSpace ;
#else
  typedef Kokkos::Serial   HostExecSpace ;
#endif

#ifdef KOKKOS_HAVE_OPENMP
   typedef Kokkos::OpenMP       MemSpace;
#elif KOKKOS_HAVE_CUDA
   typedef Kokkos::CudaSpace    MemSpace;
#else
   typedef Kokkos::HostSpace    MemSpace;
#endif

#ifdef KOKKOS_HAVE_OPENMP
typedef Kokkos::OpenMP       UVMMemSpace;
#elif KOKKOS_HAVE_CUDA
typedef Kokkos::CudaUVMSpace UVMMemSpace;
#else
typedef Kokkos::HostSpace    UVMMemSpace;
#endif

typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSPACES_H_ */
