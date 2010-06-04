#include "Tpetra_Vector.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

// #include "Tpetra_ExplicitInstantiationHelpers.hpp"

#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOS_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
#  include <Kokkos_TPINode.hpp>
#endif
#if defined(HAVE_KOKKOS_THRUST)
#  include <Kokkos_ThrustGPUNode.hpp>
#endif

#include "Tpetra_Vector_def.hpp"

namespace Tpetra {

  TPETRA_VECTOR_INSTANT(int,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  TPETRA_VECTOR_INSTANT(int,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
  TPETRA_VECTOR_INSTANT(int,int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_FLOAT)
  TPETRA_VECTOR_INSTANT(int,int,int,Kokkos::ThrustGPUNode)
#endif

#if defined(HAVE_TPETRA_INST_FLOAT)
  TPETRA_VECTOR_INSTANT(float,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  TPETRA_VECTOR_INSTANT(float,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_VECTOR_INSTANT(float,int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_FLOAT)
    TPETRA_VECTOR_INSTANT(float,int,int,Kokkos::ThrustGPUNode)
#endif
#endif

#if defined(HAVE_TPETRA_INST_DOUBLE)
  TPETRA_VECTOR_INSTANT(double,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  TPETRA_VECTOR_INSTANT(double,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_VECTOR_INSTANT(double,int,int,Kokkos::TPINode)
#endif
#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_DOUBLE)
    TPETRA_VECTOR_INSTANT(double,int,int,Kokkos::ThrustGPUNode)
#endif
#endif

#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  TPETRA_VECTOR_INSTANT(std::complex<float>,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  TPETRA_VECTOR_INSTANT(std::complex<float>,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_VECTOR_INSTANT(std::complex<float>,int,int,Kokkos::TPINode)
#endif
// no complex on GPU support for now
//#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_DOUBLE)
//    TPETRA_VECTOR_INSTANT(double,int,int,Kokkos::ThrustGPUNode)
//#endif
#endif

#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  TPETRA_VECTOR_INSTANT(std::complex<double>,int,int,Kokkos::SerialNode)
#if defined(HAVE_KOKKOS_TBB)
  TPETRA_VECTOR_INSTANT(std::complex<double>,int,int,Kokkos::TBBNode)
#endif
#if defined(HAVE_KOKKOS_THREADPOOL)
    TPETRA_VECTOR_INSTANT(std::complex<double>,int,int,Kokkos::TPINode)
#endif
// no complex on GPU support for now
//#if defined(HAVE_KOKKOS_THRUST) && defined(HAVE_KOKKOS_CUDA_DOUBLE)
//    TPETRA_VECTOR_INSTANT(double,int,int,Kokkos::ThrustGPUNode)
//#endif
#endif


} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
