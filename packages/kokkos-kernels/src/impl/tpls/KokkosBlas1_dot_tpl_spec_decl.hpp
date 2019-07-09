/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS1_DOT_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_HPP_


namespace KokkosBlas {
namespace Impl {

namespace {
  template<class RV, class XV, class YV>
  inline void dot_print_specialization() {
      #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
        printf("KokkosBlas1::dot<> TPL Blas specialization for < %s , %s , %s >\n",typeid(RV).name(),typeid(XV).name(),typeid(YV).name());
      #endif
  }
}

}
}

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      int N = numElems; \
      int one = 1; \
      R() = HostBlas<double>::dot(N,X.data(),one,Y.data(),one);    \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      int N = numElems; \
      int one = 1; \
      R() = HostBlas<float>::dot(N,X.data(),one,Y.data(),one);    \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      int N = numElems; \
      int one = 1; \
      R() = HostBlas<std::complex<double> >::dot                       \
        (N,                                                             \
         reinterpret_cast<const std::complex<double>* >(X.data()),one, \
         reinterpret_cast<const std::complex<double>* >(Y.data()),one); \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      int N = numElems; \
      int one = 1; \
      R() = HostBlas<std::complex<float> >::dot                       \
        (N,                                                             \
         reinterpret_cast<const std::complex<float>* >(X.data()),one, \
         reinterpret_cast<const std::complex<float>* >(Y.data()),one); \
    } else {                                                           \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}
}

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasDdot(s.handle, N, X.data(), one, Y.data(), one, &R()); \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasSdot(s.handle, N, X.data(), one, Y.data(), one, &R()); \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasZdotc(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), one, reinterpret_cast<const cuDoubleComplex*>(Y.data()), one, reinterpret_cast<cuDoubleComplex*>(&R())); \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Dot< \
Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void dot (RV& R, const XV& X, const XV& Y) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      dot_print_specialization<RV,XV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasCdotc(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), one, reinterpret_cast<const cuComplex*>(Y.data()), one, reinterpret_cast<cuComplex*>(&R())); \
    } else { \
      Dot<RV,XV,XV,1,1,false,ETI_SPEC_AVAIL>::dot(R,X,Y); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

}
}

#endif

#endif
