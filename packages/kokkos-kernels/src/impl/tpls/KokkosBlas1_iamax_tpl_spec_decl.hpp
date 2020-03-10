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

#ifndef KOKKOSBLAS1_IAMAX_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_DECL_HPP_


namespace KokkosBlas {
namespace Impl {
  template<class RV, class XV>
  inline void iamax_print_specialization() {
      #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
        #ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
          printf("KokkosBlas1::iamax<> TPL cuBLAS specialization for < %s , %s >\n",typeid(RV).name(),typeid(XV).name());
        #else
          #ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
            printf("KokkosBlas1::iamax<> TPL Blas specialization for < %s , %s >\n",typeid(RV).name(),typeid(XV).name());
          #endif        
        #endif
      #endif
  }
}
}

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {


#define KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_BLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx = HostBlas<double>::iamax(N,X.data(),LDX);  \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_BLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx = HostBlas<float>::iamax(N,X.data(),LDX);  \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_BLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx = HostBlas<std::complex<double> >::iamax(N,reinterpret_cast<const std::complex<double>*>(X.data()),LDX); \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_BLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx = HostBlas<std::complex<float> >::iamax(N,reinterpret_cast<const std::complex<float>*>(X.data()),LDX); \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}
}

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS( INDEX_TYPE, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIdamax(s.handle, N, X.data(), LDX, &idx); \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
}; \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_DEVICE); \
      cublasIdamax(s.handle, N, X.data(), LDX, reinterpret_cast<int*>(R.data())); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_HOST); \
      Kokkos::fence(); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS( INDEX_TYPE, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0);; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIsamax(s.handle, N, X.data(), LDX, &idx); \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
}; \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0);; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_DEVICE); \
      cublasIsamax(s.handle, N, X.data(), LDX, reinterpret_cast<int*>(R.data())); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_HOST); \
      Kokkos::fence(); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS( INDEX_TYPE, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIzamax(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), LDX, &idx); \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
}; \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_DEVICE); \
      cublasIzamax(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), LDX, reinterpret_cast<int*>(R.data())); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_HOST); \
      Kokkos::fence(); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS( INDEX_TYPE, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIcamax(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), LDX, &idx); \
      R() = static_cast<size_type>(idx); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
}; \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_DEVICE); \
      cublasIcamax(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), LDX, reinterpret_cast<int*>(R.data())); \
      cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_HOST); \
      Kokkos::fence(); \
    } else { \
      Iamax<RV,XV,1,false,ETI_SPEC_AVAIL>::iamax(R,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

}
}

#endif

#endif
