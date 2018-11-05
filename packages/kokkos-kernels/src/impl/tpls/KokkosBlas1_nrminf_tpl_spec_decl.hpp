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

#ifndef KOKKOSBLAS1_NRMINF_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_NRMINF_TPL_SPEC_DECL_HPP_

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

extern "C" int idamax_( const int* N, const double* x, const int* x_inc);
extern "C" int isamax_( const int* N, const float* x, const int* x_inc);
extern "C" int izamax_( const int* N, const std::complex<double>* x, const int* x_inc);
extern "C" int icamax_( const int* N, const std::complex<float>* x, const int* x_inc);

namespace KokkosBlas {
namespace Impl {


namespace {
  template<class RV, class XV>
  inline void nrminf_print_specialization() {
      #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
        printf("KokkosBlas1::nrminf<> TPL Blas specialization for < %s , %s >\n",typeid(RV).name(),typeid(XV).name());
      #endif
  }
}

#define KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0.0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      int N = numElems; \
      int one = 1; \
      int idx = idamax_(&N,X.data(),&one)-1; \
      R() = X(idx); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

#define KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0.0f; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      int N = numElems; \
      int one = 1; \
      int idx = isamax_(&N,X.data(),&one)-1; \
      R() = X(idx); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

#define KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  typedef Kokkos::Details::InnerProductSpaceTraits<Kokkos::complex<double>> IPT; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0.0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      int N = numElems; \
      int one = 1; \
      int idx = izamax_(&N,reinterpret_cast<const std::complex<double>*>(X.data()),&one)-1; \
      R() = IPT::norm(X(idx)); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

#define KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  typedef Kokkos::Details::InnerProductSpaceTraits<Kokkos::complex<float>> IPT; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0.0f; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      int N = numElems; \
      int one = 1; \
      int idx = icamax_(&N,reinterpret_cast<const std::complex<float>*>(X.data()),&one)-1; \
      R() = IPT::norm(X(idx)); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}
}

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0.0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIdamax(s.handle, N, X.data(), one, &idx); \
      Kokkos::deep_copy(R, subview(X,idx-1)); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

#define KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0.0f);; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIsamax(s.handle, N, X.data(), one, &idx); \
      Kokkos::deep_copy(R, subview(X,idx-1)); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

#define KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  typedef Kokkos::Details::InnerProductSpaceTraits<Kokkos::complex<double>> IPT; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0.0); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIzamax(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), one, &idx); \
      Kokkos::complex<double> R_cplx_val {0.0, 0.0}; \
      Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > R_cplx (&R_cplx_val); \
      Kokkos::deep_copy(R_cplx, subview(X,idx-1)); \
      R() = IPT::norm(R_cplx()); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

#define KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct NrmInf< \
Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  typedef Kokkos::Details::InnerProductSpaceTraits<Kokkos::complex<float>> IPT; \
  \
  static void nrminf (RV& R, const XV& X) \
  { \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { Kokkos::deep_copy (R, 0.0f); return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      nrminf_print_specialization<RV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      int idx; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasIcamax(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), one, &idx); \
      Kokkos::complex<float> R_cplx_val {0.0f, 0.0f}; \
      Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > R_cplx (&R_cplx_val); \
      Kokkos::deep_copy(R_cplx, subview(X,idx-1)); \
      R() = IPT::norm(R_cplx()); \
    } else { \
      NrmInf<RV,XV,1,false,ETI_SPEC_AVAIL>::nrminf(R,X); \
    } \
  } \
};

KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

}
}

#endif

#endif
