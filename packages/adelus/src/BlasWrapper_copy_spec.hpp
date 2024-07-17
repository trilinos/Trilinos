/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef BLASWRAPPER_COPY_SPEC_HPP_
#define BLASWRAPPER_COPY_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace BlasWrapper {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class XMV, class YMV, int rank = YMV::rank>
struct copy_eti_spec_avail {
  enum : bool { value = false };
};

// Specialization struct which defines whether a tpl exists
template<class XMV, class YMV, int rank = YMV::rank>
struct copy_tpl_spec_avail {
  enum : bool { value = true };
};

// Unification layer
template<class XMV, class YMV, int rank = YMV::rank,
         bool tpl_spec_avail = copy_tpl_spec_avail<XMV,YMV>::value,
         bool eti_spec_avail = copy_eti_spec_avail<XMV,YMV>::value>
struct Copy {
  static void copy (const XMV& X, const YMV& Y);
};
}
}


namespace BlasWrapper {
namespace Impl {

template<class XV, class YV>
inline void copy_print_specialization() {
    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
      #if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
        printf("BlasWrapper::copy<> TPL cuBLAS specialization for < %s , %s >\n",typeid(XV).name(),typeid(YV).name());
      #elif defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
        printf("BlasWrapper::copy<> TPL rocBLAS specialization for < %s , %s >\n",typeid(XV).name(),typeid(YV).name());
      #else
        #ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
          printf("BlasWrapper::copy<> TPL BLAS specialization for < %s , %s >\n",typeid(XV).name(),typeid(YV).name());
        #endif
      #endif
    #endif
}

}
}

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

extern "C" {
  void dcopy_( const int* N, const double* x, const int* x_inc, double* y, const int* y_inc);
  void scopy_( const int* N, const float* x, const int* x_inc,  float* y, const int* y_inc);
  void zcopy_( const int* N, const std::complex<double>* x, const int* x_inc, std::complex<double>* y, const int* y_inc);
  void ccopy_( const int* N, const std::complex<float>* x, const int* x_inc,  std::complex<float>* y, const int* y_inc);
}

namespace BlasWrapper {
namespace Impl {

#define BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const double*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const double*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_BLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int YST = Y.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      const int LDY = (YST == 0) ? 1 : YST; \
      dcopy_(&N, X.data(), &LDX, Y.data(), &LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const float*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const float*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_BLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int YST = Y.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      const int LDY = (YST == 0) ? 1 : YST; \
      scopy_(&N, X.data(), &LDX, Y.data(), &LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const Kokkos::complex<double>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<Kokkos::complex<double>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true,ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<Kokkos::complex<double>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_BLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int YST = Y.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      const int LDY = (YST == 0) ? 1 : YST; \
      zcopy_(&N, reinterpret_cast<const std::complex<double>*>(X.data()), &LDX, reinterpret_cast<std::complex<double>*>(Y.data()), &LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const Kokkos::complex<float>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<Kokkos::complex<float>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<Kokkos::complex<float>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_BLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int YST = Y.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      const int LDY = (YST == 0) ? 1 : YST; \
      ccopy_(&N, reinterpret_cast<const std::complex<float>*>(X.data()), &LDX, reinterpret_cast<std::complex<float>*>(Y.data()), &LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HostSpace, false)

BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HostSpace, false)

BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HostSpace, false)

BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::HostSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_BLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HostSpace, false)

}//namespace Impl
}//namespace BlasWrapper

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace BlasWrapper {
namespace Impl {

#define BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const double*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const double*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_CUBLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasDcopy(s.handle, N, X.data(), LDX, Y.data(), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const float*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const float*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_CUBLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasScopy(s.handle, N, X.data(), LDX, Y.data(), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const Kokkos::complex<double>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<Kokkos::complex<double>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true,ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<Kokkos::complex<double>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_CUBLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasZcopy(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(X.data()), LDX, reinterpret_cast<cuDoubleComplex*>(Y.data()), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const Kokkos::complex<float>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<Kokkos::complex<float>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<Kokkos::complex<float>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_CUBLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasCcopy(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), LDX, reinterpret_cast<cuComplex*>(Y.data()), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaSpace, false)

BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaSpace, false)

BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaSpace, false)

BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaSpace, false)

#if defined(KOKKOS_ENABLE_CUDA_UVM)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)

BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)

BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)

BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft,  Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::CudaUVMSpace, false)
#endif

}//namespace Impl
}//namespace BlasWrapper

#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace BlasWrapper {
namespace Impl {

#define BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const double*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const double*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<double*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_ROCBLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::RocBlasSingleton & s = KokkosBlas::Impl::RocBlasSingleton::singleton(); \
      rocblas_dcopy(s.handle, N, X.data(), LDX, Y.data(), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const float*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const float*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<float*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_ROCBLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::RocBlasSingleton & s = KokkosBlas::Impl::RocBlasSingleton::singleton(); \
      rocblas_scopy(s.handle, N, X.data(), LDX, Y.data(), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const Kokkos::complex<double>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<Kokkos::complex<double>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true,ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<Kokkos::complex<double>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_ROCBLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::RocBlasSingleton & s = KokkosBlas::Impl::RocBlasSingleton::singleton(); \
      rocblas_zcopy(s.handle, N, reinterpret_cast<const rocblas_double_complex*>(X.data()), LDX, reinterpret_cast<rocblas_double_complex*>(Y.data()), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( LAYOUTX, LAYOUTY, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Copy< \
Kokkos::View<const Kokkos::complex<float>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<Kokkos::complex<float>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUTX, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef Kokkos::View<Kokkos::complex<float>*, LAYOUTY, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YV; \
  typedef typename XV::size_type size_type; \
  \
  static void copy (const XV& X, const YV& Y) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::copy[TPL_ROCBLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      copy_print_specialization<XV,YV>(); \
      const int N   = static_cast<int> (numElems); \
      const int XST = X.stride(0), LDX = (XST == 0) ? 1 : XST; \
      const int YST = Y.stride(0), LDY = (YST == 0) ? 1 : YST; \
      KokkosBlas::Impl::RocBlasSingleton & s = KokkosBlas::Impl::RocBlasSingleton::singleton(); \
      rocblas_ccopy(s.handle, N, reinterpret_cast<const rocblas_float_complex*>(X.data()), LDX, reinterpret_cast<rocblas_float_complex*>(Y.data()), LDY); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPSpace, false)

BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPSpace, false)

BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPSpace, false)

BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPSpace, false)

BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_DCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)

BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_SCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)

BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_ZCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)

BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutLeft,   Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutRight,  Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutLeft,   Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutRight,  Kokkos::HIPManagedSpace, false)
BLASWRAPPER_CCOPY_TPL_SPEC_DECL_ROCBLAS( Kokkos::LayoutStride, Kokkos::LayoutStride, Kokkos::HIPManagedSpace, false)

}//namespace Impl
}//namespace BlasWrapper

#endif

#endif // BLASWRAPPER_COPY_SPEC_HPP_
