#ifndef BLASWRAPPER_IAMAX_SPEC_HPP_
#define BLASWRAPPER_IAMAX_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace BlasWrapper {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class RMV, class XMV, int rank = XMV::rank>
struct iamax_eti_spec_avail {
  enum : bool { value = false };
};

// Specialization struct which defines whether a tpl exists
template<class RV, class XMV, int Xrank = XMV::Rank>
struct iamax_tpl_spec_avail {
  enum : bool { value = true };
};

// Unification layer
template<class RMV, class XMV, int rank = XMV::rank,
         bool tpl_spec_avail = iamax_tpl_spec_avail<RMV,XMV>::value,
         bool eti_spec_avail = iamax_eti_spec_avail<RMV,XMV>::value>
struct Iamax {
  static void iamax (const RMV& R, const XMV& X);
};
}
}

namespace BlasWrapper {
namespace Impl {

template<class RV, class XV>
inline void iamax_print_specialization() {
    #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
      #ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
        printf("BlasWrapper::iamax<> TPL cuBLAS specialization for < %s , %s >\n",typeid(RV).name(),typeid(XV).name());
      #else
        #ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
          printf("BlasWrapper::iamax<> TPL Blas specialization for < %s , %s >\n",typeid(RV).name(),typeid(XV).name());
        #endif        
      #endif
    #endif
}

}
}

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

extern "C" {
  //int isamax_( const int* N, const float* x, const int* x_inc);
  //int idamax_( const int* N, const double* x, const int* x_inc);
  //int icamax_( const int* N, const std::complex<float>* x, const int* x_inc);
  int izamax_( const int* N, const std::complex<double>* x, const int* x_inc);
}

namespace BlasWrapper {
namespace Impl {

#define BLASWRAPPER_ZIAMAX_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true,ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<unsigned long, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::iamax[TPL_BLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if (numElems == 0) { R() = 0; return; } \
    if (numElems < static_cast<size_type> (INT_MAX)) { \
      iamax_print_specialization<RV,XV>(); \
      int N = static_cast<int> (numElems); \
      const int XST = X.stride(0); \
      const int LDX = (XST == 0) ? 1 : XST; \
      int idx = izamax_(&N,reinterpret_cast<const std::complex<double>*>(X.data()),&LDX); \
      R() = static_cast<size_type>(idx-1); \
    } else { \
      R() = 0; return; \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

BLASWRAPPER_ZIAMAX_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}
}

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace BlasWrapper {
namespace Impl {

#define BLASWRAPPER_ZIAMAX_TPL_SPEC_DECL_CUBLAS( INDEX_TYPE, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Iamax< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true,ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void iamax (RV& R, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("BlasWrapper::iamax[TPL_CUBLAS,complex<double>]"); \
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
      R() = static_cast<size_type>(idx-1); \
    } else { \
      Kokkos::deep_copy (R, 0); return; \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

BLASWRAPPER_ZIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned long, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
BLASWRAPPER_ZIAMAX_TPL_SPEC_DECL_CUBLAS( unsigned int, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

}
}

#endif

#endif // BLASWRAPPER_IAMAX_SPEC_HPP_
