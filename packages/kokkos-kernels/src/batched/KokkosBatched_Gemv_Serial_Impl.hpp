#ifndef __KOKKOSBATCHED_GEMV_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_GEMV_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemv_Serial_Internal.hpp"

namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Serial Impl
    /// ===========

    ///
    /// Implemented:
    /// NT, T
    ///
    /// Not yet implemented
    /// CT

    ///
    /// NT
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemv<Trans::NoTranspose,Algo::Gemv::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      typedef typename yViewType::value_type vector_type;
      typedef typename vector_type::value_type value_type;

      const int
        m = A.dimension(0),
        n = 1,
        k = A.dimension(1),
        vl = vector_type::vector_length;

      // no error check
      int r_val = 0;
      if (A.stride_0() == 1) {
        cblas_dgemm_compact(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_1(), 
                            (const double*)x.data(), x.stride_0(), 
                            beta,
                            (double*)y.data(), y.stride_0(),
                            (MKL_INT)vl, (MKL_INT)1);
      } else if (A.stride_1() == 1) {
        cblas_dgemm_compact(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (const double*)x.data(), x.stride_0(), 
                            beta,
                            (double*)y.data(), y.stride_0(),
                            (MKL_INT)vl, (MKL_INT)1);
      } else {
        r_val = -1;
      }
      return r_val;
    }
#endif

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemv<Trans::NoTranspose,Algo::Gemv::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      return SerialGemvInternal<Algo::Gemv::Unblocked>::
        invoke(A.extent(0), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1(),
               x.data(), x.stride_0(),
               beta,
               y.data(), y.stride_0());
    }
    
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemv<Trans::NoTranspose,Algo::Gemv::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      return SerialGemvInternal<Algo::Gemv::Blocked>::
        invoke(A.extent(0), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1(),
               x.data(), x.stride_0(),
               beta,
               y.data(), y.stride_0());
    }
    ///
    /// T
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemv<Trans::Transpose,Algo::Gemv::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      typedef typename yViewType::value_type vector_type;
      typedef typename vector_type::value_type value_type;

      const int
        m = A.dimension(0),
        n = 1,
        k = A.dimension(1),
        vl = vector_type::vector_length;

      // no error check
      int r_val = 0;
      if (A.stride_0() == 1) {
        cblas_dgemm_compact(CblasColMajor, CblasTrans, CblasNoTrans,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_1(), 
                            (const double*)x.data(), x.stride_0(), 
                            beta,
                            (double*)y.data(), y.stride_0(),
                            (MKL_INT)vl, (MKL_INT)1);
      } else if (A.stride_1() == 1) {
        cblas_dgemm_compact(CblasRowMajor, CblasTrans, CblasNoTrans,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (const double*)x.data(), x.stride_0(), 
                            beta,
                            (double*)y.data(), y.stride_0(),
                            (MKL_INT)vl, (MKL_INT)1);
      } else {
        r_val = -1;
      }
    }
#endif

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemv<Trans::Transpose,Algo::Gemv::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      return SerialGemvInternal<Algo::Gemv::Unblocked>::
        invoke(A.extent(1), A.extent(0),
               alpha, 
               A.data(), A.stride_1(), A.stride_0(),
               x.data(), x.stride_0(),
               beta,
               y.data(), y.stride_0());
    }
    
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemv<Trans::Transpose,Algo::Gemv::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      return SerialGemvInternal<Algo::Gemv::Blocked>::
        invoke(A.extent(1), A.extent(0),
               alpha, 
               A.data(), A.stride_1(), A.stride_0(),
               x.data(), x.stride_0(),
               beta,
               y.data(), y.stride_0());
    }

  }
}

#endif
