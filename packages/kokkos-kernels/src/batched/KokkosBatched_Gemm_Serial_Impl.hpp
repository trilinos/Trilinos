#ifndef __KOKKOSBATCHED_GEMM_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_GEMM_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Serial Impl
    /// ===========

    ///
    /// Implemented: 
    /// NT/NT, T/NT, NT/T, T/T
    ///
    /// Not yet immplemented (ConjTranspose):
    /// CT/NT, NT/CT, CT/CT
    ///

    ///
    /// NT/NT
    ///
    
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      typedef typename CViewType::value_type vector_type;
      //typedef typename vector_type::value_type value_type;

      const int
        m = C.dimension(0),
        n = C.dimension(1),
        k = A.dimension(1);

      static_assert(is_vector<vector_type>::value, "value type is not vector type");      
      static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8, 
                    "AVX, AVX2 and AVX512 is supported");
      const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ?  MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

      // no error check
      int r_val = 0;
      if (A.stride_0() == 1 && B.stride_0() == 1 && C.stride_0() == 1) {
        mkl_dgemm_compact(MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_1(), 
                            (const double*)B.data(), B.stride_1(), 
                            beta,
                            (double*)C.data(), C.stride_1(),
                            format, (MKL_INT)vector_type::vector_length);
      } else if (A.stride_1() == 1 && B.stride_1() == 1 && C.stride_1() == 1) {
        mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_NOTRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (const double*)B.data(), B.stride_0(), 
                            beta,
                            (double*)C.data(), C.stride_0(),
                            format, (MKL_INT)vector_type::vector_length);
      } else {
        r_val = -1;
      }
      return r_val;
    }
#endif

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Unblocked>::
        invoke(C.extent(0), C.extent(1), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1(),
               B.data(), B.stride_0(), B.stride_1(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Blocked>::
        invoke(C.extent(0), C.extent(1), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1(),
               B.data(), B.stride_0(), B.stride_1(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    ///
    /// T/NT
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::Transpose,Trans::NoTranspose,Algo::Gemm::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      typedef typename CViewType::value_type vector_type;
      //typedef typename vector_type::value_type value_type;

      const int
        m = C.dimension(0),
        n = C.dimension(1),
        k = A.dimension(0);

      static_assert(is_vector<vector_type>::value, "value type is not vector type");      
      static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                    "AVX, AVX2 and AVX512 is supported");
      const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ?  MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

      // no error check
      int r_val = 0;
      if (A.stride_0() == 1 && B.stride_0() == 1 && C.stride_0() == 1) {
        mkl_dgemm_compact(MKL_COL_MAJOR, MKL_TRANS, MKL_NOTRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_1(), 
                            (const double*)B.data(), B.stride_1(), 
                            beta,
                            (double*)C.data(), C.stride_1(),
                            format, (MKL_INT)vector_type::vector_length);
      } else if (A.stride_1() == 1 && B.stride_1() == 1 && C.stride_1() == 1) {
        mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_TRANS, MKL_NOTRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (const double*)B.data(), B.stride_0(), 
                            beta,
                            (double*)C.data(), C.stride_0(),
                            format, (MKL_INT)vector_type::vector_length);
      } else {
        r_val = -1;
      }
      return r_val;
    }
#endif

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::Transpose,Trans::NoTranspose,Algo::Gemm::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Unblocked>::
        invoke(C.extent(0), C.extent(1), A.extent(0),
               alpha, 
               A.data(), A.stride_1(), A.stride_0(),
               B.data(), B.stride_0(), B.stride_1(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::Transpose,Trans::NoTranspose,Algo::Gemm::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Blocked>::
        invoke(C.extent(0), C.extent(1), A.extent(0),
               alpha, 
               A.data(), A.stride_1(), A.stride_0(),
               B.data(), B.stride_0(), B.stride_1(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    ///
    /// NT/T
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::NoTranspose,Trans::Transpose,Algo::Gemm::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      typedef typename CViewType::value_type vector_type;
      //typedef typename vector_type::value_type value_type;

      const int
        m = C.dimension(0),
        n = C.dimension(1),
        k = A.dimension(1);

      static_assert(is_vector<vector_type>::value, "value type is not vector type");      
      static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                    "AVX, AVX2 and AVX512 is supported");
      const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ?  MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

      // no error check
      int r_val = 0;
      if (A.stride_0() == 1 && B.stride_0() == 1 && C.stride_0() == 1) {
        mkl_dgemm_compact(MKL_COL_MAJOR, MKL_NOTRANS, MKL_TRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_1(), 
                            (const double*)B.data(), B.stride_1(), 
                            beta,
                            (double*)C.data(), C.stride_1(),
                            format, (MKL_INT)vector_type::vector_length);
      } else if (A.stride_1() == 1 && B.stride_1() == 1 && C.stride_1() == 1) {
        mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_TRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (const double*)B.data(), B.stride_0(), 
                            beta,
                            (double*)C.data(), C.stride_0(),
                            format, (MKL_INT)vector_type::vector_length);
      } else {
        r_val = -1;
      }
      return r_val;
    }
#endif

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::NoTranspose,Trans::Transpose,Algo::Gemm::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Unblocked>::
        invoke(C.extent(0), C.extent(1), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1(),
               B.data(), B.stride_1(), B.stride_0(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::NoTranspose,Trans::Transpose,Algo::Gemm::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Blocked>::
        invoke(C.extent(0), C.extent(1), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1(),
               B.data(), B.stride_1(), B.stride_0(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    ///
    /// T/T
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::Transpose,Trans::Transpose,Algo::Gemm::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      typedef typename CViewType::value_type vector_type;
      //typedef typename vector_type::value_type value_type;

      const int
        m = C.dimension(0),
        n = C.dimension(1),
        k = A.dimension(0);

      static_assert(is_vector<vector_type>::value, "value type is not vector type");      
      static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                    "AVX, AVX2 and AVX512 is supported");
      const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ?  MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

      // no error check
      int r_val = 0;
      if (A.stride_0() == 1 && B.stride_0() == 1 && C.stride_0() == 1) {
        mkl_dgemm_compact(MKL_COL_MAJOR, MKL_TRANS, MKL_TRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_1(), 
                            (const double*)B.data(), B.stride_1(), 
                            beta,
                            (double*)C.data(), C.stride_1(),
                            format, (MKL_INT)vector_type::vector_length);
      } else if (A.stride_1() == 1 && B.stride_1() == 1 && C.stride_1() == 1) {
        mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_TRANS, MKL_TRANS,
                            m, n, k, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (const double*)B.data(), B.stride_0(), 
                            beta,
                            (double*)C.data(), C.stride_0(),
                            format, (MKL_INT)vector_type::vector_length);
      } else {
        r_val = -1;
      }
      return r_val;
    }
#endif

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::Transpose,Trans::Transpose,Algo::Gemm::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Unblocked>::
        invoke(C.extent(0), C.extent(1), A.extent(0),
               alpha, 
               A.data(), A.stride_1(), A.stride_0(),
               B.data(), B.stride_1(), B.stride_0(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialGemm<Trans::Transpose,Trans::Transpose,Algo::Gemm::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      return SerialGemmInternal<Algo::Gemm::Blocked>::
        invoke(C.extent(0), C.extent(1), A.extent(0),
               alpha, 
               A.data(), A.stride_1(), A.stride_0(),
               B.data(), B.stride_1(), B.stride_0(),
               beta,
               C.data(), C.stride_0(), C.stride_1());
    }
        
  }
}


#endif
