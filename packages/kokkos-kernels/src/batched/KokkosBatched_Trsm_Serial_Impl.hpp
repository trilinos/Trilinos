#ifndef __KOKKOSBATCHED_TRSM_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_TRSM_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsm_Serial_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// L/L/NT
    ///
    /// B := inv(tril(A)) (alpha*B)
    /// A(m x m), B(m x n)

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)    
    template<typename ArgDiag>
    struct SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsm::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        typedef typename BViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = B.dimension(0),
          n = B.dimension(1),
          vl = vector_type::vector_length;

        // no error check
        int r_val = 0;
        if (A.stride_0() == 1 && B.stride_0() == 1) {
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasLeft, CblasLower, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)B.data(), B.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1 && B.stride_1() == 1) {
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasLeft, CblasLower, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)B.data(), B.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif    

    template<typename ArgDiag>
    struct SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsm::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        return SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(ArgDiag::use_unit_diag,
                                                                    B.extent(0), B.extent(1),
                                                                    alpha, 
                                                                    A.data(), A.stride_0(), A.stride_1(),
                                                                    B.data(), B.stride_0(), B.stride_1());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsm::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        return SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(ArgDiag::use_unit_diag,
                                                                  B.extent(0), B.extent(1),
                                                                  alpha, 
                                                                  A.data(), A.stride_0(), A.stride_1(),
                                                                  B.data(), B.stride_0(), B.stride_1());
      }
    };

    ///
    /// R/U/NT
    ///
    /// B := (alpha*B) inv(triu(A))
    /// A(n x n), B(m x n)
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<typename ArgDiag>
    struct SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        typedef typename BViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = B.dimension(0),
          n = B.dimension(1),
          vl = vector_type::vector_length;

        // no error check
        int r_val = 0;
        if (A.stride_0() == 1 && B.stride_0() == 1) {
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasRight, CblasUpper, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)B.data(), B.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1 && B.stride_1() == 1) {
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasRight, CblasUpper, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)B.data(), B.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif

    template<typename ArgDiag>
    struct SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        return SerialTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(ArgDiag::use_unit_diag,
                                                                          B.extent(1), B.extent(0),
                                                                          alpha, 
                                                                          A.data(), A.stride_1(), A.stride_0(),
                                                                          B.data(), B.stride_1(), B.stride_0());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        return SerialTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(ArgDiag::use_unit_diag,
                                                                        B.extent(1), B.extent(0),
                                                                        alpha, 
                                                                        A.data(), A.stride_1(), A.stride_0(),
                                                                        B.data(), B.stride_1(), B.stride_0());
      }      
    };
    
    ///
    /// L/U/NT
    ///
    /// B := inv(triu(A)) (alpha*B) 
    /// A(m x m), B(m x n)
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<typename ArgDiag>
    struct SerialTrsm<Side::Left,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        typedef typename BViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = B.dimension(0),
          n = B.dimension(1),
          vl = vector_type::vector_length;

        // no error check
        int r_val = 0;
        if (A.stride_0() == 1 && B.stride_0() == 1) {
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasLeft, CblasUpper, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)B.data(), B.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1 && B.stride_1() == 1) {
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasLeft, CblasUpper, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)B.data(), B.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif

    template<typename ArgDiag>
    struct SerialTrsm<Side::Left,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        return SerialTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(ArgDiag::use_unit_diag,
                                                                    B.extent(0), B.extent(1),
                                                                    alpha, 
                                                                    A.data(), A.stride_0(), A.stride_1(),
                                                                    B.data(), B.stride_0(), B.stride_1());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsm<Side::Left,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        return SerialTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(ArgDiag::use_unit_diag,
                                                                  B.extent(0), B.extent(1),
                                                                  alpha, 
                                                                  A.data(), A.stride_0(), A.stride_1(),
                                                                  B.data(), B.stride_0(), B.stride_1());
      }
    };

  }
}

#endif
