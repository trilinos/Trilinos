#ifndef __KOKKOSBATCHED_TRSV_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_TRSV_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsv_Serial_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Serial Impl
    /// ===========

    ///
    /// Implemented:
    /// L/NT, U/NT, L/T, U/T
    /// 
    /// Not yet implemented
    /// L/CT, U/CT 

    ///
    /// L/NT
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        typedef typename bViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = b.extent(0),
          n = 1,
          vl = vector_type::vector_length;
        
        // no error check
        int r_val = 0;
        if (A.stride_0() == 1) { 
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasLeft, CblasLower, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1) {  
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasLeft, CblasLower, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif
    
    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalLower<Algo::Trsv::Unblocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(0), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalLower<Algo::Trsv::Blocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(0), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    ///
    /// L/T
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Lower,Trans::Transpose,ArgDiag,Algo::Trsv::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        typedef typename bViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = b.extent(0),
          n = 1,
          vl = vector_type::vector_length;
        
        // no error check
        int r_val = 0;
        if (A.stride_0() == 1) { 
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasLeft, CblasLower, CblasTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1) {  
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasLeft, CblasLower, CblasTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif
    
    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Lower,Trans::Transpose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(1), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Lower,Trans::Transpose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalUpper<Algo::Trsv::Blocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(1), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };

    ///
    /// U/NT
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        typedef typename bViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = b.extent(0),
          n = 1,
          vl = vector_type::vector_length;

        // no error check
        int r_val = 0;
        if (A.stride_0() == 1) { 
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasLeft, CblasUpper, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1) {  
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasLeft, CblasUpper, CblasNoTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif

    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(0), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalUpper<Algo::Trsv::Blocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(0), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    ///
    /// U/T
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Upper,Trans::Transpose,ArgDiag,Algo::Trsv::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        typedef typename bViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = b.extent(0),
          n = 1,
          vl = vector_type::vector_length;

        // no error check
        int r_val = 0;
        if (A.stride_0() == 1) { 
          cblas_dtrsm_compact(CblasColMajor, 
                              CblasLeft, CblasUpper, CblasTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else if (A.stride_1() == 1) {  
          cblas_dtrsm_compact(CblasRowMajor, 
                              CblasLeft, CblasUpper, CblasTrans, 
                              ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                              m, n, 
                              alpha, 
                              (const double*)A.data(), A.stride_0(), 
                              (double*)b.data(), b.stride_0(), 
                              (MKL_INT)vl, (MKL_INT)1);
        } else {
          r_val = -1;
        }
        return r_val;
      }
    };
#endif

    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Upper,Trans::Transpose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalLower<Algo::Trsv::Unblocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(1), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };

    template<typename ArgDiag>
    struct SerialTrsv<Uplo::Upper,Trans::Transpose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return SerialTrsvInternalLower<Algo::Trsv::Blocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.extent(1), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };
    
  }
}

#endif
