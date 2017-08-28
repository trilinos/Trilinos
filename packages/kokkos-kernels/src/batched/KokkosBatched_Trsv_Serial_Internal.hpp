#ifndef __KOKKOSBATCHED_TRSV_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRSV_SERIAL_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"
#include "KokkosBatched_Gemv_Team_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Serial Internal Impl
    /// ====================

    ///
    /// Lower
    ///

    template<typename AlgoType>
    struct SerialTrsvInternalLower {
      template<typename ScalarType,
               typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const bool use_unit_diag,
             const int m,
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             /**/  ValueType *__restrict__ b, const int bs0);
    };

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialTrsvInternalLower<Algo::Trsv::Unblocked>::
    invoke(const bool use_unit_diag,
           const int m,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ b, const int bs0) {
      typedef ValueType value_type;

      if (alpha == 0)   SerialSetInternal::invoke(m, 0, b, bs0);
      else {
        if (alpha != 1) SerialScaleInternal::invoke(m, value_type(alpha), b, bs0);
        if (m <= 0) return 0;

        for (int p=0;p<m;++p) {
          const int iend = m-p-1;

          const value_type
            *__restrict__ a21   = iend ? A+(p+1)*as0+p*as1 : NULL;

          value_type
            *__restrict__ beta1 =        b+p*bs0,
            *__restrict__ b2    = iend ? beta1+bs0 : NULL;

          if (!use_unit_diag)
            *beta1 /= A[p*as0+p*as1];

          for (int i=0;i<iend;++i)
            b2[i*bs0] -= a21[i*as0] * (*beta1);
        }
      }
      return 0;
    }

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialTrsvInternalLower<Algo::Trsv::Blocked>::
    invoke(const bool use_unit_diag,
           const int m,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ b, const int bs0) {
      typedef ValueType value_type;
      enum : int {
        mbAlgo = Algo::Trsv::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };

      if (alpha == 0)   SerialSetInternal::invoke(m, 0, b, bs0);
      else {
        if (alpha != 1) SerialScaleInternal::invoke(m, value_type(alpha), b, bs0);
        if (m <= 0) return 0;

        /// case cuda: team size is large and blocksize (mb,nb) is small
        InnerTrsmLeftLowerUnitDiag<mbAlgo>    trsm_u(as0, as1, 1, 1);
        InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_n(as0, as1, 1, 1);

        const int mb = mbAlgo;
        for (int p=0;p<m;p+=mb) {
          const int pb = ((p+mb) > m ? (m-p) : mb);

          // trsm update
          const value_type *__restrict__ Ap = A+p*as0+p*as1;
          /**/  value_type *__restrict__ bp = b+p*bs0;

          if (use_unit_diag) trsm_u.serial_invoke(Ap, pb, 1, bp);
          else               trsm_n.serial_invoke(Ap, pb, 1, bp);

          // gemv update
          GemvInternal<Algo::Gemv::Blocked>
            ::invoke(m-p-pb, pb,
                     -1,
                     Ap+pb*as0, as0, as1,
                     bp, 1,
                     1,
                     bp+pb*bs0, bs0);
        }
      }
      return 0;
    }

    ///
    /// Upper
    ///

    template<typename AlgoType>
    struct SerialTrsvInternalUpper {
      template<typename ScalarType,
               typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const bool use_unit_diag,
             const int m,
             const ScalarType alpha,
             const ValueType *__restrict__ A, const int as0, const int as1,
             /**/  ValueType *__restrict__ b, const int bs0);

    };

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialTrsvInternalUpper<Algo::Trsv::Unblocked>::
    invoke(const bool use_unit_diag,
           const int m,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ b, const int bs0) {
      typedef ValueType value_type;

      if (alpha == 0)   SerialSetInternal::invoke(m, 0, b, bs0);
      else {
        if (alpha != 1) SerialScaleInternal::invoke(m, value_type(alpha), b, bs0);
        if (m <= 0) return 0;

        value_type *__restrict__ b0 = b;
        for (int p=(m-1);p>=0;--p) {
          const int iend = p;

          const value_type *__restrict__ a01   = A+p*as1;
          /**/  value_type *__restrict__ beta1 = b+p*bs0;

          if (!use_unit_diag)
            *beta1 /= A[p*as0+p*as1];

          for (int i=0;i<iend;++i)
            b0[i*bs0] -= a01[i*as0] * (*beta1);
        }
      }
      return 0;
    }

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialTrsvInternalUpper<Algo::Trsv::Blocked>::
    invoke(const bool use_unit_diag,
           const int m,
           const ScalarType alpha,
           const ValueType *__restrict__ A, const int as0, const int as1,
           /**/  ValueType *__restrict__ b, const int bs0) {
      typedef ValueType value_type;
      enum : int {
        mbAlgo = Algo::Trsm::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };

      // note that parallel range is different ( m*n vs m-1*n);
      if (alpha == 0)   SerialSetInternal::invoke(m, 0, b, bs0);
      else {
        if (alpha != 1) SerialScaleInternal::invoke(m, value_type(alpha), b, bs0);
        if (m <= 0) return 0;

        InnerTrsmLeftUpperUnitDiag<mbAlgo>    trsm_u(as0, as1, 1, 1);
        InnerTrsmLeftUpperNonUnitDiag<mbAlgo> trsm_n(as0, as1, 1, 1);

        const int mb = mbAlgo;
        for (int pp=0;pp<m;pp+=mb) {
          const int
            ptmp = (m - pp - mb),
            p = (ptmp < 0 ? 0 : ptmp),
            pb = (mb + (ptmp < 0)*ptmp);

          // trsm update
          const value_type *__restrict__ Ap = A+p*as0+p*as1;
          /**/  value_type *__restrict__ bp = b+p*bs0;

          if (use_unit_diag) trsm_u.serial_invoke(Ap, pb, 1, bp);
          else               trsm_n.serial_invoke(Ap, pb, 1, bp);

          // gemv update
          GemvInternal<Algo::Gemv::Blocked>
            ::invoke(p, pb,
                     -1,
                     Ap-p*as0, as0, as1,
                     bp, 1,
                     1,
                     b, bs0);
        }
      }
      return 0;
    }

  }
}

#endif
