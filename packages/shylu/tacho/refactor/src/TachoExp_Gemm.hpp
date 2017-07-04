#ifndef __TACHOEXP_GEMM_HPP__
#define __TACHOEXP_GEMM_HPP__

/// \file TachoExp_Herk.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    ///
    /// Gemm:
    ///

    /// various implementation for different uplo and algo parameters
    template<typename ArgTransA, typename ArgTransB, typename ArgAlgo>
    struct Gemm;

    /// task construction for the above chol implementation
    /// Gemm<ArgTransA,ArgTransB,ArgAlgo>::invoke(_sched, member, _alpha, _A, _B, _beta, _C);
    template<typename SchedType,
             typename ScalarType,
             typename DenseMatrixViewType,
             typename ArgTransA,
             typename ArgTransB,
             typename ArgAlgo>
    struct TaskFunctor_Gemm {
    public:
      typedef SchedType sched_type;
      typedef typename sched_type::member_type member_type;

      typedef ScalarType scalar_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;

    private:
      sched_type _sched;
      scalar_type _alpha, _beta;
      dense_block_type _A, _B, _C;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Gemm() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Gemm(const sched_type &sched,
                       const scalar_type alpha,
                       const dense_block_type &A,
                       const dense_block_type &B,
                       const scalar_type beta,
                       const dense_block_type &C)
        : _sched(sched),
          _alpha(alpha),
          _beta(beta),
          _A(A),
          _B(B),
          _C(C) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Gemm<ArgTransA,ArgTransB,ArgAlgo>
          ::invoke(_sched, member, _alpha, _A, _B, _beta, _C);

        if (member.team_rank() == 0) {
          _C.set_future();
          r_val = ierr;
        }
      }
    };

  }
}

#endif
