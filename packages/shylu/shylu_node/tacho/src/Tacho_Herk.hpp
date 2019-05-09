#ifndef __TACHO_HERK_HPP__
#define __TACHO_HERK_HPP__

/// \file Tacho_Herk.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

    ///
    /// Herk:
    /// 
    
    /// various implementation for different uplo and algo parameters    
    template<typename ArgUplo, typename ArgTrans, typename ArgAlgo>
    struct Herk;
    
    /// task construction for the above chol implementation
    /// Herk<ArgUplo,ArgTrans,ArgAlgo> ::invoke(_sched, member, _alpha, _A, _beta, _C);
    template<typename SchedulerType,
             typename ScalarType,
             typename DenseMatrixViewType,
             typename ArgUplo,
             typename ArgTrans,
             typename ArgAlgo>
    struct TaskFunctor_Herk {
    public:
      typedef SchedulerType scheduler_type;
      typedef typename scheduler_type::member_type member_type;

      typedef ScalarType scalar_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;

    private:
      scheduler_type _sched;
      scalar_type _alpha, _beta;
      dense_block_type _A, _C;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Herk() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Herk(const scheduler_type &sched,
                       const scalar_type alpha,
                       const dense_block_type &A,
                       const scalar_type beta,
                       const dense_block_type &C)
        : _sched(sched),
          _alpha(alpha),
          _beta(beta),
          _A(A),
          _C(C) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Herk<ArgUplo,ArgTrans,ArgAlgo>
          ::invoke(_sched, member, _alpha, _A, _beta, _C);

        Kokkos::single(Kokkos::PerThread(member), [&] () {
            _C.set_future();
            r_val = ierr;
          });
      }
    };
    
}

#endif
