#ifndef __TACHO_LDL_HPP__
#define __TACHO_LDL_HPP__

/// \file Tacho_LDL.hpp
/// \brief Front interface for LDL^t dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {
  
    ///
    /// LDL:
    /// 
    /// 

    /// various implementation for different uplo and algo parameters
    template<typename ArgUplo, 
             typename ArgAlgo>
    struct LDL;
    
    /// task construction for the above chol implementation
    // LDL<ArgUplo,ArgAlgo>::invoke(_sched, member, _A);
    template<typename SchedulerType,
             typename DenseMatrixViewType,
             typename ArgUplo,
             typename ArgAlgo>
    struct TaskFunctor_LDL {
    public:
      typedef SchedulerType scheduler_type;
      typedef typename scheduler_type::member_type member_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;
      
    private:
      scheduler_type _sched;
      dense_block_type _A;
      dense_block_type _ipiv;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_LDL() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_LDL(const scheduler_type &sched,
                      const dense_block_type &A,
                      const dense_block_type &ipiv)
        : _sched(sched), 
          _A(A),
          _ipiv(ipiv) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = LDL<ArgUplo,ArgAlgo>
          ::invoke(_sched, member, _A, _ipiv);

        Kokkos::single(Kokkos::PerTeam(member), [&] () {
            _A.set_future();
            r_val = ierr;
          });
      }
    };

}

#endif
