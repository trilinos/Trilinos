#ifndef __TACHO_CHOL_HPP__
#define __TACHO_CHOL_HPP__

/// \file Tacho_Chol.hpp
/// \brief Front interface for Cholesky dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {
  
    ///
    /// Chol:
    /// 
    /// 

    /// various implementation for different uplo and algo parameters
    template<typename ArgUplo, 
             typename ArgAlgo>
    struct Chol;
    
    /// task construction for the above chol implementation
    // Chol<ArgUplo,ArgAlgo>::invoke(_sched, member, _A);
    template<typename SchedulerType,
             typename DenseMatrixViewType,
             typename ArgUplo,
             typename ArgAlgo>
    struct TaskFunctor_Chol {
    public:
      typedef SchedulerType scheduler_type;
      typedef typename scheduler_type::member_type member_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;
      
    private:
      scheduler_type _sched;
      dense_block_type _A;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Chol() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Chol(const scheduler_type &sched,
                       const dense_block_type &A)
        : _sched(sched), 
          _A(A) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Chol<ArgUplo,ArgAlgo>
          ::invoke(_sched, member, _A);

        Kokkos::single(Kokkos::PerTeam(member), [&]() {
            _A.set_future();
            r_val = ierr;
          });
      }
    };

}

#endif
