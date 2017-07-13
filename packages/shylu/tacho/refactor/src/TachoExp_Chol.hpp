#ifndef __TACHOEXP_CHOL_HPP__
#define __TACHOEXP_CHOL_HPP__

/// \file TachoExp_Chol.hpp
/// \brief Front interface for Cholesky dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {
  
  namespace Experimental {
    
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
    template<typename SchedType,
             typename DenseMatrixViewType,
             typename ArgUplo,
             typename ArgAlgo>
    struct TaskFunctor_Chol {
    public:
      typedef SchedType sched_type;
      typedef typename sched_type::member_type member_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;
      
    private:
      sched_type _sched;
      dense_block_type _A;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Chol() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Chol(const sched_type &sched,
                       const dense_block_type &A)
        : _sched(sched), 
          _A(A) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Chol<ArgUplo,ArgAlgo>
          ::invoke(_sched, member, _A);

        if (member.team_rank() == 0) {
          _A.set_future();
          r_val = ierr;
        }
      }
    };

  }
}

#endif
