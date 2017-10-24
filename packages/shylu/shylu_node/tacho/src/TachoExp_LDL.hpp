#ifndef __TACHOEXP_LDL_HPP__
#define __TACHOEXP_LDL_HPP__

/// \file TachoExp_LDL.hpp
/// \brief Front interface for LDL^t dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {
  
  namespace Experimental {
    
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
    template<typename SchedType,
             typename DenseMatrixViewType,
             typename ArgUplo,
             typename ArgAlgo>
    struct TaskFunctor_LDL {
    public:
      typedef SchedType sched_type;
      typedef typename sched_type::member_type member_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;
      
    private:
      sched_type _sched;
      dense_block_type _A;
      dense_block_type _ipiv;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_LDL() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_LDL(const sched_type &sched,
                      const dense_block_type &A,
                      const dense_block_type &ipiv)
        : _sched(sched), 
          _A(A),
          _ipiv(ipiv) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = LDL<ArgUplo,ArgAlgo>
          ::invoke(_sched, member, _A, _ipiv);

        if (member.team_rank() == 0) {
          _A.set_future();
          r_val = ierr;
        }
      }
    };

  }
}

#endif
