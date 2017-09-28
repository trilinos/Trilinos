#ifndef __TACHOEXP_TRSM_HPP__
#define __TACHOEXP_TRSM_HPP__

/// \file TachoExp_Trsm.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {
  
  namespace Experimental {

    ///
    /// Trsm:
    ///

    /// various implementation for different uplo and algo parameters
    template<typename ArgSide,
             typename ArgUplo, 
             typename ArgTrans, 
             typename ArgAlgo>
    struct Trsm;

    /// task construction for the above chol implementation
    /// Trsm<ArgSide,ArgUplo,ArgTrans,ArgAlgo>::invoke(_sched, member, ArgDiag(), _alpha, _A, _B);
    template<typename SchedType,
             typename ScalarType,
             typename DenseMatrixViewType,
             typename ArgSide,
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct TaskFunctor_Trsm {
    public:
      typedef SchedType sched_type;
      typedef typename sched_type::member_type member_type;

      typedef ScalarType scalar_type;

      typedef DenseMatrixViewType dense_block_type;
      typedef typename dense_block_type::future_type future_type;
      typedef typename future_type::value_type value_type;
      
    private:
      sched_type _sched;
      scalar_type _alpha;
      dense_block_type _A, _B;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Trsm() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_Trsm(const sched_type &sched,
                       const scalar_type alpha,
                       const dense_block_type &A,
                       const dense_block_type &B) 
        : _sched(sched), 
          _alpha(alpha),
          _A(A), 
          _B(B) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Trsm<ArgSide,ArgUplo,ArgTrans,ArgAlgo>
          ::invoke(_sched, member, ArgDiag(), _alpha, _A, _B);

        if (member.team_rank() == 0) {
          _B.set_future();
          r_val = ierr;
        }
      }
    };
    
  }  
}

#endif
