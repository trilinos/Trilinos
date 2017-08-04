#ifndef __TACHOEXP_HERK_BYBLOCKS_HPP__
#define __TACHOEXP_HERK_BYBLOCKS_HPP__


/// \file  Tacho_Herk_ByBlocks.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Partition.hpp"

#include "TachoExp_Herk.hpp"
#include "TachoExp_Herk_External.hpp"

#include "TachoExp_Gemm.hpp"
#include "TachoExp_Gemm_External.hpp"

namespace Tacho {

  namespace Experimental {

    template<>
    struct Herk<Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks> {
      template<typename SchedType,
               typename MemberType,
               typename ScalarType,
               typename MatrixOfDenseBlocksType>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const ScalarType alpha,
             const MatrixOfDenseBlocksType &A,
             const ScalarType beta,
             const MatrixOfDenseBlocksType &C) {
        typedef SchedType sched_type;
        typedef ScalarType scalar_type;
        typedef typename MatrixOfDenseBlocksType::value_type dense_block_type;
        //typedef typename dense_block_type::value_type value_type;
        typedef typename dense_block_type::future_type future_type;

        if (get_team_rank(member) == 0) {
          const ordinal_type pend = A.dimension_0();
          for (ordinal_type p=0;p<pend;++p) {
            const auto beta_select = (p > 0 ? ScalarType(1.0) : beta);
            const ordinal_type k2end = C.dimension_1();
            for (ordinal_type k2=0;k2<k2end;++k2) {
              for (ordinal_type k1=0;k1<(k2+1);++k1) {
                auto &aa = A(p,  k1);
                auto &cc = C(k1, k2);
                if (k1 == k2) {
                  const future_type dep[2] = { aa.future(), cc.future() };
                  future_type f =
                    Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 2), Kokkos::TaskPriority::High),
                                       TaskFunctor_Herk
                                       <sched_type,scalar_type,dense_block_type,
                                       Uplo::Upper,Trans::ConjTranspose,Algo::External>
                                       (sched, alpha, aa, beta_select, cc));
                  TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
                  cc.set_future(f);
                } else {
                  auto &bb = A(p, k2);
                  const future_type dep[3] = { aa.future(), bb.future(), cc.future() };
                  future_type f =
                    Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 3), Kokkos::TaskPriority::High),
                                       TaskFunctor_Gemm
                                       <sched_type,scalar_type,dense_block_type,
                                       Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
                                       (sched, alpha, aa, bb, beta_select, cc));
                  TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
                  cc.set_future(f);
                }
              }
            }
          }
        }
        return 0;
      }
    };
  }
}
#endif
