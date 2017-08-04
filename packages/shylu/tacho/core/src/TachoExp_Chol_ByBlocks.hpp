#ifndef __TACHOEXP_CHOL_BYBLOCKS_HPP__
#define __TACHOEXP_CHOL_BYBLOCKS_HPP__

/// \file  TachoExp_Chol_ByBlocks.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Partition.hpp"

#include "TachoExp_Chol.hpp"
#include "TachoExp_Chol_External.hpp"

#include "TachoExp_Trsm.hpp"
#include "TachoExp_Trsm_External.hpp"

#include "TachoExp_Herk.hpp"
#include "TachoExp_Herk_External.hpp"

#include "TachoExp_Gemm.hpp"
#include "TachoExp_Gemm_External.hpp"

namespace Tacho {

  namespace Experimental {

    /// LAPACK Chol
    /// ===========
    template<>
    struct Chol<Uplo::Upper,Algo::ByBlocks> {
      template<typename SchedType,
               typename MemberType,
               typename MatrixOfDenseBlocksType>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const MatrixOfDenseBlocksType &A) {
        typedef SchedType sched_type;

        typedef typename MatrixOfDenseBlocksType::value_type dense_block_type;
        typedef typename dense_block_type::value_type value_type;
        typedef typename dense_block_type::future_type future_type;

        typedef typename TypeTraits<value_type>::magnitude_type scalar_type;
        
        int r_val = 0;      

        if (get_team_rank(member) == 0) {
          MatrixOfDenseBlocksType ATL, ATR,      A00, A01, A02,
            /**/                  ABL, ABR,      A10, A11, A12,
            /**/                                 A20, A21, A22;
          
          Part_2x2(A,  ATL, ATR,
                   /**/ABL, ABR,
                   0, 0, Partition::TopLeft);
          while (ATL.dimension_0() < A.dimension_0()) {
            Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                            /*******/ /**/  A10, A11, A12,
                            ABL, ABR, /**/  A20, A21, A22,
                            1, 1, Partition::BottomRight);
            // -----------------------------------------------------
            
            // A11 = chol(A11)
            {
              auto &aa = A11(0,0);
              future_type f = 
                Kokkos::task_spawn(Kokkos::TaskTeam(sched, aa.future(), Kokkos::TaskPriority::High),
                                   TaskFunctor_Chol
                                   <sched_type,dense_block_type,
                                   Uplo::Upper,Algo::External>(sched, aa));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
              aa.set_future(f);
            }
            
            // A12 = inv(triu(A11)') * A12
            {
              auto &aa = A11(0, 0); 
              const ordinal_type n = A12.dimension_1();
              for (auto j=0;j<n;++j) {
                auto &bb = A12(0, j); 

                const future_type dep[2] = { aa.future(), bb.future() };
                future_type f = 
                  Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 2), Kokkos::TaskPriority::High),
                                     TaskFunctor_Trsm
                                     <sched_type,scalar_type,dense_block_type,
                                     Side::Left,Uplo::Upper,Trans::ConjTranspose,Diag::NonUnit,Algo::External>
                                     (sched, 1.0, aa, bb));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
                bb.set_future(f);
              }
            }

            // A22 = A22 - A12' * A12
            {
              const ordinal_type n = A22.dimension_1();
              for (auto j=0;j<n;++j) {
                {
                  auto &aa = A12(0, j);
                  auto &cc = A22(j, j);

                  const future_type dep[] = { aa.future(), cc.future() };
                  future_type f = 
                    Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 2), Kokkos::TaskPriority::High),
                                       TaskFunctor_Herk
                                       <sched_type,scalar_type,dense_block_type,
                                       Uplo::Upper,Trans::ConjTranspose,Algo::External>(sched, -1.0, aa, 1.0, cc));
                  TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
                  cc.set_future(f);
                }
                for (auto i=0;i<j;++i) {
                  auto &aa = A12(0, i);
                  auto &bb = A12(0, j);
                  auto &cc = A22(i, j);

                  const future_type dep[] = { aa.future(), bb.future(), cc.future() };
                  future_type f = 
                    Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 3), Kokkos::TaskPriority::High),
                                       TaskFunctor_Gemm
                                       <sched_type,scalar_type,dense_block_type,
                                       Trans::ConjTranspose,Trans::NoTranspose,Algo::External>(sched, -1.0, aa, bb, 1.0, cc));
                  TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
                  cc.set_future(f);
                }
              }
            }
            // -----------------------------------------------------
            Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                             A10, A11, A12, /**/ /******/
                             A20, A21, A22, /**/ ABL, ABR,
                             Partition::TopLeft);
          }
        }
        
        return r_val;
      }
    };
  }
}

#endif
