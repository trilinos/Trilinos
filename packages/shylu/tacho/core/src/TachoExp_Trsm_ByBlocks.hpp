#ifndef __TACHOEXP_TRSM_BYBLOCKS_HPP__
#define __TACHOEXP_TRSM_BYBLOCKS_HPP__


/// \file  TachoExp_Trsm_ByBlocks.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Partition.hpp"

#include "TachoExp_Trsm.hpp"
#include "TachoExp_Trsm_External.hpp"

#include "TachoExp_Gemm.hpp"
#include "TachoExp_Gemm_External.hpp"
#include "TachoExp_Gemm_ByBlocks.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<>
    struct Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks> {
      template<typename SchedType,
               typename MemberType,
               typename DiagType,
               typename ScalarType,
               typename MatrixOfDenseBlocksType>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const DiagType diagA,
             const ScalarType alpha,
             const MatrixOfDenseBlocksType &A,
             const MatrixOfDenseBlocksType &B) {
        typedef SchedType sched_type;
        typedef ScalarType scalar_type;
        typedef typename MatrixOfDenseBlocksType::value_type dense_block_type;
        typedef typename dense_block_type::future_type future_type;

        if (get_team_rank(member) == 0) {
          MatrixOfDenseBlocksType ATL, ATR,      A00, A01, A02,
            /**/                  ABL, ABR,      A10, A11, A12,
            /**/                                 A20, A21, A22;

          MatrixOfDenseBlocksType BT,            B0,
            /**/                  BB,            B1,
            /**/                                 B2;
        
          Part_2x2(A,  ATL, ATR,
                   /**/ABL, ABR,
                   0, 0, Partition::TopLeft);

          Part_2x1(B,  BT,
                   /**/BB,
                   0, Partition::Top);
        
          while (ATL.dimension_0() < A.dimension_0()) {
            const auto alpha_select = (ATL.dimension_0() > 0 ? ScalarType(1.0) : alpha);

            Part_2x2_to_3x3(ATL, ATR, /**/ A00, A01, A02,
                            /*******/ /**/ A10, A11, A12,
                            ABL, ABR, /**/ A20, A21, A22,
                            1, 1, Partition::BottomRight);
          
            Part_2x1_to_3x1(BT,  /**/ B0,
                            /**/ /**/ B1,
                            BB,  /**/ B2,
                            1, Partition::Bottom);
          
            //------------------------------------------------------------
            auto &aa = A11(0, 0);
            const auto n = B1.dimension_1();
            for (auto j=0;j<n;++j) {
              auto &bb = B1(0, j);

              const future_type dep[2] = { aa.future(), bb.future() };            
              future_type f =
                Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 2), Kokkos::TaskPriority::High),
                                   TaskFunctor_Trsm
                                   <sched_type,scalar_type,dense_block_type,
                                   Side::Left,Uplo::Upper,Trans::ConjTranspose,DiagType,Algo::External>
                                   (sched, alpha_select, aa, bb));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
              bb.set_future(f);
            }

            Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::ByBlocks>
              ::invoke(sched, member, -1.0, A12, B1, alpha_select, B2);
          
            //------------------------------------------------------------
            Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                             A10, A11, A12, /**/ /******/
                             A20, A21, A22, /**/ ABL, ABR,
                             Partition::TopLeft);
            
            Merge_3x1_to_2x1(B0, /**/ BT,
                             B1, /**/ /**/
                             B2, /**/ BB,
                             Partition::Top);
          }
        }
        return 0;
      }
    };

    template<>
    struct Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::ByBlocks> {
      template<typename SchedType,
               typename MemberType,
               typename DiagType,
               typename ScalarType,
               typename MatrixOfDenseBlocksType>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const DiagType diagA,
             const ScalarType alpha,
             const MatrixOfDenseBlocksType &A,
             const MatrixOfDenseBlocksType &B) {
        typedef SchedType sched_type;
        typedef ScalarType scalar_type;
        typedef typename MatrixOfDenseBlocksType::value_type dense_block_type;
        typedef typename dense_block_type::future_type future_type;

        if (get_team_rank(member) == 0) {
          MatrixOfDenseBlocksType ATL, ATR,      A00, A01, A02,
            /**/                  ABL, ABR,      A10, A11, A12,
            /**/                                 A20, A21, A22;

          MatrixOfDenseBlocksType BT,            B0,
            /**/                  BB,            B1,
            /**/                                 B2;
        
          Part_2x2(A,  ATL, ATR,
                   /**/ABL, ABR,
                   0, 0, Partition::BottomRight);

          Part_2x1(B,  BT,
                   /**/BB,
                   0, Partition::Bottom);
        
          while (ABR.dimension_0() < A.dimension_0()) {
            const auto alpha_select = (ABR.dimension_0() > 0 ? ScalarType(1.0) : alpha);

            Part_2x2_to_3x3(ATL, ATR, /**/ A00, A01, A02,
                            /*******/ /**/ A10, A11, A12,
                            ABL, ABR, /**/ A20, A21, A22,
                            1, 1, Partition::TopLeft);
          
            Part_2x1_to_3x1(BT,  /**/ B0,
                            /**/ /**/ B1,
                            BB,  /**/ B2,
                            1, Partition::Top);
          
            //------------------------------------------------------------
            Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::ByBlocks>
              ::invoke(sched, member, -1.0, A12, B2, alpha_select, B1);
          
            auto &aa = A11(0, 0);
            const auto n = B1.dimension_1();
            for (auto j=0;j<n;++j) {
              auto &bb = B1(0, j);

              const future_type dep[2] = { aa.future(), bb.future() };            
              future_type f =
                Kokkos::task_spawn(Kokkos::TaskTeam(sched, Kokkos::when_all(dep, 2), Kokkos::TaskPriority::High),
                                   TaskFunctor_Trsm
                                   <sched_type,scalar_type,dense_block_type,
                                   Side::Left,Uplo::Upper,Trans::NoTranspose,DiagType,Algo::External>
                                   (sched, alpha_select, aa, bb));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task_spawn return a null future");
              bb.set_future(f);
            }

            //------------------------------------------------------------
            Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                             A10, A11, A12, /**/ /******/
                             A20, A21, A22, /**/ ABL, ABR,
                             Partition::BottomRight);
            
            Merge_3x1_to_2x1(B0, /**/ BT,
                             B1, /**/ /**/
                             B2, /**/ BB,
                             Partition::Bottom);
          }
        }
        return 0;
      }
    };



  }
}
#endif
