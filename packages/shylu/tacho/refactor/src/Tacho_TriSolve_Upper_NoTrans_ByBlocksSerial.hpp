#ifndef __TACHO_TRI_SOLVE_UPPER_NOTRANS_BY_BLOCKS_SERIAL_HPP__
#define __TACHO_TRI_SOLVE_UPPER_NOTRANS_BY_BLOCKS_SERIAL_HPP__

/// \file Tacho_TriSolve_Upper_NoTrans_ByBlocksSerial.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Tacho {

  template<int ArgVariant,
           template<int,int> class ControlType>
  class TriSolve<Uplo::Upper,Trans::NoTranspose,
                 AlgoTriSolve::ByBlocksSerial,ArgVariant,ControlType>{
  public:

    // data-parallel interface
    // =======================
    template<typename PolicyType,
             typename MemberType,
             typename CrsTaskViewTypeA,
             typename DenseTaskViewTypeB,
             typename DenseTaskViewTypeW>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const int diagA,
                      CrsTaskViewTypeA &A,
                      DenseTaskViewTypeB &B,
                      DenseTaskViewTypeW &W) {
      typedef typename CrsTaskViewTypeA::row_view_type row_view_type;
      
      // ---------------------------------------------
      if (member.team_rank() == 0) {
        CrsTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/           ABL, ABR,      A10, A11, A12,
          /**/                          A20, A21, A22;

        DenseTaskViewTypeB BT,      B0,
          /**/             BB,      B1,
          /**/                      B2;

        DenseTaskViewTypeW WT,      W0,
          /**/             WB,      W1,
          /**/                      W2;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 0, 0, Partition::BottomRight);

        Part_2x1(B,  BT,
                 /**/BB,
                 0, Partition::Bottom);

        Part_2x1(W,  WT,
                 /**/WB,
                 0, Partition::Bottom);

        while (ABR.NumRows() < A.NumRows()) {
          Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                          /*******/ /**/  A10, A11, A12,
                          ABL, ABR, /**/  A20, A21, A22,
                          1, 1, Partition::TopLeft);

          Part_2x1_to_3x1(BT,  /**/  B0,
                          /**/ /**/  B1,
                          BB,  /**/  B2,
                          1, Partition::Top);

          Part_2x1_to_3x1(WT,  /**/  W0,
                          /**/ /**/  W1,
                          WB,  /**/  W2,
                          1, Partition::Top);

          // -----------------------------------------------------

          // B1 = B1 - A12*B2;
          {
            row_view_type a(A12,0);
            const auto nnz = a.NumNonZeros();

            const auto jend = B1.NumCols();
            for (auto j=0;j<jend;++j) {
              const auto col_at_j = j;
              
              auto &cc = W1.Value(0, col_at_j);
              Impl::DenseMatrixTools::Serial::set(cc, 0.0);
              for (auto i=0;i<nnz;++i) {
                const auto row_at_i = a.Col(i);
                auto &aa = a.Value(i);
                auto &bb = B2.Value(row_at_i, col_at_j);                
                Gemm<Trans::NoTranspose,Trans::NoTranspose,
                  CtrlDetail(ControlType,AlgoTriSolve::ByBlocksSerial,ArgVariant,Gemm)>
                  ::invoke(policy, member, /**/ -1.0, aa, bb, 1.0, cc);
              }
              auto &dd = B1.Value(0, col_at_j);
              Impl::DenseMatrixTools::Serial::axpy(dd, cc, 1.0);
            }
          }

          // B1 = inv(triu(A11))*B1
          {
            row_view_type a(A11,0);
            auto &aa = a.Value(0);

            const auto jend = B1.NumCols();
            for (auto j=0;j<jend;++j) {
              auto &bb = B1.Value(0, j);
              Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
                CtrlDetail(ControlType,AlgoTriSolve::ByBlocksSerial,ArgVariant,Trsm)>
                ::invoke(policy, member, /**/ diagA, 1.0, aa, bb);
            }
          }
          
          // -----------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::BottomRight);
          
          Merge_3x1_to_2x1(B0, /**/   BT,
                           B1, /**/  /**/
                           B2, /**/   BB,
                           Partition::Bottom);
          
          Merge_3x1_to_2x1(W0, /**/   WT,
                           W1, /**/  /**/
                           W2, /**/   WB,
                           Partition::Bottom);
        }
      }
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeW>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      int _diagA;

      ExecViewTypeA _A;
      ExecViewTypeB _B;
      ExecViewTypeW _W;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const int diagA,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B,
                  const ExecViewTypeW &W)
        : _diagA(diagA),
          _A(A),
          _B(B),
          _W(W),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "TriSolve::ByBlocks::Serial"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = TriSolve::invoke(_policy, _policy.member_single(),
                                 _diagA, _A, _B, _W);
        _B.setFuture(typename ExecViewTypeB::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {

        if (member.team_rank() == 0) {
          _policy.clear_dependence(this);
          
          const int ierr = TriSolve::invoke(_policy, member,
                                            _diagA, _A, _B, _W);

          _B.setFuture(typename ExecViewTypeB::future_type());
        }
      }

    };

    template<typename PolicyType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeW>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeW>
    createTaskFunctor(const PolicyType &policy,
                      const int diagA,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B,
                      const ExecViewTypeW &W) {
      return TaskFunctor<PolicyType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeW>
        (policy, diagA, A, B, W);
    }

  };
}
#endif
