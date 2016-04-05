#ifndef __TACHO_TRI_SOLVE_UPPER_NOTRANS_BY_BLOCKS_HPP__
#define __TACHO_TRI_SOLVE_UPPER_NOTRANS_BY_BLOCKS_HPP__

/// \file Tacho_TriSolve_Upper_NoTrans_ByBlocks.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Tacho {

  template<int ArgVariant,
           template<int,int> class ControlType = Control>
  class TriSolve<Uplo::Upper,Trans::NoTrans,
                 AlgoTriSolve::ByBlocks,ArgVariant,ControlType>{
  public:

    // data-parallel interface
    // =======================
    template<typename PolicyType,
             typename MemberType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const int diagA,
                      CrsTaskViewTypeA &A,
                      DenseTaskViewTypeB &B) {

      typedef typename CrsTaskViewTypeA::ordinal_type      ordinal_type;

      typedef typename CrsTaskViewTypeA::value_type        crs_value_type;
      typedef typename CrsTaskViewTypeA::row_view_type     row_view_type;

      typedef typename DenseTaskViewTypeB::value_type      dense_value_type;

      typedef typename CrsTaskViewTypeA::future_type       future_type;

      TaskFactory factory;

      // ---------------------------------------------
      if (member.team_rank() == 0) {
        const ordinal_type ntasks_window = 4096;
        ordinal_type ntasks_spawned = 0;

        CrsTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/           ABL, ABR,      A10, A11, A12,
          /**/                          A20, A21, A22;

        DenseTaskViewTypeB BT,      B0,
          /**/             BB,      B1,
          /**/                      B2;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 0, 0, Partition::BottomRight);

        Part_2x1(B,  BT,
                 /**/BB,
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

          // -----------------------------------------------------

          // B1 = B1 - A12*B2;
          genGemmTasks_TriSolveUpperNoTransposeByBlocks(policy, A12, B2, B1);
          row_view_type a(A12,0);
          const ordinal_type nnz = a.NumNonZeros();

          for (ordinal_type i=0;i<nnz;++i) {
            const ordinal_type row_at_i = a.Col(i);
            crs_value_type &aa = a.Value(i);

            for (ordinal_type j=0;j<B1.NumCols();++j) {
              const ordinal_type col_at_j = j;
              dense_value_type &bb = B2.Value(row_at_i, col_at_j);
              dense_value_type &cc = B1.Value(0, col_at_j);

              const future_type f = factory.create<future_type>
                (policy,
                 Gemm<Trans::NoTranspose,Trans::NoTranspose,
                 CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Gemm)>
                 ::createTaskFunctor(policy, -1.0, aa, bb, 1.0, cc), 3);

              // dependence
              factory.depend(policy, f, aa.Future());
              factory.depend(policy, f, bb.Future());
              factory.depend(policy, f, cc.Future());

              // place task signature on y
              cc.setFuture(f);

              // spawn a task
              factory.spawn(policy, f);
              ++ntasks_spawned;
            }
          }

          // B1 = inv(triu(A11))*B1
          genTrsmTasks_TriSolveUpperNoTransposeByBlocks(policy, diagA, A11, B1);
          row_view_type a(A11,0);
          crs_value_type &aa = a.Value(0);

          for (ordinal_type j=0;j<B1.NumCols();++j) {
            dense_value_type &bb = B1.Value(0, j);

            const future_type f = factory.create<future_type>
              (policy,
               Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
               CtrlDetail(ControlType,AlgoTriSolve::ByBlocks,ArgVariant,Trsm)>
               ::createTaskFunctor(policy, diagA, 1.0, aa, bb), 2);
            
            // trsm dependence
            factory.depend(policy, f, aa.Future());
            factory.depend(policy, f, bb.Future());

            // place task signature on b
            bb.setFuture(f);

            // spawn a task
            factory.spawn(policy, f);
            ++ntasks_spawned;
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
          
          if (ntasks_spawned > ntasks_window)
            break;
        }
        
        A = ATL;
        B = BT;
      }
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      int _diagA;
      ScalarType _alpha;

      ExecViewTypeA _A;
      ExecViewTypeB _B;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const int diagA,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B)
        : _diagA(diagA),
          _alpha(alpha),
          _A(A),
          _B(B),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Trsm"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Trsm::invoke(_policy, _policy.member_single(),
                             _diagA, _alpha, _A, _B);
        _B.setFuture(typename ExecViewTypeB::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {

        if (member.team_rank() == 0) {
          _policy.clear_dependence(this);

          const int ierr = Trsm::invoke(_policy, member,
                                        _diagA, _alpha, _A, _B);

          if (_A.NumRows()) {
            _policy.respawn_needing_memory(this);
          } else {
            _B.setFuture(typename ExecViewTypeB::future_type());
          }
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB>
    createTaskFunctor(const PolicyType &policy,
                      const int diagA,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB>
        (policy, diagA, alpha, A, B);
    }

  }


#endif
