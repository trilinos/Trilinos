#ifndef __TACHO_CHOL_UPPER_BY_BLOCKS_HPP__
#define __TACHO_CHOL_UPPER_BY_BLOCKS_HPP__

/// \file Tacho_Chol_Upper_ByBlocks.hpp
/// \brief Cholesky factorization by-blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 
  
  // Chol-By-Blocks
  // =================================================================================
  template<int ArgVariant, template<int,int> class ControlType>
  class Chol<Uplo::Upper,
             AlgoChol::ByBlocks,ArgVariant,ControlType> {
  public:

    // function interface
    // ==================
    template<typename PolicyType,
             typename MemberType,
             typename CrsTaskViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      MemberType &member,
                      CrsTaskViewTypeA &A,
                      unsigned int &part) {

      typedef typename CrsTaskViewTypeA::row_view_type row_view_type;
#ifdef TACHO_EXECUTE_TASKS_SERIAL
#else
      typedef typename CrsTaskViewTypeA::future_type future_type;
#endif

      if (member.team_rank() == 0) {
        const unsigned int ntasks_window  = TaskWindow::CholByBlocks;
        /**/  unsigned int ntasks_spawned = 0;

        CrsTaskViewTypeA ATL, ATR,      A00, A01, A02,
          /**/           ABL, ABR,      A10, A11, A12,
          /**/                          A20, A21, A22;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 part, part, Partition::TopLeft);

        while (ATL.NumRows() < A.NumRows()) {
          Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                          /*******/ /**/  A10, A11, A12,
                          ABL, ABR, /**/  A20, A21, A22,
                          1, 1, Partition::BottomRight);
          // -----------------------------------------------------

          row_view_type diag(A11, 0); 
          const bool diagNull = diag.Value(0).isNull();

          if (!diagNull) {
            // A11 = chol(A11)
            {
              auto &aa = diag.Value(0);
              
#ifdef TACHO_EXECUTE_TASKS_SERIAL
              Chol<Uplo::Upper,
                CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Chol)>
                ::invoke(policy, member, aa);
#else
              switch (ArgVariant) {
              case Variant::Three: // SuperNodes-ByBlocks with DenseByBlocks
                {
                  Chol<Uplo::Upper,
                    CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Chol)>
                    ::invoke(policy, member, aa);
                  break;
                }
              case Variant::One:   // SuperNodes-ByBlocks
              case Variant::Two:   // Sparse-ByBlocks
                { 
                  const auto task_type     = Kokkos::TaskTeam;
                  const auto task_priority = Kokkos::TaskHighPriority;
                  
                  // construct a task
                  const auto f = 
                    policy.task_spawn(Chol<Uplo::Upper,
                                      CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Chol)>
                                      ::createTaskFunctor(policy, aa), 
                                      aa.Future(),
                                      task_type, task_priority);
                  TACHO_TEST_FOR_ABORT(f.is_null(),
                                       ">> Tacho::CholByBlocks(Upper) returns a null future (out of memory)");
                  aa.setFuture(f);
                  break;
                }
              }
              ++ntasks_spawned;
#endif
            }
            
            // A12 = inv(triu(A11)') * A12
            {
              row_view_type b(A12, 0); 
              auto &aa = diag.Value(0);
              
              const auto nnz = b.NumNonZeros();
              for (auto j=0;j<nnz;++j) {
                auto &bb = b.Value(j);
                
                if (!bb.isNull()) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
                    CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Trsm)>
                    ::invoke(policy, member, Diag::NonUnit, 1.0, aa, bb);
#else
                  switch (ArgVariant) {
                  case Variant::Three:  // SuperNodes-ByBlocks with DenseByBlocks
                    {
                      Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
                        CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Trsm)>
                        ::invoke(policy, member, Diag::NonUnit, 1.0, aa, bb);
                      break;
                    }
                  case Variant::One:    // Sparse-ByBlocks
                  case Variant::Two:    // SuperNodes-ByBlocks
                    {
                      const auto task_type     = Kokkos::TaskTeam;
                      const auto task_priority = Kokkos::TaskRegularPriority;
                      
                      const future_type dep[] = { aa.Future(), bb.Future() };
                      const auto f = 
                        policy.task_spawn(Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
                                          CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Trsm)>
                                          ::createTaskFunctor(policy, Diag::NonUnit, 1.0, aa, bb),
                                          policy.when_all(2,dep),
                                          task_type, task_priority);
                      TACHO_TEST_FOR_ABORT(f.is_null(),
                                           ">> Tacho::DenseGemmByBlocks(NoTrans,NoTrans) returns a null future (out of memory)");
                      bb.setFuture(f);
                      break;
                    }
                  }
                  ++ntasks_spawned;
#endif
                } 
              }
            }
            
            // A22 = A22 - A12' * A12
            {
              // case that X.transpose, A.no_transpose, Y.no_transpose
              row_view_type a(A12,0);
              
              const auto nnz = a.NumNonZeros();
              
              // update herk
              for (auto i=0;i<nnz;++i) {
                const auto row_at_i = a.Col(i);
                auto &aa = a.Value(i);
                
                if (!aa.isNull()) {
                  row_view_type c(A22, row_at_i);
                  
                  auto idx = 0;
                  for (auto j=i;j<nnz && (idx > -2);++j) {
                    const auto col_at_j = a.Col(j);
                    auto &bb = a.Value(j);
                    
                    if (!bb.isNull()) {
                      if (row_at_i == col_at_j) {
                        idx = c.Index(row_at_i, idx);
                        if (idx >= 0) {
                          auto &cc = c.Value(idx);
                          
                          if (!cc.isNull()) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                            Herk<Uplo::Upper,Trans::ConjTranspose,
                              CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Herk)>
                              ::invoke(policy, member, -1.0, aa, 1.0, cc);
#else
                            switch (ArgVariant) {
                            case Variant::Three:
                              {
                                Herk<Uplo::Upper,Trans::ConjTranspose,
                                  CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Herk)>
                                  ::invoke(policy, member, -1.0, aa, 1.0, cc);
                                break;
                              }
                            case Variant::One:
                            case Variant::Two:
                              {
                                const auto task_type     = Kokkos::TaskTeam;
                                const auto task_priority = Kokkos::TaskHighPriority;
                                
                                const future_type dep[] = { aa.Future(), cc.Future() };
                                const auto f = 
                                  policy.task_spawn(Herk<Uplo::Upper,Trans::ConjTranspose,
                                                    CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Herk)>
                                                    ::createTaskFunctor(policy, -1.0, aa, 1.0, cc), 
                                                    policy.when_all(2,dep),
                                                    task_type, task_priority);
                                TACHO_TEST_FOR_ABORT(f.is_null(),
                                                     ">> Tacho::DenseGemmByBlocks(NoTrans,NoTrans) returns a null future (out of memory)");
                                cc.setFuture(f);
                                break;
                              }
                            }
                            ++ntasks_spawned;
#endif
                          } 
                        }
                      } else {
                        idx = c.Index(col_at_j, idx);
                        if (idx >= 0) {
                          auto &cc = c.Value(idx);
                          
                          if (!cc.isNull()) {
#ifdef TACHO_EXECUTE_TASKS_SERIAL
                            Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                              CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Gemm)>
                              ::invoke(policy, member, -1.0, aa, bb, 1.0, cc);
#else
                            switch (ArgVariant) {
                            case Variant::Three: 
                              {
                                Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                                  CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Gemm)>
                                  ::invoke(policy, member, -1.0, aa, bb, 1.0, cc);
                                break;
                              }
                            case Variant::One:
                            case Variant::Two: 
                              {
                                const auto task_type     = Kokkos::TaskTeam;
                                const auto task_priority = Kokkos::TaskRegularPriority;

                                const future_type dep[] = { aa.Future(), bb.Future(), cc.Future() };
                                const auto f = 
                                  policy.task_spawn(Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                                                    CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Gemm)>
                                                    ::createTaskFunctor(policy, -1.0, aa, bb, 1.0, cc), 
                                                    policy.when_all(3,dep),
                                                    task_type, task_priority);
                                TACHO_TEST_FOR_ABORT(f.is_null(),
                                                     ">> Tacho::DenseGemmByBlocks(NoTrans,NoTrans) returns a null future (out of memory)");                                
                                cc.setFuture(f);
                                break;
                              }
                            }
                            ++ntasks_spawned;
#endif
                          } 
                        }
                      }
                    } 
                  }
                }
              }
            }
          }
          
          // -----------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::TopLeft);

          if (ntasks_spawned > ntasks_window) 
            break;
        }

        part = ATL.NumRows();
      }
      
      return 0;
    }
    
    // task-data parallel interface
    // ============================
    template<typename PolicyType,
             typename ExecViewTypeA>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;
      
    private:
      ExecViewTypeA _A;
      
      PolicyType _policy;

      unsigned int _part;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() {}
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ExecViewTypeA A)
        : _A(A),
          _policy(policy),
          _part(0)
      { } 
      
      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Chol::ByBlocks"; }
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {

        // return for only team leader
        if (member.team_rank() == 0) {
          //_policy.clear_dependence(this);

          const int ierr = Chol::invoke(_policy, member,
                                        _A, _part);

          if (_part < _A.NumRows()) {
            _policy.respawn(this, Kokkos::TaskLowPriority);
            //_policy.respawn_needing_memory(this); 
          } else {
            _A.setFuture(typename ExecViewTypeA::future_type());
          }
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ExecViewTypeA>
    createTaskFunctor(const PolicyType &policy,
                      const ExecViewTypeA &A) {
      return TaskFunctor<PolicyType,ExecViewTypeA>
        (policy, A);
    }

  };
}

#endif
