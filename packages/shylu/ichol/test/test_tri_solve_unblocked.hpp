#pragma once
#ifndef __TEST_TRI_SOLVE_UNBLOCKED_HPP__
#define __TEST_TRI_SOLVE_UNBLOCKED_HPP__

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "parallel_for.hpp"

#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "tri_solve.hpp"

#include "tmg_dense_matrix_base_simple.hpp"

namespace Example {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testTriSolveUnblocked(const string file_input, 
                            const OrdinalType nrhs) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;
    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    
    typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType>,
      Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;

    typedef ParallelFor ForType;
    
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    typedef Tmg_DenseMatrixBase_Simple<DenseMatrixBaseType> TmgType;
    
    int r_val = 0;

    __DOT_LINE__;
    cout << "testTriSolveUnblocked:: input = " << file_input << ", nrhs = " << nrhs << endl;        
    __DOT_LINE__;

    CrsMatrixBaseType AA("AA"), UU("UU");    
    {
      ifstream in;
      in.open(file_input);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_input << endl;
        return ++r_val;
      }
      AA.importMatrixMarket(in);

      UU.copy(Uplo::Upper, AA);
    }

    TmgType tmg(UU.NumRows(), nrhs);

    DenseMatrixBaseType 
      BB("BB", UU.NumRows(), nrhs), 
      CC("CC", UU.NumRows(), nrhs);

    cout << "testTriSolveUnblocked::Begin - " << r_val << endl;
    CrsTaskViewType U(&UU);
    DenseTaskViewType B(&BB), C(&CC);
    
    cout << UU << endl;
    U.fillRowViewArray();
    {
      r_val += tmg.fill(BB);

      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::Unblocked>
                                       ::TaskFunctor<ForType,CrsTaskViewType,DenseTaskViewType>
                                       (Diag::NonUnit, U, B), 0);
      
      policy.spawn(future);
      Kokkos::Experimental::wait(policy);

      cout << BB << endl;
    }
    {
      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(Gemm<Trans::ConjTranspose,Trans::NoTranspose,AlgoGemm::ForTriSolveBlocked>
                                       ::TaskFunctor<ForType,double,
                                       CrsTaskViewType,DenseTaskViewType,DenseTaskViewType>
                                       (1.0, U, B, 0.0, C), 0);
      
      policy.spawn(future);
      Kokkos::Experimental::wait(policy);

      r_val += tmg.check(CC);
    }
    {
      r_val += tmg.fill(BB);

      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::Unblocked>
                                       ::TaskFunctor<ForType,CrsTaskViewType,DenseTaskViewType>
                                       (Diag::NonUnit, U, B), 0);

      policy.spawn(future);
      Kokkos::Experimental::wait(policy);
      
      cout << BB << endl;
    }
    {
      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(Gemm<Trans::NoTranspose,Trans::NoTranspose,AlgoGemm::ForTriSolveBlocked>
                                       ::TaskFunctor<ForType,double,
                                       CrsTaskViewType,DenseTaskViewType,DenseTaskViewType>
                                       (1.0, U, B, 0.0, C), 0);
      
      policy.spawn(future);
      Kokkos::Experimental::wait(policy);

      r_val += tmg.check(CC);
    }
    cout << "testTriSolveUnblocked::End - " << r_val << endl;    

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testTriSolveUnblocked::Eval - " << eval << endl;
    
    __DOT_LINE__;

    return r_val;
  }
}

#endif
