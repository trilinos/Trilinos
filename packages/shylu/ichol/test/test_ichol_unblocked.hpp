#pragma once
#ifndef __TEST_ICHOL_UNBLOCKED_HPP__
#define __TEST_ICHOL_UNBLOCKED_HPP__

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "parallel_for.hpp"

#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "ichol.hpp"

namespace Example {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testICholUnblocked(const string file_input, 
                         const string file_factored) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

    typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType>,
      Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;

    typedef ParallelFor ForType;
    
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    
    int r_val = 0;

    __DOT_LINE__;
    cout << "testICholUnblocked:: input = " << file_input << ", factored = " << file_factored << endl;        
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

    cout << "testICholUnblocked::Begin - " << r_val << endl;
    CrsTaskViewType U(&UU);
    {
      U.fillRowViewArray();
    
      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
                                       ::TaskFunctor<ForType,CrsTaskViewType>(U), 0);
      policy.spawn(future);
      Kokkos::Experimental::wait(policy);
    
      cout << UU << endl;
    }   

    CrsMatrixBaseType FF("FF");    
    {
      ifstream in;
      in.open(file_factored);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_factored << endl;
        return ++r_val;
      }
      FF.importMatrixMarket(in);
      cout << FF << endl;
    }

    {
      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      for (ordinal_type k=0;k<UU.NumNonZeros();++k) {
        auto tmp = abs(UU.Value(k) - FF.Value(k));
        __ASSERT_TRUE__(tmp < epsilon);
      }
    }
    cout << "testICholUnblocked::End - " << r_val << endl; 

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testICholUnblocked::Eval - " << eval << endl;
    
    __DOT_LINE__;

    return r_val;
  }
}

#endif
