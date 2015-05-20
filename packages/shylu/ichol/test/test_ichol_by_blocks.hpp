#pragma once
#ifndef __TEST_ICHOL_BY_BLOCKS_HPP__
#define __TEST_ICHOL_BY_BLOCKS_HPP__

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "crs_matrix_helper.hpp"

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
  int testICholByBlocks(const string file_input) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType>,
      Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;
    typedef ParallelFor ForType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    
    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    int r_val = 0;

    __DOT_LINE__;
    cout << "testICholByBlocks:: input = " << file_input << endl;        
    __DOT_LINE__;

    CrsMatrixBaseType AA("AA");
    {
      ifstream in;
      in.open(file_input);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_input << endl;
        return ++r_val;
      }
      AA.importMatrixMarket(in);
    }


    CrsMatrixBaseType UU_Unblocked("UU_Unblocked"), UU_ByBlocks("UU_ByBlocks");
    CrsHierMatrixBaseType HU("HU");
    {
      GraphHelperType S(AA);
      S.computeOrdering();

      CrsMatrixBaseType PA("Permuted AA");
      PA.copy(S.PermVector(), S.InvPermVector(), AA);

      UU_Unblocked.copy(Uplo::Upper, PA);
      UU_ByBlocks.copy(Uplo::Upper, PA);

      CrsMatrixHelper::flat2hier(Uplo::Upper, UU_ByBlocks, HU,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 S.TreeVector());
    }

    cout << "testICholByBlocks::Begin - " << r_val << endl;
    {
      CrsTaskViewType U(&UU_Unblocked);
      U.fillRowViewArray();
    
      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
                                       ::TaskFunctor<ForType,CrsTaskViewType>(U), 0);
      policy.spawn(future);
      Kokkos::Experimental::wait(policy);
    
      cout << UU_Unblocked << endl;
    }
    {
      CrsHierTaskViewType H(&HU);
      for (ordinal_type k=0;k<HU.NumNonZeros();++k)
        HU.Value(k).fillRowViewArray();

      typedef typename CrsTaskViewType::policy_type policy_type;
      policy_type policy;
      auto future = policy.create_team(IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
                                       TaskFunctor<ForType,CrsHierTaskViewType>(H), 0);
      policy.spawn(future);
      Kokkos::Experimental::wait(policy);

      cout << UU_ByBlocks << endl;
    }  
    
    {
      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      for (ordinal_type k=0;k<UU_Unblocked.NumNonZeros();++k) {
        auto tmp = abs(UU_Unblocked.Value(k) - UU_ByBlocks.Value(k));
        __ASSERT_TRUE__(tmp < epsilon);
      }
    }
    cout << "testICholByBlocks::End - " << r_val << endl;  

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testICholByBlocks::Eval - " << eval << endl;
    
    __DOT_LINE__;

    return r_val;
  }
}

#endif
