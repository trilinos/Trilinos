#pragma once
#ifndef __TEST_CHOL_BY_BLOCKS_HPP__
#define __TEST_CHOL_BY_BLOCKS_HPP__

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "crs_matrix_helper.hpp"

#include "task_view.hpp"

#include "task_factory.hpp"

#include "chol.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testCholByBlocks(const string file_input) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    
    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    int r_val = 0;

    __DOT_LINE__;
    cout << "testCholByBlocks:: input = " << file_input << endl;        
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

    cout << "testCholByBlocks::Begin - " << r_val << endl;
    typename TaskFactoryType::policy_type policy;
    TaskFactoryType::setPolicy(&policy);

    {
      CrsTaskViewType U(&UU_Unblocked);
      U.fillRowViewArray();
    
      auto future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
                                                          ::TaskFunctor<CrsTaskViewType>(U), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());
    
      cout << UU_Unblocked << endl;
    }
    {
      CrsHierTaskViewType H(&HU);
      for (size_type k=0;k<HU.NumNonZeros();++k)
        HU.Value(k).fillRowViewArray();

      auto future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks>::
                                                          TaskFunctor<CrsHierTaskViewType>(H), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());

      cout << UU_ByBlocks << endl;
    }  
    
    {
      const auto epsilon = sqrt(NumericTraits<value_type>::epsilon());
      for (size_type k=0;k<UU_Unblocked.NumNonZeros();++k) {
        auto tmp = abs(UU_Unblocked.Value(k) - UU_ByBlocks.Value(k));
        __ASSERT_TRUE__(tmp < epsilon);
      }
    }
    cout << "testCholByBlocks::End - " << r_val << endl;  

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCholByBlocks::Eval - " << eval << endl;
    
    __DOT_LINE__;

    return r_val;
  }
}

#endif
