#pragma once
#ifndef __TEST_CHOL_BY_BLOCKS_GRAPHVIZ_HPP__
#define __TEST_CHOL_BY_BLOCKS_GRAPHVIZ_HPP__

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "crs_matrix_helper.hpp"

#include "task_view.hpp"

#include "task_policy_graphviz.hpp"

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
  int testCholByBlocksGraphviz(const string file_input,
                                const string file_output) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskFactory<TaskPolicy,Future> TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    
    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    int r_val = 0;

    __DOT_LINE__;
    cout << "testCholByBlocksGraphviz:: input = " << file_input 
         << ", output = " << file_output << endl;        
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

    CrsMatrixBaseType UU("UU");
    CrsHierMatrixBaseType HU("HU");
    {
      GraphHelperType S(AA);
      S.computeOrdering();

      CrsMatrixBaseType PA("Permuted AA");
      PA.copy(S.PermVector(), S.InvPermVector(), AA);

      UU.copy(Uplo::Upper, PA);

      CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 S.TreeVector());

      cout << S << endl;
      cout << HU << endl;
    }

    cout << "testCholByBlocksGraphviz::Begin - " << r_val << endl;
    typename TaskFactoryType::policy_type policy;
    TaskFactoryType::setPolicy(&policy);

    {
      CrsHierTaskViewType H(&HU);
      for (size_type k=0;k<HU.NumNonZeros();++k)
        HU.Value(k).fillRowViewArray();

      int r_val_chol = 0;

      Chol<Uplo::Upper,AlgoChol::ByBlocks>::
        TaskFunctor<CrsHierTaskViewType>(H).apply(r_val_chol);
    }  
    
    {
      ofstream out;
      out.open(file_output);
      if (!out.good()) {
        cout << "Error in open the file: " << file_output << endl;
        return -1;
      }
      
      TaskFactoryType::Policy().graphviz(out);
      TaskFactoryType::Policy().clear();
    }
    cout << "testCholByBlocksGraphviz::End - " << r_val << endl;  

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCholByBlocksGraphviz::Eval - " << eval << endl;
    
    __DOT_LINE__;

    return r_val;
  }
}

#endif
