#pragma once
#ifndef __TEST_DENSE_TASK_VIEW_HPP__
#define __TEST_DENSE_TASK_VIEW_HPP__

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"
#include "graph_helper_scotch.hpp" 

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "dense_matrix_helper.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

namespace Example {

  using namespace std;
  
  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testDenseTaskView(const string filename,
                        const OrdinalType nb,
                        const OrdinalType nrhs) {
    typedef double value_type;
    typedef int    ordinal_type;
    typedef int    size_type;
    
    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;    

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    typedef DenseMatrixBase<DenseTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> DenseHierTaskBaseType;

    typedef DenseMatrixView<DenseHierTaskBaseType> DenseHierViewType;
    typedef TaskView<DenseHierViewType,TaskFactoryType> DenseHierTaskViewType;

    __DOT_LINE__;
    cout << "testDenseTaskView:: filename = " << filename << endl;
    __DOT_LINE__;

    int r_val = 0;
    
    CrsMatrixBaseType AA("AA");
    {
      ifstream in;
      in.open(filename);
      if (!in.good()) {
        cout << "Failed in open the file: " << filename << endl;
        return ++r_val;
      }
      AA.importMatrixMarket(in);
    }

    GraphHelperType S(AA);
    {
      S.computeOrdering();
    }

    DenseMatrixBaseType BB("BB", AA.NumRows(), nrhs);
    {
      DenseHierTaskBaseType HB("HB");
      DenseMatrixHelper::flat2hier(BB, HB,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   nb);
      cout << HB << endl;

      DenseHierTaskViewType H(&HB);
    }

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCrsTaskView::Eval - " << eval << endl;

    __DOT_LINE__;
   
    return r_val;
  }
}

#endif
