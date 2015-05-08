#pragma once
#ifndef __TEST_CRS_TASK_VIEW_HPP__
#define __TEST_CRS_TASK_VIEW_HPP__

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"
#include "graph_helper_scotch.hpp" 

namespace Example {

  using namespace std;
  
  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testCrsTaskView(const string filename) {
    typedef double value_type;
    typedef int    ordinal_type;
    typedef int    size_type;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;
    
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
     
    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierTaskBaseType;

    typedef CrsMatrixView<CrsHierTaskBaseType> CrsHierViewType;
    typedef TaskView<CrsHierViewType,TaskFactoryType> CrsHierTaskViewType;

    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;    

    __DOT_LINE__;
    cout << "testCrsTaskView:: filename = " << filename << endl;
    __DOT_LINE__;

    int r_val = 0; //, r_val_prev = 0;
    
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

    CrsMatrixBaseType UU("UU"), LL("LL");
    GraphHelperType S(AA);
    {
      S.computeOrdering();

      CrsMatrixBaseType PA("Permuted AA");
      PA.copy(S.PermVector(), S.InvPermVector(), AA);

      UU.copy(Uplo::Upper, PA);
      LL.copy(Uplo::Lower, PA);
    }

    {
      CrsHierTaskBaseType HU("HU");
      CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                 S.NumBlocks(), 
                                 S.RangeVector(),
                                 S.TreeVector());
      cout << HU << endl;

      CrsHierTaskViewType H(&HU);
    }

    {
      CrsHierTaskBaseType HL("HL");
      CrsMatrixHelper::flat2hier(Uplo::Lower, LL, HL, 
                                 S.NumBlocks(), 
                                 S.RangeVector(),
                                 S.TreeVector());
      cout << HL << endl;

      CrsHierTaskViewType H(&HL);
    }

    {
      CrsHierTaskBaseType HH("HH");
      CrsMatrixHelper::flat2hier(AA, HH);
      cout << HH << endl;

      CrsHierTaskViewType H(&HH);
    }

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCrsTaskView::Eval - " << eval << endl;

    __DOT_LINE__;
   
    return r_val;
  }
}

#endif
