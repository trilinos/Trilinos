#pragma once
#ifndef __TEST_CRS_HIER_BASE_HPP__
#define __TEST_CRS_HIER_BASE_HPP__

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
  int testCrsHierBase(const string filename) {
    typedef double value_type;
    typedef int    ordinal_type;
    typedef int    size_type;
    
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    
    typedef CrsMatrixBase<CrsMatrixViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;    

    __DOT_LINE__;
    cout << "testCrsHierBase:: filename = " << filename << endl;
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
      CrsHierBaseType HU("HU");
      CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                 S.NumBlocks(), 
                                 S.RangeVector(),
                                 S.TreeVector());
      cout << HU << endl;
    }

    {
      CrsHierBaseType HL("HL");
      CrsMatrixHelper::flat2hier(Uplo::Lower, LL, HL, 
                                 S.NumBlocks(), 
                                 S.RangeVector(),
                                 S.TreeVector());
      cout << HL << endl;
    }

    {
      CrsHierBaseType HH("HH");
      CrsMatrixHelper::flat2hier(AA, HH);
      cout << HH << endl;
    }

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCrsHierBase::Eval - " << eval << endl;

    __DOT_LINE__;
   
    return r_val;
  }
}

#endif
