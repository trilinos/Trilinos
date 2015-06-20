#pragma once
#ifndef __TEST_DENSE_HIER_BASE_HPP__
#define __TEST_DENSE_HIER_BASE_HPP__

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"
#include "graph_helper_scotch.hpp" 

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "dense_matrix_helper.hpp"

namespace Example {

  using namespace std;
  
  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testDenseHierBase(const string filename,
                        const OrdinalType nb,
                        const OrdinalType nrhs) {
    typedef double value_type;
    typedef int    ordinal_type;
    typedef int    size_type;
    
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;    

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;
    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;

    typedef DenseMatrixBase<DenseMatrixViewType,ordinal_type,size_type,SpaceType,MemoryTraits> DenseHierBaseType;

    __DOT_LINE__;
    cout << "testDenseHierBase:: filename = " << filename << endl;
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
      DenseHierBaseType HB("HB");
      DenseMatrixHelper::flat2hier(BB, HB,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   nb);
      cout << HB << endl;
    }

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCrsHierBase::Eval - " << eval << endl;

    __DOT_LINE__;
   
    return r_val;
  }
}

#endif
