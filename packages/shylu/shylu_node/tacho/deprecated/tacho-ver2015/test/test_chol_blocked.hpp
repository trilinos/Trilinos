#pragma once
#ifndef __TEST_CHOL_BLOCKED_HPP__
#define __TEST_CHOL_BLOCKED_HPP__

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

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
  int testCholBlocked(const string file_input,
                       const OrdinalType blocksize,
                       const string file_factored) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;
  
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    
    int r_val = 0;

    __DOT_LINE__;
    cout << "testCholBlocked:: input = " << file_input 
         << ", blocksize = " << blocksize 
         << ", factored = " << file_factored 
         << endl;        
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

    cout << "testCholBlocked::Begin - " << r_val << endl;
    typename TaskFactoryType::policy_type policy;
    TaskFactoryType::setPolicy(&policy);
    
    CrsTaskViewType U(&UU);
    {
      Chol<Uplo::Upper,AlgoChol::Blocked>::blocksize = blocksize;
      
      Chol<Uplo::Upper,AlgoChol::Blocked>
        ::invoke(TaskFactoryType::Policy(), 
                          TaskFactoryType::Policy().member_single(),
                          U);
      
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
      for (size_type k=0;k<UU.NumNonZeros();++k) {
        auto tmp = abs(UU.Value(k) - FF.Value(k));
        __ASSERT_TRUE__(tmp < epsilon);
      }
    }
    cout << "testCholBlocked::End - " << r_val << endl;    

    string eval;
    __EVAL_STRING__(r_val, eval);
    cout << "testCholBlocked::Eval - " << eval << endl;
    
    __DOT_LINE__;

    return r_val;
  }
}

#endif
