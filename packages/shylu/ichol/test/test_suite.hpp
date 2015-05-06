#pragma once
#ifndef __TEST_SUITE_HPP__
#define __TEST_SUITE_HPP__

/// \file test_suite.hpp
/// \brief Simple test suite
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "test_crs_matrix_base.hpp"
#include "test_crs_matrix_base_io.hpp"

namespace Example { 
  
  using namespace std;
  
  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  class TestSuite : public Disp {
  public:
    static string label;

    static int doUnitTests() {
      int r_val = 0;

      cout << label << "::doUnitTests::Begin" << endl;

      r_val += testCrsMatrixBase<ValueType,OrdinalType,SizeType,SpaceType,MemoryTraits>(0,0);
      r_val += testCrsMatrixBase<ValueType,OrdinalType,SizeType,SpaceType,MemoryTraits>(3,3);

      r_val += testCrsMatrixBaseIO<ValueType,OrdinalType,SizeType,SpaceType,MemoryTraits>("test_crs_input.mtx");

      cout << label << "::doUnitTests::End" << endl;

      string eval;
      __EVAL_STRING__(r_val, eval);

      __DOT_LINE__;
      cout << label << "::doUnitTests::Eval::" << eval << endl;
      return r_val;
    }
  };

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType,
           typename SpaceType,
           typename MemoryTraits> 
  string TestSuite<ValueType,
                   OrdinalType,
                   SizeType,
                   SpaceType,
                   MemoryTraits>::label = "NoName";
}

#endif
