#pragma once
#ifndef __TEST_CRS_MATRIX_BASE_IO_HPP__
#define __TEST_CRS_MATRIX_BASE_IO_HPP__

/// \file test_crs_matrix_base_io.hpp
/// \brief Test CrsMatrixBase class file IO interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "crs_matrix_base.hpp"
#include "tmg_crs_matrix_base_simple.hpp"

namespace Example {
  
  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int testCrsMatrixBaseIO(const string filename) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef Tmg_CrsMatrixBase_Simple<CrsMatrixBaseType> TmgType;

    __DOT_LINE__;
    cout << "testCrsMatrixBaseIO:: filename = " << filename << endl;
    __DOT_LINE__;

    int r_val = 0, r_val_prev = 0;
    string eval;

    { // Read matrix
      r_val_prev = r_val;
      cout << "testCrsMatrixBaseIO::ReadFile::Begin - " << r_val << endl;

      ifstream in;
      in.open(filename);
      if (!in.good()) {
        cout << "Error in open the file: " << filename << endl;
        ++r_val;
      }

      CrsMatrixBaseType A("A, Imported from a file");
      A.importMatrixMarket(in);

      cout << A << endl;

      cout << "testCrsMatrixBaseIO::ReadFile::End - " << r_val << endl;

      __EVAL_STRING__(r_val - r_val_prev, eval);
      cout << "testCrsMatrixBaseIO::ReadFile::Eval - " << eval << endl;

    }
    __DOT_LINE__;
    
    return r_val;
  }
}

#endif
