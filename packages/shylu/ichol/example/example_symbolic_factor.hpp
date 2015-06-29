#pragma once
#ifndef __EXAMPLE_SYMBOLIC_FACTOR_HPP__
#define __EXAMPLE_SYMBOLIC_FACTOR_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"

#include "symbolic_factor_helper.hpp"

namespace Example {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleSymbolicFactor(const string file_input,
                            const int fill_level,
                            const int league_size,
                            const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef SymbolicFactorHelper<CrsMatrixBaseType> SymbolicFactorHelperType;
    
    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "SymbolicFactor:: import input file = " << file_input << endl;        
    CrsMatrixBaseType AA("AA");
    {
      timer.reset();

      ifstream in;
      in.open(file_input);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_input << endl;
        return ++r_val;
      }
      AA.importMatrixMarket(in);
      t = timer.seconds();
      
      if (verbose)
        cout << AA << endl;
    }
    cout << "SymbolicFactor:: import input file::time = " << t << endl;        

    CrsMatrixBaseType UU("UU");    
    {
      timer.reset();

      SymbolicFactorHelperType symbolic(AA, league_size);
      symbolic.createNonZeroPattern(fill_level, Uplo::Upper, UU);

      t = timer.seconds();

      if (verbose) {
        cout << symbolic << endl;
        cout << UU << endl;
      }
    }
    cout << "SymbolicFactor:: factorize the matrix::time = " << t << endl; 
    
    return r_val;
  }
}

#endif
