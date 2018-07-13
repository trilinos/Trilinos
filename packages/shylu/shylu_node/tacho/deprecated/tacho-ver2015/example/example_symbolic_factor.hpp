#pragma once
#ifndef __EXAMPLE_SYMBOLIC_FACTOR_HPP__
#define __EXAMPLE_SYMBOLIC_FACTOR_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "graph_helper_scotch.hpp"
#include "graph_helper_camd.hpp"
#include "crs_matrix_base.hpp"

#include "symbolic_factor_helper.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleSymbolicFactor(const string file_input,
                            const int treecut,
                            const int minblksize,
                            const int seed,
                            const int fill_level,
                            const int league_size,
                            const bool scotch,
                            const bool camd,
                            const bool symbolic,
                            const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType_Scotch;
    typedef GraphHelper_CAMD<CrsMatrixBaseType> GraphHelperType_CAMD;
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

      cout << "SymbolicFactor:: AA nnz = " << AA.NumNonZeros() << endl;
      if (verbose)
        cout << AA << endl;
    }
    cout << "SymbolicFactor:: import input file::time = " << t << endl;

    CrsMatrixBaseType PA("AA after Scotch");

    typename GraphHelperType_Scotch::size_type_array rptr(AA.Label()+"Graph::RowPtrArray", AA.NumRows() + 1);      
    typename GraphHelperType_Scotch::ordinal_type_array cidx(AA.Label()+"Graph::ColIndexArray", AA.NumNonZeros()); 
                                                                                                              
    AA.convertGraph(rptr, cidx);                                                                            
    GraphHelperType_Scotch S(AA.Label()+"ScotchHelper",                                                            
                             AA.NumRows(),                                                                         
                             rptr,                                                                                 
                             cidx,
                             seed);  

    S.setStratGraph(SCOTCH_STRATLEVELMAX   |
                    SCOTCH_STRATLEVELMIN   |
                    SCOTCH_STRATLEAFSIMPLE |
                    SCOTCH_STRATSEPASIMPLE);
    ordinal_type level = int(log2(double(AA.NumRows()))) - treecut;
    cout << "SymbolicFactor:: Scotch level = " << level<< endl;
    S.setTreeLevel(max(level, 1));

    if (scotch) {
      timer.reset();
      S.computeOrdering(treecut, minblksize);
      PA.copy(S.PermVector(), S.InvPermVector(), AA);
      t = timer.seconds();

      if (verbose)
        cout << S << endl
             << PA << endl;

      cout << "SymbolicFactor:: Scotch reorder the matrix::time = " << t << endl;
    } else {
      PA = AA;
      t = 0.0;
    }


    CrsMatrixBaseType CA("AA after CAMD");
    GraphHelperType_CAMD C("AfterScotch",
                           S.NumRows(),
                           S.NumNonZeros(),
                           S.RowPtrVector(),
                           S.ColIndexVector(),
                           S.NumBlocks(), 
                           S.RangeVector());
    if (scotch && camd) {
      timer.reset();
      C.computeOrdering();
      CA.copy(C.PermVector(), C.InvPermVector(), PA);
      t = timer.seconds();

      if (verbose)
        cout << C << endl
             << CA << endl;

      cout << "SymbolicFactor:: CAMD reorder the matrix::time = " << t << endl;
    } else {
      CA = PA;
      t = 0.0;
    }


    CrsMatrixBaseType UU("UU");
    if (symbolic) {
      timer.reset();
      SymbolicFactorHelperType symbolic(CA, league_size);
      symbolic.createNonZeroPattern(fill_level, Uplo::Upper, UU);
      t = timer.seconds();
      cout << "SymbolicFactor:: UU nnz = " << UU.NumNonZeros() << endl;

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
