#pragma once
#ifndef __EXAMPLE_STAT_BY_BLOCKS_HPP__
#define __EXAMPLE_STAT_BY_BLOCKS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "symbolic_factor_helper.hpp"
#include "crs_matrix_helper.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleStatByBlocks(const string file_input,
                          const int treecut,
                          const int minblksize,
                          const int prunecut,
                          const int seed,
                          const int fill_level,
                          const int league_size,
                          const int histogram_size,
                          const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;
    typedef SymbolicFactorHelper<CrsMatrixBaseType> SymbolicFactorHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef CrsMatrixBase<CrsMatrixViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "StatByBlocks:: import input file = " << file_input << endl;
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
    cout << "StatByBlocks:: import input file::time = " << t << endl;


    CrsMatrixBaseType UU("UU");
    CrsHierMatrixBaseType HU("HU");
    {
      CrsMatrixBaseType PA("Permuted AA");
      typename GraphHelperType::size_type_array rptr(AA.Label()+"Graph::RowPtrArray", AA.NumRows() + 1);
      typename GraphHelperType::ordinal_type_array cidx(AA.Label()+"Graph::ColIndexArray", AA.NumNonZeros());

      AA.convertGraph(rptr, cidx);
      GraphHelperType S(AA.Label()+"ScotchHelper",
                        AA.NumRows(),
                        rptr,
                        cidx,
                        seed);
      {
        timer.reset();

        S.computeOrdering(treecut, minblksize);
        S.pruneTree(prunecut);

        PA.copy(S.PermVector(), S.InvPermVector(), AA);

        t = timer.seconds();

        if (verbose)
          cout << S << endl;
      }
      cout << "StatByBlocks:: reorder the matrix::time = " << t << endl;

      {
        SymbolicFactorHelperType F(PA, league_size);
        timer.reset();

        F.createNonZeroPattern(fill_level, Uplo::Upper, UU);

        t = timer.seconds();
        cout << "StatByBlocks:: AA (nnz) = " << AA.NumNonZeros() << ", UU (nnz) = " << UU.NumNonZeros() << endl;
      }
      cout << "StatByBlocks:: symbolic factorization::time = " << t << endl;

      {
        timer.reset();

        CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());

        for (ordinal_type k=0;k<HU.NumNonZeros();++k)
          HU.Value(k).fillRowViewArray();

        t = timer.seconds();
        cout << "StatByBlocks:: Hier (dof, nnz) = " << HU.NumRows() << ", " << HU.NumNonZeros() << endl;
      }
      cout << "StatByBlocks:: construct hierarchical matrix::time = " << t << endl;
    }

    {
      cout << endl;
      cout << " -- Flat matrix: UU --" << endl;
      cout << "    # of Rows     = " << UU.NumRows() << endl;
      cout << "    # of Cols     = " << UU.NumCols() << endl;
      cout << "    # of Nonzeros = " << UU.NumNonZeros() << endl;
      cout << endl;
      cout << " -- Hierarchical matrix: HU --" << endl;
      cout << "    # of Rows     = " << HU.NumRows() << endl;
      cout << "    # of Cols     = " << HU.NumCols() << endl;
      cout << "    # of Nonzeros = " << HU.NumNonZeros() << endl;
      cout << endl;
      cout << " -- Blocks of HU --" << endl;

      map<size_type,size_type> histogram;

      if (HU.NumNonZeros()) {

        size_type
          nnz_min = HU.Value(0).countNumNonZeros(),
          nnz_max = nnz_min,
          nnz_sum = 0,
          nnz_ave = 0;

        size_type nnz_cnt = 0;
        for (ordinal_type k=0;k<HU.NumNonZeros();++k) {
          const auto nnz_blk = HU.Value(k).countNumNonZeros();
          if (nnz_blk) {
            nnz_min  = min(nnz_min, nnz_blk);
            nnz_max  = max(nnz_max, nnz_blk);
            nnz_sum += nnz_blk;
            ++nnz_cnt;

            if (histogram_size)
              ++histogram[nnz_blk/histogram_size];
          }
        }
        nnz_ave = nnz_sum/nnz_cnt;

        cout << "    Min # of Nonzeros = " << nnz_min << endl;
        cout << "    Max # of Nonzeros = " << nnz_max << endl;
        cout << "    Ave # of Nonzeros = " << nnz_ave << endl;
        cout << "    Sum # of Nonzeros = " << nnz_sum << endl;
        cout << "    # of empty Blocks = " << (HU.NumNonZeros() - nnz_cnt) << endl;

        if (histogram_size) {
          cout << "    Histogram" << endl;
          for (auto it=histogram.begin();it!=histogram.end();++it)
            cout << (it->first*histogram_size) << " , " << it->second << endl;
        }
      } else {
        cout << "    No registered blocks" << endl;
      }

    }

    return r_val;
  }
}

#endif
