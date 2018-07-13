#pragma once
#ifndef __EXAMPLE_CHOL_PERFORMANCE_DEVICE_HPP__
#define __EXAMPLE_CHOL_PERFORMANCE_DEVICE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLU_NodeTacho_config.h"

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "symbolic_factor_helper.hpp"
#include "crs_matrix_helper.hpp"

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
  int exampleCholPerformanceDevice(const string file_input,
                                   const int treecut,
                                   const int prunecut,
                                   const int seed,
                                   const int nthreads,
                                   const int max_task_dependence,
                                   const int max_concurrency,
                                   const int team_size,
                                   const int fill_level,
                                   const int league_size,
                                   const bool skip_serial,
                                   const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;
    typedef SymbolicFactorHelper<CrsMatrixBaseType> SymbolicFactorHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double
      t_import = 0.0,
      t_reorder = 0.0,
      t_symbolic = 0.0,
      t_flat2hier = 0.0,
      t_factor_seq = 0.0,
      t_factor_task = 0.0;

    cout << "CholPerformanceDevice:: import input file = " << file_input << endl;
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

      t_import = timer.seconds();

      if (verbose)
        cout << AA << endl;
    }
    cout << "CholPerformanceDevice:: import input file::time = " << t_import << endl;

    cout << "CholPerformanceDevice:: reorder the matrix" << endl;
    CrsMatrixBaseType PA("Permuted AA");
    CrsMatrixBaseType UU("UU");     // permuted base upper triangular matrix
    CrsHierMatrixBaseType HU("HU"); // hierarchical matrix of views

    {
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

        S.computeOrdering(treecut, 0);
        S.pruneTree(prunecut);

        PA.copy(S.PermVector(), S.InvPermVector(), AA);

        t_reorder = timer.seconds();

        if (verbose)
          cout << S << endl
               << PA << endl;
      }

      cout << "CholPerformanceDevice:: reorder the matrix::time = " << t_reorder << endl;
      {
        SymbolicFactorHelperType F(PA, league_size);
        timer.reset();
        F.createNonZeroPattern(fill_level, Uplo::Upper, UU);
        t_symbolic = timer.seconds();
        cout << "CholPerformanceDevice:: AA (nnz) = " << AA.NumNonZeros() << ", UU (nnz) = " << UU.NumNonZeros() << endl;

        if (verbose)
          cout << F << endl
               << UU << endl;
      }
      cout << "CholPerformanceDevice:: symbolic factorization::time = " << t_symbolic << endl;
      {
        timer.reset();
        CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());

        for (ordinal_type k=0;k<HU.NumNonZeros();++k)
          HU.Value(k).fillRowViewArray();

        t_flat2hier = timer.seconds();

        cout << "CholPerformanceDevice:: Hier (dof, nnz) = " << HU.NumRows() << ", " << HU.NumNonZeros() << endl;
      }
      cout << "CholPerformanceDevice:: construct hierarchical matrix::time = " << t_flat2hier << endl;
    }

    // copy of UU
    CrsMatrixBaseType RR("RR");
    RR.copy(UU);

    cout << "CholPerformanceDevice:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "CholPerformanceDevice:: max task size   = " << max_task_size << endl;

    if (!skip_serial) {
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   1);

      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType U(&UU);
      U.fillRowViewArray();

      cout << "CholPerformanceDevice:: Serial factorize the matrix" << endl;
      {
        UU.copy(RR);
        timer.reset();
        {
          Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
            ::invoke(TaskFactoryType::Policy(),
                              TaskFactoryType::Policy().member_single(),
                              U);
        }
        t_factor_seq = timer.seconds();
        if (verbose)
          cout << UU << endl;
      }
      cout << "CholPerformanceDevice:: Serial factorize the matrix::time = " << t_factor_seq << endl;
    }

    {
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   team_size);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      cout << "CholPerformanceDevice:: ByBlocks factorize the matrix:: team_size = " << team_size << endl;
      CrsHierTaskViewType H(&HU);
      {
        UU.copy(RR);
        timer.reset();
        {
          auto future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks>::
                                                              TaskFunctor<CrsHierTaskViewType>(H), 0);
          TaskFactoryType::Policy().spawn(future);
          Kokkos::Experimental::wait(TaskFactoryType::Policy());
        }
        t_factor_task += timer.seconds();

        if (verbose)
          cout << UU << endl;
      }
      cout << "CholPerformanceDevice:: ByBlocks factorize the matrix::time = " << t_factor_task << endl;
    }

    if (!skip_serial) {
      cout << "CholPerformanceDevice:: task scale [seq/task] = " << t_factor_seq/t_factor_task << endl;
    }

    return r_val;
  }
}

#endif
