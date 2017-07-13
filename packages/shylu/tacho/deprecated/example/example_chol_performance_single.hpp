#pragma once
#ifndef __EXAMPLE_CHOL_PERFORMANCE_SINGLE_HPP__
#define __EXAMPLE_CHOL_PERFORMANCE_SINGLE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

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

#ifdef HAVE_SHYLUTACHO_VTUNE
#include "ittnotify.h"
#endif

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleCholPerformanceSingle(const string file_input,
                                   const int treecut,
                                   const int minblksize,
                                   const int prunecut,
                                   const int seed,
                                   const int nthreads,
                                   const int max_task_dependence,
                                   const int team_size,
                                   const int fill_level,
                                   const int league_size,
                                   const bool team_interface,
                                   const bool skip_serial,
                                   const bool vtune_symbolic,
                                   const bool vtune_serial,
                                   const bool vtune_task,
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

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double
      t_import = 0.0,
      t_reorder = 0.0,
      t_symbolic = 0.0,
      t_flat2hier = 0.0,
      t_factor_seq = 0.0,
      t_factor_task = 0.0;

    cout << "CholPerformance:: import input file = " << file_input << endl;
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
    cout << "CholPerformance:: import input file::time = " << t_import << endl;

    cout << "CholPerformance:: reorder the matrix" << endl;
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

        S.computeOrdering(treecut, minblksize);
        S.pruneTree(prunecut);

        PA.copy(S.PermVector(), S.InvPermVector(), AA);

        t_reorder = timer.seconds();

        if (verbose)
          cout << S << endl
               << PA << endl;
      }

      cout << "CholPerformance:: reorder the matrix::time = " << t_reorder << endl;
      {
        SymbolicFactorHelperType F(PA, league_size);
        if (vtune_symbolic) {
#ifdef HAVE_SHYLUTACHO_VTUNE
          __itt_resume();
#endif
        }
        timer.reset();
        F.createNonZeroPattern(fill_level, Uplo::Upper, UU);
        t_symbolic = timer.seconds();
        if (vtune_symbolic) {
#ifdef HAVE_SHYLUTACHO_VTUNE
          __itt_pause();
#endif
        }

        cout << "CholPerformance:: AA (nnz) = " << AA.NumNonZeros() << ", UU (nnz) = " << UU.NumNonZeros() << endl;

        if (verbose)
          cout << F << endl
               << UU << endl;
      }
      cout << "CholPerformance:: symbolic factorization::time = " << t_symbolic << endl;
      {
        timer.reset();
        CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());

        for (ordinal_type k=0;k<HU.NumNonZeros();++k)
          HU.Value(k).fillRowViewArray();

        t_flat2hier = timer.seconds();

        cout << "CholPerformance:: Hier (dof, nnz) = " << HU.NumRows() << ", " << HU.NumNonZeros() << endl;
      }
      cout << "CholPerformance:: construct hierarchical matrix::time = " << t_flat2hier << endl;
    }

    // copy of UU
    CrsMatrixBaseType RR("RR");
    RR.copy(UU);

    const size_t max_concurrency = 16384;
    cout << "CholPerformance:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "CholPerformance:: max task size   = " << max_task_size << endl;

    if (!skip_serial) {
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   1);

      TaskFactoryType::setUseTeamInterface(team_interface);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType U(&UU);
      U.fillRowViewArray();

      cout << "CholPerformance:: Serial factorize the matrix" << endl;
      {
        UU.copy(RR);
        if (vtune_serial) {
#ifdef HAVE_SHYLUTACHO_VTUNE
          __itt_resume();
#endif
        }
        timer.reset();
        {
          Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
            ::invoke(TaskFactoryType::Policy(),
                              TaskFactoryType::Policy().member_single(),
                              U);
        }
        t_factor_seq = timer.seconds();
        if (vtune_serial) {
#ifdef HAVE_SHYLUTACHO_VTUNE
          __itt_pause();
#endif
        }
        if (verbose)
          cout << UU << endl;
      }
      cout << "CholPerformance:: Serial factorize the matrix::time = " << t_factor_seq << endl;
    }

    {
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   team_size);
      TaskFactoryType::setUseTeamInterface(team_interface);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      cout << "CholPerformance:: ByBlocks factorize the matrix:: team_size = " << team_size << endl;
      CrsHierTaskViewType H(&HU);
      {
        UU.copy(RR);
        if (vtune_task) {
#ifdef HAVE_SHYLUTACHO_VTUNE
          __itt_resume();
#endif
        }
        timer.reset();
        {
          auto future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks>::
                                                              TaskFunctor<CrsHierTaskViewType>(H), 0);
          TaskFactoryType::Policy().spawn(future);
          Kokkos::Experimental::wait(TaskFactoryType::Policy());
        }
        t_factor_task += timer.seconds();
        if (vtune_task) {
#ifdef HAVE_SHYLUTACHO_VTUNE
        __itt_pause();
#endif
        }
        if (verbose)
          cout << UU << endl;
      }
      cout << "CholPerformance:: ByBlocks factorize the matrix::time = " << t_factor_task << endl;
    }

    if (!skip_serial) {
      cout << "CholPerformance:: task scale [seq/task] = " << t_factor_seq/t_factor_task << endl;
    }

    return r_val;
  }
}

#endif
