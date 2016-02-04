#pragma once
#ifndef __EXAMPLE_CHOL_PERFORMANCE_HPP__
#define __EXAMPLE_CHOL_PERFORMANCE_HPP__

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

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_rci.h"
#endif

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleCholPerformance(const string file_input,
                              const int treecut,
                              const int minblksize,
                              const int prunecut,
                              const int seed,
                              const int niter,
                              const int nthreads,
                              const int max_task_dependence,
                              const int team_size,
                              const int fill_level,
                              const int league_size,
                              const bool team_interface,
                              const bool skip_serial,
                              const bool mkl_interface,
                              const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;
    typedef SymbolicFactorHelper<CrsMatrixBaseType> SymbolicFactorHelperType;

#ifdef HAVE_SHYLUTACHO_MKL
    typedef typename CrsMatrixBaseType::value_type_array value_type_array;
#endif

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
#ifdef HAVE_SHYLUTACHO_MKL
      t_mkl_seq = 0.0,
#endif
      t_factor_seq = 0.0,
      t_factor_task = 0.0;
    const int start = -2;

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
        for (int i=start;i<niter;++i) {
          timer.reset();

          F.createNonZeroPattern(fill_level, Uplo::Upper, UU);

          // UU.copy(Uplo::Upper, PA);

          t_symbolic += timer.seconds() * (i>=0);
        }
        t_symbolic /= niter;

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

#ifdef HAVE_SHYLUTACHO_MKL
    if (!skip_serial && mkl_interface) {
      cout << "CholPerformance:: MKL factorize the matrix" << endl;
      CrsMatrixBaseType MM("MM");
      for (int i=start;i<niter;++i) {
        MM.copy(RR);
        MM.hermitianize(Uplo::Upper);

        MKL_INT  n  = static_cast<MKL_INT>(MM.NumRows());
        double  *a  = static_cast<double*>(MM.ValuePtr());
        MKL_INT *ia = static_cast<MKL_INT*>(MM.RowPtr());
        MKL_INT *ja = static_cast<MKL_INT*>(MM.ColPtr());

        // convert to 1-based matrix
        {
          for (ordinal_type k=0;k<(MM.NumRows()+1);++k)
            ++ia[k];

          for (size_type k=0;k<MM.NumNonZeros();++k)
            ++ja[k];
        }
        value_type_array mkl_result = value_type_array("mkl-ilu-values", MM.NumNonZeros());
        double  *bilu0 = static_cast<double*>(&mkl_result[0]);
        MKL_INT ipar[128];
        double dpar[128];
        MKL_INT ierr;

        // we provide ilu-k pattern
        timer.reset();
        dcsrilu0(&n, a, ia, ja, bilu0, ipar, dpar, &ierr);
        t_mkl_seq += timer.seconds() * (i>=0) * 0.5;

        if (ierr != 0)
          cout << " MKL Error = " << ierr << endl;
      }
      t_mkl_seq /= niter;
      cout << "CholPerformance:: MKL factorize the matrix::time = " << t_mkl_seq << endl;
    }
#endif

    const size_t max_concurrency = 16384;
    cout << "CholPerformance:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "CholPerformance:: max task size   = " << max_task_size << endl;

    if (!skip_serial) {
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   team_size);

      TaskFactoryType::setUseTeamInterface(team_interface);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType U(&UU);
      U.fillRowViewArray();

      cout << "CholPerformance:: Serial factorize the matrix" << endl;
      {
        for (int i=start;i<niter;++i) {
          UU.copy(RR);
          timer.reset();
          // {
          //   auto future = TaskFactoryType::Policy().create(Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
          //                                                  ::TaskFunctor<CrsTaskViewType>(U), 0);
          //   TaskFactoryType::Policy().spawn(future);
          //   Kokkos::Experimental::wait(TaskFactoryType::Policy());
          // }
          {
            Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
              ::invoke(TaskFactoryType::Policy(),
                                TaskFactoryType::Policy().member_single(),
                                U);
          }
          t_factor_seq += timer.seconds() * (i>=0);
        }
        t_factor_seq /= niter;

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
        for (int i=start;i<niter;++i) {
          UU.copy(RR);
          timer.reset();
          {
            auto future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks>::
                                                                TaskFunctor<CrsHierTaskViewType>(H), 0);
            TaskFactoryType::Policy().spawn(future);
            Kokkos::Experimental::wait(TaskFactoryType::Policy());
          }
          t_factor_task += timer.seconds() * (i>=0);
        }
        t_factor_task /= niter;

        if (verbose)
          cout << UU << endl;
      }
      cout << "CholPerformance:: ByBlocks factorize the matrix::time = " << t_factor_task << endl;
    }

    if (!skip_serial) {
#ifdef HAVE_SHYLUTACHO_MKL
      cout << "CholPerformance:: mkl/chol scale [mkl/chol] = " << t_mkl_seq/t_factor_seq << endl;
      cout << "CholPerformance:: mkl/task  scale [mkl/task]  = " << t_mkl_seq/t_factor_task << endl;
#else
      cout << "CholPerformance:: task scale [seq/task] = " << t_factor_seq/t_factor_task << endl;
#endif
    }

    return r_val;
  }
}

#endif
