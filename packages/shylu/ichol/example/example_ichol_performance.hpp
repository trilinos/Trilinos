#pragma once
#ifndef __EXAMPLE_ICHOL_PERFORMANCE_HPP__
#define __EXAMPLE_ICHOL_PERFORMANCE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "symbolic_factor_helper.hpp"
#include "crs_matrix_helper.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "parallel_for.hpp"

#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "ichol.hpp"

namespace Example {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleICholPerformance(const string file_input,
                              const int seed,
                              const int niter,
                              const int nthreads,
                              const int max_task_dependence,
                              const int team_size,
                              const int fill_level,
                              const int league_size,
                              const bool team_interface,
                              const bool skip_serial,
                              const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType>,
      Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;
    typedef ParallelFor ForType;

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
    const int start = -2;
    
    cout << "ICholPerformance:: import input file = " << file_input << endl;        
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
    cout << "ICholPerformance:: import input file::time = " << t_import << endl;

    cout << "ICholPerformance:: reorder the matrix" << endl;        
    CrsMatrixBaseType PA("Permuted AA");
    CrsMatrixBaseType UU("UU");     // permuted base upper triangular matrix
    CrsHierMatrixBaseType HU("HU"); // hierarchical matrix of views

    {
      GraphHelperType S(AA, seed);
      {
        timer.reset();
        
        S.computeOrdering();
        
        PA.copy(S.PermVector(), S.InvPermVector(), AA);
        
        t_reorder = timer.seconds();
      }
      cout << "ICholPerformance:: reorder the matrix::time = " << t_reorder << endl;            
      {
        for (int i=start;i<niter;++i) {
          timer.reset();
        
          SymbolicFactorHelperType F(PA, league_size);
          F.createNonZeroPattern(fill_level, Uplo::Upper, UU);

          // UU.copy(Uplo::Upper, PA);

          t_symbolic += timer.seconds() * (i>=0);
        }
        t_symbolic /= niter;

        cout << "ICholPerformance:: AA (nnz) = " << AA.NumNonZeros() << ", UU (nnz) = " << UU.NumNonZeros() << endl;

        if (verbose)
          cout << UU << endl;
      }
      cout << "ICholPerformance:: symbolic factorization::time = " << t_symbolic << endl;            
      {
        timer.reset();

        CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());
        
        for (ordinal_type k=0;k<HU.NumNonZeros();++k)
          HU.Value(k).fillRowViewArray();
        
        t_flat2hier = timer.seconds();
        
        cout << "ICholPerformance:: Hier (dof, nnz) = " << HU.NumRows() << ", " << HU.NumNonZeros() << endl;
      }
      cout << "ICholPerformance:: construct hierarchical matrix::time = " << t_flat2hier << endl;            
    }

    // copy of UU
    CrsMatrixBaseType RR("RR");
    RR.copy(UU);

    if (!skip_serial) {
#ifdef __USE_FIXED_TEAM_SIZE__
      typename TaskFactoryType::policy_type policy(max_task_dependence);
#else
      typename TaskFactoryType::policy_type policy(max_task_dependence, 1);
#endif
      TaskFactoryType::setUseTeamInterface(team_interface);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType U(&UU);
      U.fillRowViewArray();

      cout << "ICholPerformance:: Serial factorize the matrix" << endl;
      {
        for (int i=start;i<niter;++i) {
          UU.copy(RR);
          timer.reset();          
          // {
          //   auto future = TaskFactoryType::Policy().create(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
          //                                                  ::TaskFunctor<ForType,CrsTaskViewType>(U), 0);
          //   TaskFactoryType::Policy().spawn(future);
          //   Kokkos::Experimental::wait(TaskFactoryType::Policy());
          // }
          {
            IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
              ::invoke<ForType>(TaskFactoryType::Policy(),
                                TaskFactoryType::Policy().member_single(),
                                U);
          }
          t_factor_seq += timer.seconds() * (i>=0);
        }
        t_factor_seq /= niter;

        if (verbose)
          cout << UU << endl;
      }
      cout << "ICholPerformance:: Serial factorize the matrix::time = " << t_factor_seq << endl;
    }

//     {
// #ifdef __USE_FIXED_TEAM_SIZE__
//       typename TaskFactoryType::policy_type policy(max_task_dependence);
// #else
//       typename TaskFactoryType::policy_type policy(max_task_dependence, nthreads);
// #endif
//       TaskFactoryType::setPolicy(&policy);

//       CrsTaskViewType U(&UU);
//       U.fillRowViewArray();

//       cout << "ICholPerformance:: Team factorize the matrix:: team_size = " << nthreads << endl;
//       {
//         timer.reset();
        
//         auto future = TaskFactoryType::Policy().create(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
//                                                        ::TaskFunctor<ForType,CrsTaskViewType>(U), 0);
//         TaskFactoryType::Policy().spawn(future);
//         Kokkos::Experimental::wait(TaskFactoryType::Policy());
        
//         t_factor_team = timer.seconds();
        
//         if (verbose)
//           cout << UU << endl;
//       }
//       cout << "ICholPerformance:: Team factorize the matrix::time = " << t_factor_team << endl;
//     }

    {
#ifdef __USE_FIXED_TEAM_SIZE__
      typename TaskFactoryType::policy_type policy(max_task_dependence);
#else
      typename TaskFactoryType::policy_type policy(max_task_dependence, team_size);
#endif
      TaskFactoryType::setUseTeamInterface(team_interface);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);
      
      cout << "ICholPerformance:: ByBlocks factorize the matrix:: team_size = " << team_size << endl;
      CrsHierTaskViewType H(&HU);
      {
        for (int i=start;i<niter;++i) {
          UU.copy(RR);
          timer.reset();
          {
            auto future = TaskFactoryType::Policy().create_team(IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
                                                                TaskFunctor<ForType,CrsHierTaskViewType>(H), 0);
            TaskFactoryType::Policy().spawn(future);
            Kokkos::Experimental::wait(TaskFactoryType::Policy());
          }
          t_factor_task += timer.seconds() * (i>=0);
        }
        t_factor_task /= niter;

        if (verbose)
          cout << UU << endl;
      }  
      cout << "ICholPerformance:: ByBlocks factorize the matrix::time = " << t_factor_task << endl;
    }

    if (!skip_serial) {
      cout << "ICholPerformance:: task scale [seq/task] = " << t_factor_seq/t_factor_task << endl;    
      //cout << "ICholPerformance:: team scale [seq/team] = " << t_factor_seq/t_factor_team << endl;    
    }

    return r_val;
  }
}

#endif
