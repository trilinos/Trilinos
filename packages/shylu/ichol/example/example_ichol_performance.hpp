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
                              const int nthreads,
                              const int max_task_dependences,
                              const int team_size,
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
      t_factor_seq = 0.0,
      //t_factor_team = 0.0,
      t_factor_task = 0.0;
    
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
      timer.reset();

      GraphHelperType S(AA);
      S.computeOrdering();

      PA.copy(S.PermVector(), S.InvPermVector(), AA);
      UU.copy(Uplo::Upper, PA);

      CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 S.TreeVector());
      
      for (ordinal_type k=0;k<HU.NumNonZeros();++k)
        HU.Value(k).fillRowViewArray();
      
      t_reorder = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }
    cout << "ICholPerformance:: reorder the matrix::time = " << t_reorder << endl;            

    { // warm up
#ifdef __USE_SERIAL_EXEC_SPACE__
      typename TaskFactoryType::policy_type policy(max_task_dependences);
#else
      typename TaskFactoryType::policy_type policy(max_task_dependences, 1);
#endif
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType U(&UU);
      U.fillRowViewArray();

      {
        auto future = TaskFactoryType::Policy().create(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
                                                       ::TaskFunctor<ForType,CrsTaskViewType>(U), 0);
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
      }
    }

    {
      UU.copy(Uplo::Upper, PA);
#ifdef __USE_SERIAL_EXEC_SPACE__
      typename TaskFactoryType::policy_type policy(max_task_dependences);
#else
      typename TaskFactoryType::policy_type policy(max_task_dependences, 1);
#endif
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType U(&UU);
      U.fillRowViewArray();

      cout << "ICholPerformance:: Serial factorize the matrix" << endl;
      {
        timer.reset();
        
        auto future = TaskFactoryType::Policy().create(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
                                                       ::TaskFunctor<ForType,CrsTaskViewType>(U), 0);
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
        
        t_factor_seq = timer.seconds();
        
        if (verbose)
          cout << UU << endl;
      }
      cout << "ICholPerformance:: Serial factorize the matrix::time = " << t_factor_seq << endl;
    }

//     {
//       UU.copy(Uplo::Upper, PA);
// #ifdef __USE_SERIAL_EXEC_SPACE__
//       typename TaskFactoryType::policy_type policy(max_task_dependences);
// #else
//       typename TaskFactoryType::policy_type policy(max_task_dependences, nthreads);
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
      UU.copy(Uplo::Upper, PA);
#ifdef __USE_SERIAL_EXEC_SPACE__
      typename TaskFactoryType::policy_type policy(max_task_dependences);
#else
      typename TaskFactoryType::policy_type policy(max_task_dependences, team_size);
#endif
      TaskFactoryType::setPolicy(&policy);
      
      cout << "ICholPerformance:: ByBlocks factorize the matrix:: team_size = " << team_size << endl;
      CrsHierTaskViewType H(&HU);
      {
        timer.reset();
        
        auto future = TaskFactoryType::Policy().create_team(IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
                                                            TaskFunctor<ForType,CrsHierTaskViewType>(H), 0);
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
        
        t_factor_task = timer.seconds();
        
        if (verbose)
          cout << UU << endl;
      }  
      cout << "ICholPerformance:: ByBlocks factorize the matrix::time = " << t_factor_task << endl;
    }

    cout << "ICholPerformance:: task scale [seq/task] = " << t_factor_seq/t_factor_task << endl;    
    //cout << "ICholPerformance:: team scale [seq/team] = " << t_factor_seq/t_factor_team << endl;    

    return r_val;
  }
}

#endif
