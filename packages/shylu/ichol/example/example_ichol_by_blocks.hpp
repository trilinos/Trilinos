#pragma once
#ifndef __EXAMPLE_ICHOL_BY_BLOCKS_HPP__
#define __EXAMPLE_ICHOL_BY_BLOCKS_HPP__

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
  int exampleICholByBlocks(const string file_input,
                           const int max_task_dependence,
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
    double t = 0.0;

    cout << "ICholByBlocks:: import input file = " << file_input << endl;        
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
    cout << "ICholByBlocks:: import input file::time = " << t << endl;


    cout << "ICholByBlocks:: reorder the matrix" << endl;        
    CrsMatrixBaseType UU("UU");     // permuted base matrix
    CrsHierMatrixBaseType HU("HU"); // hierarchical matrix of views
    {
      timer.reset();

      GraphHelperType S(AA);
      S.computeOrdering();

      CrsMatrixBaseType PA("Permuted AA");
      PA.copy(S.PermVector(), S.InvPermVector(), AA);

      UU.copy(Uplo::Upper, PA);

      CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 S.TreeVector());
      
      for (ordinal_type k=0;k<HU.NumNonZeros();++k)
        HU.Value(k).fillRowViewArray();
      
      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }
    cout << "ICholByBlocks:: reorder the matrix::time = " << t << endl;            

#ifdef __USE_FIXED_TEAM_SIZE__
    typename TaskFactoryType::policy_type policy(max_task_dependence);
#else
    typename TaskFactoryType::policy_type policy(max_task_dependence, team_size);
#endif
    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    cout << "ICholByBlocks:: factorize the matrix" << endl;
    CrsHierTaskViewType H(&HU);
    {
      timer.reset();

      auto future = TaskFactoryType::Policy().create_team(IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
                                                          TaskFunctor<ForType,CrsHierTaskViewType>(H), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());

      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }  
    cout << "ICholByBlocks:: factorize the matrix::time = " << t << endl;

    return r_val;
  }
}

#endif
