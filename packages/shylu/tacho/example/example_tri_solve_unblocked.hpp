#pragma once
#ifndef __EXAMPLE_TRI_SOLVE_UNBLOCKED_HPP__
#define __EXAMPLE_TRI_SOLVE_UNBLOCKED_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "task_view.hpp"

#include "task_factory.hpp"

#include "tri_solve.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleTriSolveUnblocked(const string file_input, 
                               const OrdinalType nrhs,
                               const int max_task_dependence,
                               const int team_size,
                               const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;
    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    
    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "TriSolveUnblocked:: import input file = " << file_input << endl;        
    CrsMatrixBaseType AA("AA"), UU("UU");    
    {
      timer.reset();

      ifstream in;
      in.open(file_input);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_input << endl;
        return ++r_val;
      }
      AA.importMatrixMarket(in);
      
      UU.copy(Uplo::Upper, AA);

      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }
    cout << "TriSolveUnblocked:: import input file::time = " << t << endl;

    cout << "TriSolveUnblocked:: create right hand side = " << nrhs << endl;        
    DenseMatrixBaseType BB("BB", UU.NumRows(), nrhs); 
    {
      timer.reset();

      const value_type one(1);
      for (ordinal_type j=0;j<BB.NumCols();++j)
        for (ordinal_type i=0;i<BB.NumRows();++i) 
          BB.Value(i, j) = one;

      t = timer.seconds();      

      if (verbose)
        cout << BB << endl;
    }
    cout << "TriSolveUnblocked:: create right hand side::time = " << t << endl;

    const size_t max_concurrency = 10;
    cout << "CholPerformance:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "CholPerformance:: max task size   = " << max_task_size << endl;

    typename TaskFactoryType::policy_type policy(max_concurrency,
                                                 max_task_size,
                                                 max_task_dependence, 
                                                 team_size);

    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    cout << "TriSolveUnblocked:: perform forward and backward solve of the matrix" << endl;
    CrsTaskViewType U(&UU);
    DenseTaskViewType B(&BB);
    U.fillRowViewArray();
    {
      timer.reset();

      {
        auto future = TaskFactoryType::Policy().create_team(TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::Unblocked>
                                                            ::TaskFunctor<CrsTaskViewType,DenseTaskViewType>
                                                            (Diag::NonUnit, U, B), 0);
        
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
      }
      {
        auto future = TaskFactoryType::Policy().create_team(TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::Unblocked>
                                                            ::TaskFunctor<CrsTaskViewType,DenseTaskViewType>
                                                            (Diag::NonUnit, U, B), 0);
        
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
      }

      t = timer.seconds();
      
      if (verbose)
        cout << BB << endl;
    }
    cout << "TriSolveUnblocked:: perform forward and backward solve of the matrix::time = " << t << endl;

    return r_val;
  }
}

#endif
