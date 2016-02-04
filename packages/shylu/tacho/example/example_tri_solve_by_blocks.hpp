#pragma once
#ifndef __EXAMPLE_TRI_SOLVE_BY_BLOCKS_HPP__
#define __EXAMPLE_TRI_SOLVE_BY_BLOCKS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "graph_helper_scotch.hpp"
#include "crs_matrix_helper.hpp"
#include "dense_matrix_helper.hpp"

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
  int exampleTriSolveByBlocks(const string file_input,
                              const OrdinalType nrhs,
                              const OrdinalType nb,
                              const int nthreads,
                              const int max_task_dependence,
                              const int team_size,
                              const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    typedef DenseMatrixBase<DenseTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> DenseHierMatrixBaseType;

    typedef DenseMatrixView<DenseHierMatrixBaseType> DenseHierMatrixViewType;
    typedef TaskView<DenseHierMatrixViewType,TaskFactoryType> DenseHierTaskViewType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "TriSolveByBlocks:: import input file = " << file_input << endl;
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
    cout << "TriSolveByBlocks:: import input file::time = " << t << endl;

    cout << "TriSolveByBlocks:: create right hand side = " << nrhs << endl;
    CrsMatrixBaseType   UU("UU");
    DenseMatrixBaseType BB("BB",  AA.NumRows(), nrhs);
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
    cout << "TriSolveByBlocks:: create right hand side::time = " << t << endl;

    cout << "TriSolveByBlocks:: reorder the matrix and partition right hand side, nb = " << nb << endl;
    CrsHierMatrixBaseType   HU("HU");
    DenseHierMatrixBaseType HB("HB");
    {
      timer.reset();

      typename GraphHelperType::size_type_array rptr(AA.Label()+"Graph::RowPtrArray", AA.NumRows() + 1);
      typename GraphHelperType::ordinal_type_array cidx(AA.Label()+"Graph::ColIndexArray", AA.NumNonZeros());

      AA.convertGraph(rptr, cidx);
      GraphHelperType S(AA.Label()+"ScotchHelper",
                        AA.NumRows(),
                        rptr,
                        cidx);
      S.computeOrdering();

      CrsMatrixBaseType PA("Permuted AA");
      PA.copy(S.PermVector(), S.InvPermVector(), AA);

      UU.copy(Uplo::Upper, PA);

      CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 S.TreeVector());

      DenseMatrixHelper::flat2hier(BB, HB,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   nb);

      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }
    cout << "TriSolveByBlocks:: reorder the matrix and partition right hand side::time = " << t << endl;

    const size_t max_concurrency = 16384;
    cout << "TriSolveByBlocks:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "TriSolveByBlocks:: max task size   = " << max_task_size << endl;

    typename TaskFactoryType::policy_type policy(max_concurrency,
                                                 max_task_size,
                                                 max_task_dependence,
                                                 team_size);

    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    // wrap the hierarchically partitioned matrix with task handler
    CrsHierTaskViewType TU(&HU);
    for (ordinal_type k=0;k<HU.NumNonZeros();++k)
      HU.Value(k).fillRowViewArray();

    DenseHierTaskViewType TB(&HB);

    cout << "TriSolveByBlocks:: perform forward and backward solve of the matrix" << endl;
    {
      timer.reset();

      auto future_forward_solve = TaskFactoryType::Policy().create_team
        (TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::ByBlocks>
         ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
         (Diag::NonUnit, TU, TB), 0);

      TaskFactoryType::Policy().spawn(future_forward_solve);

      auto future_backward_solve = TaskFactoryType::Policy().create_team
        (TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::ByBlocks>
         ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
         (Diag::NonUnit, TU, TB), 1);

      TaskFactoryType::Policy().add_dependence(future_backward_solve, future_forward_solve);
      TaskFactoryType::Policy().spawn(future_backward_solve);

      Kokkos::Experimental::wait(TaskFactoryType::Policy());

      t = timer.seconds();

      if (verbose)
        cout << BB << endl;
    }
    cout << "TriSolveByBlocks:: perform forward and backward solve of the matrix::time = " << t << endl;

    return r_val;
  }
}

#endif
