#pragma once
#ifndef __EXAMPLE_CHOL_BY_BLOCKS_HPP__
#define __EXAMPLE_CHOL_BY_BLOCKS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
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
  int exampleCholByBlocks(const string file_input,
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

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "CholByBlocks:: import input file = " << file_input << endl;
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
    cout << "CholByBlocks:: import input file::time = " << t << endl;


    cout << "CholByBlocks:: reorder the matrix" << endl;
    CrsMatrixBaseType UU("UU");     // permuted base matrix
    CrsHierMatrixBaseType HU("HU"); // hierarchical matrix of views
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

      for (ordinal_type k=0;k<HU.NumNonZeros();++k)
        HU.Value(k).fillRowViewArray();

      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }
    cout << "CholByBlocks:: reorder the matrix::time = " << t << endl;

    const size_t max_concurrency = 16384;
    cout << "CholByBlocks:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "CholByBlocks:: max task size   = " << max_task_size << endl;

    typename TaskFactoryType::policy_type policy(max_concurrency,
                                                 max_task_size,
                                                 max_task_dependence,
                                                 team_size);

    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    cout << "CholByBlocks:: factorize the matrix" << endl;
    CrsHierTaskViewType H(&HU);
    {
      timer.reset();

      auto future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks>::
                                                          TaskFunctor<CrsHierTaskViewType>(H), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());

      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }
    cout << "CholByBlocks:: factorize the matrix::time = " << t << endl;

    return r_val;
  }
}

#endif
