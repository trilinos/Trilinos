#pragma once
#ifndef __EXAMPLE_CHOL_DIRECT_PLAIN_HPP__
#define __EXAMPLE_CHOL_DIRECT_PLAIN_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "graph_helper_scotch.hpp"
#include "symbolic_factor_helper.hpp"
#include "crs_matrix_helper.hpp"
#include "dense_matrix_helper.hpp"

#include "task_view.hpp"
#include "task_factory.hpp"

#include "chol.hpp"
#include "tri_solve.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleCholDirectPlain(const string file_input,
                             const int prunecut,
                             const int seed,
                             const int nrhs,
                             const int nb,
                             const int nthreads,
                             const int max_concurrency,
                             const int max_task_dependence,
                             const int team_size,
                             const int league_size,
                             const bool team_interface,
                             const bool serial,
                             const bool solve,
                             const bool check,
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

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    typedef DenseMatrixBase<DenseTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> DenseHierMatrixBaseType;

    typedef DenseMatrixView<DenseHierMatrixBaseType> DenseHierMatrixViewType;
    typedef TaskView<DenseHierMatrixViewType,TaskFactoryType> DenseHierTaskViewType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double
      t_import = 0.0,
      t_reorder = 0.0,
      t_symbolic = 0.0,
      t_flat2hier = 0.0,
      t_factor = 0.0,
      t_solve = 0.0;

    cout << "CholDirectPlain:: import input file = " << file_input << endl;
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

      cout << "CholDirectPlain:: input nnz = " << AA.NumNonZeros() << endl;
      CrsMatrixHelper::filterZeros(AA);
      cout << "CholDirectPlain:: resized nnz = " << AA.NumNonZeros() << endl;
    }
    cout << "CholDirectPlain:: import input file::time = " << t_import << endl;

    // matrix A and its upper triangular factors U
    CrsMatrixBaseType PA("Permuted AA");
    CrsMatrixBaseType UU("UU"); // permuted base upper triangular matrix
    CrsHierMatrixBaseType HU("HU"); // hierarchical matrix of views

    // right hand side and solution matrix
    DenseMatrixBaseType BB("BB", AA.NumRows(), nrhs), XX("XX", AA.NumRows(), nrhs);
    DenseHierMatrixBaseType HB("HB"), HX("HX");

    {
      cout << "CholDirectPlain:: reorder the matrix" << endl;

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
        S.computeOrdering();
        S.pruneTree(prunecut);
        PA.copy(S.PermVector(), S.InvPermVector(), AA);
        t_reorder = timer.seconds();
      }
      cout << "CholDirectPlain:: reorder the matrix::time = " << t_reorder << endl;

      {
        SymbolicFactorHelperType F(PA, league_size);
        timer.reset();
        F.createNonZeroPattern(Uplo::Upper, UU);
        t_symbolic = timer.seconds();
        cout << "CholDirectPlain:: AA (nrow, nnz) = " << AA.NumRows() << ", " << AA.NumNonZeros() << ", UU (nnz) = " << UU.NumNonZeros() << endl;
      }
      cout << "CholDirectPlain:: symbolic factorization::time = " << t_symbolic << endl;

      {
        timer.reset();
        CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());

        // fill internal meta data for sparse blocs
        size_type nnz = 0;
        for (ordinal_type k=0;k<HU.NumNonZeros();++k) {
          HU.Value(k).fillRowViewArray();
          nnz += HU.Value(k).countNumNonZeros();
        }
        if (nnz != UU.NumNonZeros())
          ERROR(">> 2D Blocks do not cover the flat matrix");

        DenseMatrixHelper::flat2hier(BB, HB,
                                     S.NumBlocks(),
                                     S.RangeVector(),
                                     nb);

        DenseMatrixHelper::flat2hier(XX, HX,
                                     S.NumBlocks(),
                                     S.RangeVector(),
                                     nb);
        t_flat2hier = timer.seconds();
        cout << "CholDirectPlain:: Hier (dof, nnz) = " << HU.NumRows() << ", " << HU.NumNonZeros() << endl;
      }
      cout << "CholDirectPlain:: construct hierarchical matrix::time = " << t_flat2hier << endl;
    }


    {
      cout << "CholDirectPlain:: max concurrency = " << max_concurrency << endl;

      const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
      cout << "CholDirectPlain:: max task size   = " << max_task_size << endl;

      // Policy setup
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   team_size);

      TaskFactoryType::setUseTeamInterface(team_interface);
      TaskFactoryType::setMaxTaskDependence(max_task_dependence);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType A(&PA), U(&UU);
      DenseTaskViewType X(&XX), B(&BB);

      A.fillRowViewArray();
      U.fillRowViewArray();

      CrsHierTaskViewType TU(&HU);
      DenseHierTaskViewType TB(&HB), TX(&HX);

      {
        // Manufacture B = AX
        const int m = A.NumRows();
        for (int j=0;j<nrhs;++j)
          for (int i=0;i<m;++i)
            X.Value(i,j) = (j+1);

        Gemm<Trans::NoTranspose,Trans::NoTranspose,AlgoGemm::ForTriSolveBlocked>
          ::invoke(TaskFactoryType::Policy(),
                   TaskFactoryType::Policy().member_single(),
                   1.0, A, X, 0.0, B);
        XX.copy(BB);
      }

      if (serial) {
        cout << "CholDirectPlain:: Serial factorize the matrix" << endl;
        timer.reset();
        Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
          ::invoke(TaskFactoryType::Policy(),
                   TaskFactoryType::Policy().member_single(),
                   U);
        t_factor = timer.seconds();
        cout << "CholDirectPlain:: Serial factorize the matrix::time = " << t_factor << endl;
      } else {
        cout << "CholDirectPlain:: ByBlocks factorize the matrix:: team_size = " << team_size << endl;
        timer.reset();
        auto future = TaskFactoryType::Policy().create_team
          (Chol<Uplo::Upper,AlgoChol::ByBlocks>
           ::TaskFunctor<CrsHierTaskViewType>(TU), 0);
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
        t_factor = timer.seconds();
        cout << "CholDirectPlain:: ByBlocks factorize the matrix::time = " << t_factor << endl;
      }

      if (solve) {
        if (serial) {
          cout << "CholDirectPlain:: Serial forward/backward solve" << endl;
          timer.reset();
          TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::Unblocked>
            ::invoke(TaskFactoryType::Policy(),
                     TaskFactoryType::Policy().member_single(),
                     Diag::NonUnit, U, X);
          TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::Unblocked>
            ::invoke(TaskFactoryType::Policy(),
                     TaskFactoryType::Policy().member_single(),
                     Diag::NonUnit, U, X);
          t_solve = timer.seconds();
          cout << "CholDirectPlain:: Serial forward/backward solve::time = " << t_solve << endl;
        } else {
          cout << "CholDirectPlain:: ByBlocks forward/backward solve" << endl;
          timer.reset();
          auto future_forward_solve = TaskFactoryType::Policy().create_team
            (TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::ByBlocks>
             ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
             (Diag::NonUnit, TU, TX), 0);

          TaskFactoryType::Policy().spawn(future_forward_solve);

          auto future_backward_solve = TaskFactoryType::Policy().create_team
            (TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::ByBlocks>
             ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
             (Diag::NonUnit, TU, TX), 1);

          TaskFactoryType::Policy().add_dependence(future_backward_solve, future_forward_solve);
          TaskFactoryType::Policy().spawn(future_backward_solve);

          Kokkos::Experimental::wait(TaskFactoryType::Policy());
          t_solve = timer.seconds();
          cout << "CholDirectPlain:: ByBlocks forward/backward solve::time = " << t_solve << endl;
        }
      }

      if (solve && check) {
        // Check manufactured solution
        double l2 = 0.0, linf = 0.0;
        const int m = A.NumRows();
        for (int j=0;j<nrhs;++j)
          for (int i=0;i<m;++i) {
            double diff = abs(X.Value(i,j) - (j+1));
            l2 += diff*diff;
            linf = max(diff, linf);
          }
        l2 = sqrt(l2);
        cout << "CholDirectPlain:: Check solution::L2 = " << l2 << ", Linf = " << linf << endl;
      }
    }

    return r_val;
  }
}

#endif
