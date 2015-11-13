#pragma once
#ifndef __EXAMPLE_CHOL_DIRECT_PERFORMANCE_HPP__
#define __EXAMPLE_CHOL_DIRECT_PERFORMANCE_HPP__

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
#include "mkl_service.h"
#endif

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleCholDirectPerformance(const string file_input,
                                   const int treecut,
                                   const int minblksize,
                                   const int prunecut,
                                   const int seed,
                                   const int niter,
                                   const int nthreads,
                                   const int max_task_dependence,
                                   const int team_size,
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
      t_mkl = 0.0,
#endif
      t_factor_seq = 0.0,   t_solve_seq = 0.0,
      t_factor_task = 0.0,  t_solve_task = 0.0;
    const int start = 0;
    
    cout << "CholDirectPerformance:: import input file = " << file_input << endl;        
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
    }
    cout << "CholDirectPerformance:: import input file::time = " << t_import << endl;

    cout << "CholDirectPerformance:: reorder the matrix" << endl;        
    CrsMatrixBaseType PA("Permuted AA");
    CrsMatrixBaseType UU("UU");     // permuted base upper triangular matrix
    CrsHierMatrixBaseType HU("HU"); // hierarchical matrix of views

    DenseMatrixBaseType BB("BB", AA.NumRows(), nrhs);
    DenseHierMatrixBaseType HB("HB");

    {
      GraphHelperType S(AA, seed);
      {
        timer.reset();
        
        S.computeOrdering(treecut, minblksize);
        S.pruneTree(prunecut);
        
        PA.copy(S.PermVector(), S.InvPermVector(), AA);
        
        t_reorder = timer.seconds();
      }

      cout << "CholDirectPerformance:: reorder the matrix::time = " << t_reorder << endl;            
      {
        SymbolicFactorHelperType F(PA, league_size);
        for (int i=start;i<niter;++i) {
          timer.reset();
        
          F.createNonZeroPattern(Uplo::Upper, UU);

          t_symbolic += timer.seconds() * (i>=0);
        }
        t_symbolic /= niter;

        cout << "CholDirectPerformance:: AA (nnz) = " << AA.NumNonZeros() << ", UU (nnz) = " << UU.NumNonZeros() << endl;
      }
      cout << "CholDirectPerformance:: symbolic factorization::time = " << t_symbolic << endl;            
      {
        timer.reset();

        CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());
        
        for (ordinal_type k=0;k<HU.NumNonZeros();++k)
          HU.Value(k).fillRowViewArray();

        DenseMatrixHelper::flat2hier(BB, HB,
                                     S.NumBlocks(),
                                     S.RangeVector(),
                                     nb);
        
        t_flat2hier = timer.seconds();
        
        cout << "CholDirectPerformance:: Hier (dof, nnz) = " << HU.NumRows() << ", " << HU.NumNonZeros() << endl;
      }
      cout << "CholDirectPerformance:: construct hierarchical matrix::time = " << t_flat2hier << endl;            
    }

    // copy of UU
    CrsMatrixBaseType RR("RR");
    RR.copy(UU);

    /////////////////////////// Serial Numeric Factorization
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

      cout << "CholDirectPerformance:: Serial factorize the matrix" << endl;
      {
        for (int i=start;i<niter;++i) {
          UU.copy(RR);
          timer.reset();          
          {
            Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
              ::invoke(TaskFactoryType::Policy(),
                       TaskFactoryType::Policy().member_single(),
                       U);
          }
          t_factor_seq += timer.seconds() * (i>=0);
        }
        t_factor_seq /= niter;
      }
      cout << "CholDirectPerformance:: Serial factorize the matrix::time = " << t_factor_seq << endl;

      cout << "CholDirectPerformance:: Serial forward/backward solve" << endl;      
      {
        for (int i=start;i<niter;++i) {
          XX.copy(BB);
          timer.reset();          
          {
            TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::Unblocked>
              ::invoke(TaskFactoryType::Policy(),
                       TaskFactoryType::Policy().member_single(),
                       Diag::NonUnit, U, X);
            TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::Unblocked>
              ::invoke(TaskFactoryType::Policy(),
                       TaskFactoryType::Policy().member_single(),
                       Diag::NonUnit, U, X);
          }
          t_factor_seq += timer.seconds() * (i>=0);
        }
        t_factor_seq /= niter;
      }
      cout << "CholDirectPerformance:: Serial forward/backward solve::time = " << t_solve_seq << endl;
    }

    // if (!skip_serial) 
    //   cout << "CholDirectPerformance:: task scale [seq/task] = " << t_factor_seq/t_factor_task << endl;    

    return r_val;
  }
}

#endif
