#pragma once
#ifndef __EXAMPLE_CHOL_DIRECT_SOLVER_HPP__
#define __EXAMPLE_CHOL_DIRECT_SOLVER_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "direct_solver.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int exampleCholDirectSolver(const string file_input,
                              const int nrhs,
                              const int nthreads,
                              const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef DirectSolver<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DirectSolverType;
    typedef typename DirectSolverType::CrsMatrixBaseType CrsMatrixBaseType;
    typedef typename DirectSolverType::DenseMatrixBaseType DenseMatrixBaseType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double
      t_import = 0.0,
      t_init = 0.0,
      t_reorder = 0.0,
      t_analyze = 0.0,
      t_factor = 0.0,
      t_solve = 0.0;

    cout << "CholDirectSolver:: import input file = " << file_input << endl;
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
      CrsMatrixHelper::filterZeros(AA);
    }
    cout << "CholDirectSolver:: import input file::time = " << t_import << endl;

    DenseMatrixBaseType BB("BB", AA.NumRows(), nrhs);
    {
      const value_type one(1);
      for (ordinal_type j=0;j<BB.NumCols();++j)
        for (ordinal_type i=0;i<BB.NumRows();++i) 
          BB.Value(i,j) = one;
    }

    DirectSolverType solver("DirectSolver", AA);
    solver.setB(BB);

    timer.reset();
    solver.setDefaultParameters(nthreads);
    solver.init();
    t_init = timer.seconds();
    cout << "CholDirectSolver:: init::time = " << t_init << endl;

    timer.reset();
    solver.reorder();
    t_reorder = timer.seconds();
    cout << "CholDirectSolver:: reorder::time = " << t_reorder << endl;

    timer.reset();
    solver.analyze();
    t_analyze = timer.seconds();
    cout << "CholDirectSolver:: analyze::time = " << t_analyze << endl;

    timer.reset();
    solver.factorize();
    t_factor = timer.seconds();
    cout << "CholDirectSolver:: factor::time = " << t_factor << endl;

    timer.reset();
    solver.solve();
    t_solve = timer.seconds();
    cout << "CholDirectSolver:: solve::time = " << t_solve << endl;

    return r_val;
  }
}

#endif
