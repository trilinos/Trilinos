#pragma once
#ifndef __EXAMPLE_PARDISO_CHOL_PERFORMANCE_HPP__
#define __EXAMPLE_PARDISO_CHOL_PERFORMANCE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "dense_matrix_base.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#include "example_pardiso.hpp"

namespace Tacho {

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  KOKKOS_INLINE_FUNCTION
  int examplePardisoCholPerformance(const string file_input,
                                    const int nrhs,
                                    const int nthreads,
                                    const bool skip_factorize,
                                    const bool skip_solve,
                                    const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    int r_val = 0;

    // mkl nthreads setting 
    mkl_set_dynamic(0);
    mkl_set_num_threads(nthreads);

    Kokkos::Impl::Timer timer;
    double t = 0.0;
    Pardiso pardiso;

    cout << "PardisoCholPerformance:: init" << endl;
    {
      timer.reset();
      r_val = pardiso.init<value_type,AlgoChol::ExternalPardiso>();
      t = timer.seconds();
      
      if (r_val) 
        cout << "PardisoCholPerformance:: Pardiso init error = " << r_val 
             << pardiso.showErrorCode(cout) << endl;
    }
    cout << "PardisoCholPerformance:: init ::time = " << t << endl;
    
    cout << "PardisoCholPerformance:: import input file = " << file_input << endl;
    CrsMatrixBaseType AA("AA");
    {
      timer.reset();

      CrsMatrixBaseType Atmp("Atmp");
      ifstream in;
      in.open(file_input);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_input << endl;
        return ++r_val;
      }
      Atmp.importMatrixMarket(in);

      AA.copy(Uplo::Upper, Atmp);

      t = timer.seconds();

      if (verbose) {
        cout << Atmp << endl;
        cout << AA << endl;
      }
    }
    cout << "PardisoCholPerformance:: import input file::time = " << t << endl;

    DenseMatrixBaseType BB("BB",  AA.NumRows(), nrhs);
    DenseMatrixBaseType XX("XX",  AA.NumRows(), nrhs);
    DenseMatrixBaseType PP("PP",  AA.NumRows(), 1);
    
    pardiso.setProblem(AA.NumRows(),
                       (double*)AA.ValuePtr(),
                       (int*)AA.RowPtr(),
                       (int*)AA.ColPtr(),
                       (int*)PP.ValuePtr(),
                       BB.NumCols(),
                       (double*)BB.ValuePtr(),
                       (double*)XX.ValuePtr());
    
    cout << "PardisoCholPerformance:: analyze matrix" << endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::Analyze);
      t = timer.seconds();
      
      if (r_val) 
        cout << "PardisoCholPerformance:: Pardiso analyze error = " << r_val 
             << pardiso.showErrorCode(cout) << endl;
      else
        cout << pardiso.showStat(cout, Pardiso::Analyze) << endl;
      cout << "PardisoCholPerformance:: analyze matrix::time = " << t << endl;
    }

    if (!skip_factorize) {
      cout << "PardisoCholPerformance:: factorize matrix" << endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Factorize);
        t = timer.seconds();
        
        if (r_val) 
          cout << "PardisoCholPerformance:: Pardiso factorize error = " << r_val 
               << pardiso.showErrorCode(cout) << endl;
        else
          cout << pardiso.showStat(cout, Pardiso::Factorize) << endl;
      }
      cout << "PardisoCholPerformance:: factorize matrix::time = " << t << endl;
    }

    if (!skip_factorize && !skip_solve) {
      cout << "PardisoCholPerformance:: solve matrix" << endl;
      {
        timer.reset();
        r_val = pardiso.run(Pardiso::Solve);
        t = timer.seconds();
        
        if (r_val) 
          cout << "PardisoCholPerformance:: Pardiso solve error = " << r_val 
               << pardiso.showErrorCode(cout) << endl;
        else
          cout << pardiso.showStat(cout, Pardiso::Solve) << endl;
      }
      cout << "PardisoCholPerformance:: solve matrix::time = " << t << endl;
    }
    
    cout << "PardisoCholPerformance:: release all" << endl;
    {
      timer.reset();
      r_val = pardiso.run(Pardiso::ReleaseAll);
      t = timer.seconds();

      if (r_val) 
        cout << "PardisoCholPerformance:: release error = " << r_val 
             << pardiso.showErrorCode(cout) << endl;
      else
        cout << pardiso.showStat(cout, Pardiso::ReleaseAll) << endl;
    }
    cout << "PardisoCholPerformance:: release all::time = " << t << endl;
    
    return r_val;
  }

}

#endif
#endif
