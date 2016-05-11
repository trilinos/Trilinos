#pragma once
#ifndef __EXAMPLE_CHOL_UNBLOCKED_HPP__
#define __EXAMPLE_CHOL_UNBLOCKED_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

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
  int exampleCholUnblocked(const string file_input,
                           const int max_task_dependence,
                           const int team_size,
                           const int algo,
                           const int variant,
                           const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> CrsMatrixBaseType;
    typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
    
    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "CholUnblocked:: import input file = " << file_input << endl;        
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
    cout << "CholUnblocked:: import input file::time = " << t << endl;        

    const size_t max_concurrency = 10;
    cout << "CholUnblocked:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 3*sizeof(CrsTaskViewType)+128;
    cout << "CholUnblocked:: max task size   = " << max_task_size << endl;


    typename TaskFactoryType::policy_type policy(max_concurrency,
                                                 max_task_size,
                                                 max_task_dependence, 
                                                 team_size);

    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    cout << "CholUnblocked:: factorize the matrix" << endl;
    CrsTaskViewType U(&UU);
    U.fillRowViewArray();
    {
      timer.reset();
    
      typename TaskFactoryType::future_type future;
      switch (algo) {
      case AlgoChol::UnblockedOpt: {
        if (variant == Variant::One)
          future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::One>
                                                         ::TaskFunctor<CrsTaskViewType>(U), 0);
        else if (variant == Variant::Two)
          future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::Two>
                                                         ::TaskFunctor<CrsTaskViewType>(U), 0);
        else {
          ERROR(">> Not supported algorithm variant");          
        }
        break;
      }
      case AlgoChol::Dummy: {
        future = TaskFactoryType::Policy().create_team(Chol<Uplo::Upper,AlgoChol::Dummy>
                                                       ::TaskFunctor<CrsTaskViewType>(U), 0);
        break;
      }
      default:
        ERROR(">> Not supported algorithm");
        break;
      }
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());

      t = timer.seconds();

      if (verbose)
        cout << UU << endl;
    }   
    cout << "CholUnblocked:: factorize the matrix::time = " << t << endl; 
    
    return r_val;
  }
}

#endif
