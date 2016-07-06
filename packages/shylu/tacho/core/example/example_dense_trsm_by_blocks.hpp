#pragma once
#ifndef __EXAMPLE_DENSE_TRSM_BY_BLOCKS_HPP__
#define __EXAMPLE_DENSE_TRSM_BY_BLOCKS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"  
#include "Teuchos_ScalarTraits.hpp" 

#include "util.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"
#include "dense_matrix_helper.hpp"

#include "task_view.hpp"
#include "task_factory.hpp"

#include "trsm.hpp"
#include "dense_flop.hpp"

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
  int exampleDenseTrsmByBlocks(const OrdinalType mmin,
                               const OrdinalType mmax,
                               const OrdinalType minc,
                               const OrdinalType k,
                               const OrdinalType mb,
                               const int max_concurrency,
                               const int max_task_dependence,
                               const int team_size,
                               const int mkl_nthreads,
                               const bool check,
                               const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryTraits> DenseMatrixBaseType;

    typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
    typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

    typedef DenseMatrixBase<DenseTaskViewType,ordinal_type,size_type,SpaceType,MemoryTraits> DenseHierMatrixBaseType;

    typedef DenseMatrixView<DenseHierMatrixBaseType> DenseHierMatrixViewType;
    typedef TaskView<DenseHierMatrixViewType,TaskFactoryType> DenseHierTaskViewType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    cout << "DenseGemmByBlocks:: test matrices "
         <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc << " , k = "<< k << endl;

    const size_t max_task_size = (3*sizeof(DenseTaskViewType)+196); // when 128 error
    typename TaskFactoryType::policy_type policy(max_concurrency,
                                                 max_task_size,
                                                 max_task_dependence, 
                                                 team_size);
    
    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    ostringstream os;
    os.precision(3);
    os << scientific;

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      os.str("");

      DenseMatrixBaseType AA("AA", m, m), BB("BB", m, k), BC("BC", m, k);
      
      // setup upper triangular
      for (ordinal_type j=0;j<AA.NumCols();++j) {
        AA.Value(j,j) = 10.0;
        for (ordinal_type i=0;i<j;++i)
          AA.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
      }

      // setup one and right hand side is going to be overwritten by the product of AB
      for (ordinal_type j=0;j<BB.NumCols();++j)
        for (ordinal_type i=0;i<BB.NumRows();++i)
          BB.Value(i,j) = 1.0;

      Teuchos::BLAS<ordinal_type,value_type> blas;

      blas.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
                m, k, m,
                1.0,
                AA.ValuePtr(), AA.ColStride(),
                BB.ValuePtr(), BB.ColStride(),
                0.0,
                BC.ValuePtr(), BC.ColStride());
      BB.copy(BC);

      const double flop = get_flop_trsm_upper<value_type>(m, k);

#ifdef HAVE_SHYLUTACHO_MKL
      mkl_set_num_threads(mkl_nthreads);
#endif
      os << "DenseTrsmByBlocks:: m = " << m << " k = " << k;
      if (check) {
        timer.reset();
        DenseTaskViewType A(&AA), B(&BC);
        Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,AlgoTrsm::ExternalBlas>::invoke
          (TaskFactoryType::Policy(),
           TaskFactoryType::Policy().member_single(),
           Diag::NonUnit, 1.0, A, B);
        t = timer.seconds();
        os << ":: Serial Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      }

      {
        DenseHierMatrixBaseType HA, HB;
        DenseMatrixHelper::flat2hier(AA, HA, mb, mb);
        DenseMatrixHelper::flat2hier(BB, HB, mb, mb);

        DenseHierTaskViewType TA(&HA), TB(&HB);
        timer.reset();
        auto future = TaskFactoryType::Policy().create_team
          (Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,AlgoTrsm::DenseByBlocks,Variant::One>
           ::TaskFunctor<value_type,DenseHierTaskViewType,DenseHierTaskViewType>
           (Diag::NonUnit, 1.0, TA, TB), 0);
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
        t = timer.seconds();       
        os << ":: Parallel Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      } 

      if (check) {
        typedef typename Teuchos::ScalarTraits<value_type>::magnitudeType real_type; 
        real_type err = 0.0, norm = 0.0;
        for (ordinal_type j=0;j<BB.NumCols();++j)
          for (ordinal_type i=0;i<BB.NumRows();++i) {
            const real_type diff = abs(BB.Value(i,j) - BC.Value(i,j));
            const real_type val  = BC.Value(i,j);
            err  += diff*diff;
            norm += val*val;
          }
        os << ":: Check result ::err = " << sqrt(err) << ", norm = " << sqrt(norm);
      }
      cout << os.str() << endl;
    }

    return r_val;
  }
}

#endif
