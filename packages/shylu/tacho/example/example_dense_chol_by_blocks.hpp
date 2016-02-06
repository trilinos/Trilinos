#pragma once
#ifndef __EXAMPLE_DENSE_CHOL_BY_BLOCKS_HPP__
#define __EXAMPLE_DENSE_CHOL_BY_BLOCKS_HPP__

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

#include "chol.hpp"
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
  int exampleDenseCholByBlocks(const OrdinalType mmin,
                               const OrdinalType mmax,
                               const OrdinalType minc,
                               const OrdinalType mb,
                               const int max_concurrency,
                               const int max_task_dependence,
                               const int team_size,
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

    cout << "DenseCholByBlocks:: test matrices "
         <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc << endl;

    const size_t max_task_size = (3*sizeof(DenseTaskViewType)+196); // when 128 error
    //cout << "max task size = "<< max_task_size << endl;
    typename TaskFactoryType::policy_type policy(max_concurrency,
                                                 max_task_size,
                                                 max_task_dependence,
                                                 team_size);

    TaskFactoryType::setMaxTaskDependence(max_task_dependence);
    TaskFactoryType::setPolicy(&policy);

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      DenseMatrixBaseType TT("TT", m, m), AA("AA", m, m), AB("AB", m, m);

      // random T matrix
      for (ordinal_type j=0;j<AA.NumCols();++j) {
        for (ordinal_type i=0;i<AA.NumRows();++i)
          TT.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        TT.Value(j,j) = abs(TT.Value(j,j));
      }

      {
        Teuchos::BLAS<ordinal_type,value_type> blas;

        // should be square
        const ordinal_type nn = AA.NumRows();
        const ordinal_type kk = TT.NumRows();

        blas.HERK(Teuchos::UPPER_TRI, Teuchos::CONJ_TRANS,
                  nn, kk,
                  1.0,
                  TT.ValuePtr(), TT.ColStride(),
                  0.0,
                  AA.ValuePtr(), AA.ColStride());
        AB.copy(AA);
      }

      const double flop = get_flop_chol<value_type>(m);

#ifdef HAVE_SHYLUTACHO_MKL
      mkl_set_num_threads(1);
#endif
      int ierr = 0;
      if (check) {
        timer.reset();
        DenseTaskViewType A(&AB);
        ierr = Chol<Uplo::Upper,AlgoChol::ExternalLapack>::invoke
          (TaskFactoryType::Policy(),
           TaskFactoryType::Policy().member_single(),
           A);
        t = timer.seconds();
        if (ierr) 
          ERROR(">> Fail to perform Cholesky (serial) : no reference solution -> no numeric error information");
        cout << "DenseCholByBlocks:: Serial Performance = " << (flop/t/1.0e9) << " [GFLOPs]" << endl;
      }

      {
        DenseHierMatrixBaseType HA;
        DenseMatrixHelper::flat2hier(AA, HA, mb, mb);

        DenseHierTaskViewType TA(&HA);
        timer.reset();
        auto future = TaskFactoryType::Policy().create_team
          (Chol<Uplo::Upper,AlgoChol::DenseByBlocks>
           ::TaskFunctor<DenseHierTaskViewType>(TA), 0);
        TaskFactoryType::Policy().spawn(future);
        Kokkos::Experimental::wait(TaskFactoryType::Policy());
        t = timer.seconds();
        cout << "DenseCholByBlocks:: Parallel Performance = " << (flop/t/1.0e9) << " [GFLOPs]" << endl;
      }

      if (!ierr && check) {
        typedef typename Teuchos::ScalarTraits<value_type>::magnitudeType real_type;
        real_type err = 0.0, norm = 0.0;
        for (ordinal_type j=0;j<AA.NumCols();++j)
          for (ordinal_type i=0;i<=j;++i) {
            const real_type diff = abs(AA.Value(i,j) - AB.Value(i,j));
            const real_type val  = AB.Value(i,j);
            err  += diff*diff;
            norm += val*val;
          }
        cout << "DenseCholByBlocks:: Check result::err = " << sqrt(err) << ", norm = " << sqrt(norm) << endl;
      }
    }

    return r_val;
  }
}

#endif
