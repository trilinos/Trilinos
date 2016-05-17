#ifndef __TACHO_EXAMPLE_DENSE_CHOL_BY_BLOCKS_HPP__
#define __TACHO_EXAMPLE_DENSE_CHOL_BY_BLOCKS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"

#include "Tacho_DenseMatrixBase.hpp"
#include "Tacho_DenseMatrixView.hpp"
#include "Tacho_DenseFlopCount.hpp"


#include "Tacho_DenseMatrixTools.hpp"

#include "Tacho_TaskView.hpp"
#include "Tacho_TaskFactory.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Herk.hpp"
#include "Tacho_Trsm.hpp"
#include "Tacho_Chol.hpp"

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif

namespace Tacho {

  template<int ArgUplo,
           int ArgVariant,
           typename DeviceSpaceType>
  int exampleDenseCholByBlocks(const ordinal_type mmin,
                               const ordinal_type mmax,
                               const ordinal_type minc,
                               const ordinal_type mb,
                               const int max_concurrency,
                               const int max_task_dependence,
                               const int team_size,
                               const int mkl_nthreads,
                               const bool check,
                               const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);

    typedef Kokkos::Experimental::TaskPolicy<DeviceSpaceType> PolicyType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> DenseMatrixBaseHostType;
    typedef DenseMatrixView<DenseMatrixBaseHostType> DenseMatrixViewHostType;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,DeviceSpaceType> DenseMatrixBaseDeviceType;
    typedef DenseMatrixView<DenseMatrixBaseDeviceType> DenseMatrixViewDeviceType;
    typedef TaskView<DenseMatrixViewDeviceType> DenseTaskViewDeviceType;

    typedef DenseMatrixBase<DenseTaskViewDeviceType,ordinal_type,size_type,DeviceSpaceType> DenseHierMatrixBaseDeviceType;

    typedef DenseMatrixView<DenseHierMatrixBaseDeviceType> DenseHierMatrixViewDeviceType;
    typedef TaskView<DenseHierMatrixViewDeviceType> DenseHierTaskViewDeviceType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    double t = 0.0;

    std::cout << "DenseCholByBlocks:: test matrices "
              <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc 
              << " , mb = " << mb << std::endl;

    const size_t max_task_size = (3*sizeof(DenseTaskViewDeviceType)+sizeof(PolicyType)+128); 
    PolicyType policy(max_concurrency,
                      max_task_size,
                      max_task_dependence,
                      team_size);

    std::ostringstream os;
    os.precision(3);
    os << std::scientific;

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      os.str("");
      
      // host matrices
      DenseMatrixBaseHostType AA_host("AA_host", m, m), AB_host("AB_host"), TT_host("TT_host");

      // random T matrix
      {
        TT_host.createConfTo(AA_host);
        for (ordinal_type j=0;j<TT_host.NumCols();++j) {
          for (ordinal_type i=0;i<TT_host.NumRows();++i)
            TT_host.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
          TT_host.Value(j,j) = std::fabs(TT_host.Value(j,j));
        }
      }
      // create SPD matrix
      {
        Teuchos::BLAS<ordinal_type,value_type> blas;

        blas.HERK(ArgUplo == Uplo::Upper ? Teuchos::UPPER_TRI : Teuchos::LOWER_TRI, 
                  Teuchos::CONJ_TRANS,
                  m, m,
                  1.0,
                  TT_host.ValuePtr(), TT_host.ColStride(),
                  0.0,
                  AA_host.ValuePtr(), AA_host.ColStride());

        // preserve a copy of A
        AB_host.createConfTo(AA_host);        
        DenseMatrixTools::copy(AB_host, AA_host);
      }

      const double flop = DenseFlopCount<value_type>::Chol(m);

#ifdef HAVE_SHYLUTACHO_MKL
      mkl_set_num_threads(mkl_nthreads);
#endif
      os << "DenseCholByBlocks:: m = " << m << "  ";
      int ierr = 0;
      if (check) {
        timer.reset();
        DenseMatrixViewHostType A_host(AB_host); 
        ierr = Chol<ArgUplo,AlgoChol::ExternalLapack,Variant::One>::invoke
          (policy, policy.member_single(),
           A_host);
        t = timer.seconds();
        TACHO_TEST_FOR_ABORT( ierr, "Fail to perform Cholesky (serial)");
        os << ":: Serial Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      }

      DenseMatrixBaseDeviceType AA_device("AA_device");
      {
        timer.reset();
        AA_device.mirror(AA_host);
        t = timer.seconds();
        os << ":: Mirror = " << t << " [sec]  ";
      }

      {
        DenseHierMatrixBaseDeviceType HA_device("HA_device");
        
        DenseMatrixTools::createHierMatrix(HA_device, AA_device, mb, mb);

        DenseHierTaskViewDeviceType TA_device(HA_device);

        timer.reset();
        auto future = policy.proc_create_team
          (Chol<ArgUplo,AlgoChol::DenseByBlocks,ArgVariant>
           ::createTaskFunctor(policy, TA_device), 
           0);
        policy.spawn(future);
        Kokkos::Experimental::wait(policy);
        t = timer.seconds();
        os << ":: Parallel Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      }

      AA_host.mirror(AA_device);
      if (!ierr && check) {
        double err = 0.0, norm = 0.0;
        for (ordinal_type j=0;j<AA_host.NumCols();++j)
          for (ordinal_type i=0;i<=j;++i) {
            const double diff = abs(AA_host.Value(i,j) - AB_host.Value(i,j));
            const double val  = AB_host.Value(i,j);
            err  += diff*diff;
            norm += val*val;
          }
        os << ":: Check result ::norm = " << sqrt(norm) << ", error = " << sqrt(err);
      }
      std::cout << os.str() << std::endl;
    }
    
    return r_val;
  }
}

#endif
