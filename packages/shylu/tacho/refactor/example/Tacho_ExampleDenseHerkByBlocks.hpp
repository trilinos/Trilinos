#ifndef __TACHO_EXAMPLE_DENSE_HERK_BY_BLOCKS_HPP__
#define __TACHO_EXAMPLE_DENSE_HERK_BY_BLOCKS_HPP__

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

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif   

namespace Tacho {

  template<int ArgUplo,
           int ArgTrans,
           int ArgVariant,
           typename DeviceSpaceType>
  int exampleDenseHerkByBlocks(const ordinal_type mmin,
                               const ordinal_type mmax,
                               const ordinal_type minc,
                               const ordinal_type k,
                               const ordinal_type mb,
                               const int max_concurrency,
                               const int memory_pool_grain_size,
                               const int mkl_nthreads,
                               const bool check,
                               const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
    std::cout << std::endl;

    typedef Kokkos::TaskPolicy<DeviceSpaceType> PolicyType;

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

    std::cout << "DenseHerkByBlocks:: test matrices "
              <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc 
              << " , k = "<< k << " , mb = " << mb << std::endl;

    const size_t max_task_size = (3*sizeof(DenseTaskViewDeviceType)+sizeof(PolicyType)+128); 

    PolicyType policy( typename PolicyType::memory_space(),
                       max_task_size*max_concurrency,
                       memory_pool_grain_size );

    typename Kokkos::TaskPolicy<Kokkos::Serial>::member_type serial_member;

    std::ostringstream os;
    os.precision(3);
    os << std::scientific;

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      os.str("");
      
      // host matrices
      DenseMatrixBaseHostType AA_host, CC_host("CC_host", m, m), CB_host;
      {
        if (ArgTrans == Trans::NoTranspose) 
          AA_host = DenseMatrixBaseHostType("AA_host", m, k); 
        else 
          AA_host = DenseMatrixBaseHostType("AA_host", k, m);
        
        for (ordinal_type j=0;j<AA_host.NumCols();++j)
          for (ordinal_type i=0;i<AA_host.NumRows();++i)
            AA_host.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        
        for (ordinal_type j=0;j<CC_host.NumCols();++j)
          for (ordinal_type i=0;i<CC_host.NumRows();++i)
            CC_host.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        
        CB_host.createConfTo(CC_host);
        DenseMatrixTools::copy(CB_host, CC_host);
      }

      const double flop = DenseFlopCount<value_type>::Syrk(k, m);

#ifdef HAVE_SHYLUTACHO_MKL
      mkl_set_num_threads(mkl_nthreads);
#endif

      os << "DenseHerkByBlocks:: m = " << m << " k = " << k << "  ";
      if (check) {
        timer.reset();
        DenseMatrixViewHostType A_host(AA_host), C_host(CB_host);
        Herk<ArgUplo,ArgTrans,AlgoHerk::ExternalBlas,Variant::One>::invoke
          (policy, serial_member,
           1.0, A_host, 1.0, C_host);
        t = timer.seconds();
        os << ":: Serial Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      }

      DenseMatrixBaseDeviceType AA_device("AA_device"), CC_device("CC_device");
      {
        timer.reset();
        AA_device.mirror(AA_host);
        CC_device.mirror(CC_host);
        t = timer.seconds();
        os << ":: Mirror = " << t << " [sec]  ";
      }

      {
        DenseHierMatrixBaseDeviceType HA_device("HA_device"), HC_device("HC_device");

        DenseMatrixTools::createHierMatrix(HA_device, AA_device, mb, mb);        
        DenseMatrixTools::createHierMatrix(HC_device, CC_device, mb, mb);        

        DenseHierTaskViewDeviceType TA_device(HA_device), TC_device(HC_device);
        timer.reset();
        auto future = policy.host_spawn
          (Herk<ArgUplo,ArgTrans,AlgoHerk::DenseByBlocks,ArgVariant>
           ::createTaskFunctor(policy, 1.0, TA_device, 1.0, TC_device));
        TACHO_TEST_FOR_EXCEPTION(future.is_null(), std::runtime_error,
                                 ">> host_spawn returns a null future");
        Kokkos::wait(policy);
        t = timer.seconds();       
        os << ":: Parallel Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      } 

      CC_host.mirror(CC_device);
      if (check) {
        double err = 0.0, norm = 0.0;
        switch (ArgUplo) {
        case Uplo::Upper: {
          for (ordinal_type j=0;j<CC_host.NumCols();++j)
            for (ordinal_type i=0;i<=j;++i) {
              const double diff = abs(CC_host.Value(i,j) - CB_host.Value(i,j));
              const double val  = CB_host.Value(i,j);
              err  += diff*diff;
              norm += val*val;
            }
          break;
        }
        case Uplo::Lower: {
          for (ordinal_type j=0;j<CC_host.NumCols();++j)
            for (ordinal_type i=j;i<CC_host.NumRows();++i) {
              const double diff = abs(CC_host.Value(i,j) - CB_host.Value(i,j));
              const double val  = CB_host.Value(i,j);
              err  += diff*diff;
              norm += val*val;
            }
          break;
        }
        }
        os << ":: Check result ::norm = " << sqrt(norm) << ", error = " << sqrt(err);
      }
      std::cout << os.str() << std::endl;
    }

    return r_val;
  }
}

#endif
