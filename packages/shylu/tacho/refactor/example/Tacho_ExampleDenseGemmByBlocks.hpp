#ifndef __TACHO_EXAMPLE_DENSE_GEMM_BY_BLOCKS_HPP__
#define __TACHO_EXAMPLE_DENSE_GEMM_BY_BLOCKS_HPP__

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

#ifdef HAVE_SHYLUTACHO_MKL
#include "mkl_service.h"
#endif   

namespace Tacho {

  template<int ArgTransA,
           int ArgTransB,
           int ArgVariant,
           typename DeviceSpaceType>
  int exampleDenseGemmByBlocks(const ordinal_type mmin,
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

    std::cout << "DenseGemmByBlocks:: test matrices "
              <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc 
              << " , k = "<< k << " , mb = " << mb << std::endl;

    const size_t max_task_size = (3*sizeof(DenseTaskViewDeviceType)+sizeof(PolicyType)+128); 

    PolicyType policy( typename PolicyType::memory_space(), 
                       max_task_size*max_concurrency, 
                       memory_pool_grain_size );

    typename Kokkos::TaskPolicy<Kokkos::Serial>::member_type serial_member;

    std::cout << std::scientific;

    for (ordinal_type m=mmin;m<=mmax;m+=minc) {
      // host matrices
      DenseMatrixBaseHostType AA_host, BB_host, CC_host("CC_host", m, m), CB_host;
      {
        if (ArgTransA == Trans::NoTranspose) 
          AA_host = DenseMatrixBaseHostType("AA_host", m, k); 
        else 
          AA_host = DenseMatrixBaseHostType("AA_host", k, m);
        
        if (ArgTransB == Trans::NoTranspose) 
          BB_host = DenseMatrixBaseHostType("BB_host", k, m);
        else 
          BB_host = DenseMatrixBaseHostType("BB_host", m, k);
        
        for (ordinal_type j=0;j<AA_host.NumCols();++j)
          for (ordinal_type i=0;i<AA_host.NumRows();++i)
            AA_host.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        
        for (ordinal_type j=0;j<BB_host.NumCols();++j)
          for (ordinal_type i=0;i<BB_host.NumRows();++i)
            BB_host.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        
        for (ordinal_type j=0;j<CC_host.NumCols();++j)
          for (ordinal_type i=0;i<CC_host.NumRows();++i)
            CC_host.Value(i,j) = 2.0*((value_type)rand()/(RAND_MAX)) - 1.0;
        
        CB_host.createConfTo(CC_host);
        DenseMatrixTools::copy(CB_host, CC_host);
      }

      const double flop = DenseFlopCount<value_type>::Gemm(m, m, k);

#ifdef HAVE_SHYLUTACHO_MKL
      mkl_set_num_threads(mkl_nthreads);
#endif

      std::cout << "DenseGemmByBlocks:: m = " << m << " n = " << m << " k = " << k << "  ";
      if (check) {
        timer.reset();
        DenseMatrixViewHostType A_host(AA_host), B_host(BB_host), C_host(CB_host);
        Gemm<ArgTransA,ArgTransB,AlgoGemm::ExternalBlas,Variant::One>::invoke
          (policy, serial_member,
           1.0, A_host, B_host, 1.0, C_host);
        t = timer.seconds();
        std::cout << ":: Serial Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";

      }

      DenseMatrixBaseDeviceType AA_device("AA_device"), BB_device("BB_device"), CC_device("CC_device");
      {
        timer.reset();
        AA_device.mirror(AA_host);
        BB_device.mirror(BB_host);
        CC_device.mirror(CC_host);
        t = timer.seconds();
        std::cout << ":: Mirror = " << t << " [sec]  ";
      }

      {
        DenseHierMatrixBaseDeviceType HA_device("HA_device"), HB_device("HB_device"), HC_device("HC_device");

        DenseMatrixTools::createHierMatrix(HA_device, AA_device, mb, mb);        
        DenseMatrixTools::createHierMatrix(HB_device, BB_device, mb, mb);        
        DenseMatrixTools::createHierMatrix(HC_device, CC_device, mb, mb);        

        DenseHierTaskViewDeviceType TA_device(HA_device), TB_device(HB_device), TC_device(HC_device);
        timer.reset();

        auto future = policy.host_spawn
          (Gemm<ArgTransA,ArgTransB,AlgoGemm::DenseByBlocks,ArgVariant>
           ::createTaskFunctor(policy, 1.0, TA_device, TB_device, 1.0, TC_device));
        TACHO_TEST_FOR_EXCEPTION(future.is_null(), std::runtime_error, 
                                 ">> host_spawn returns a null future");
        Kokkos::wait(policy);

        t = timer.seconds();       
        std::cout << ":: Parallel Performance = " << (flop/t/1.0e9) << " [GFLOPs]  ";
      } 

      CC_host.mirror(CC_device);
      if (check) {
        double err = 0.0, norm = 0.0;
        for (ordinal_type j=0;j<CC_host.NumCols();++j)
          for (ordinal_type i=0;i<CC_host.NumRows();++i) {
            const double diff = abs(CC_host.Value(i,j) - CB_host.Value(i,j));
            const double val  = CB_host.Value(i,j);
            err  += diff*diff;
            norm += val*val;
          }
        std::cout << ":: Check result ::norm = " << sqrt(norm) << ", error = " << sqrt(err);
      }

      std::cout << std::endl;
    }

    return r_val;
  }
}

#endif
