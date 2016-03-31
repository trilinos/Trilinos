#ifndef __TACHO_EXAMPLE_DENSE_MATRIX_BASE_HPP__
#define __TACHO_EXAMPLE_DENSE_MATRIX_BASE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_DenseMatrixBase.hpp"
#include "Tacho_DenseMatrixTools.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleDenseMatrixBase(const ordinal_type mmin,
                             const ordinal_type mmax,
                             const ordinal_type minc,
                             const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
  
    typedef DenseMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>   DenseMatrixBaseHostType;
    typedef DenseMatrixBase<value_type,ordinal_type,size_type,DeviceSpaceType> DenseMatrixBaseDeviceType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;

    std::cout << "DenseMatrixBase:: test matrices "
              <<":: mmin = " << mmin << " , mmax = " << mmax << " , minc = " << minc << std::endl;
    
    for (auto m=mmin;m<=mmax;m+=minc) {
      // random test matrix on host 
      DenseMatrixBaseHostType TT("TT", m, m);
      for (ordinal_type j=0;j<TT.NumCols();++j) {
        for (ordinal_type i=0;i<TT.NumRows();++i)
          TT.Value(i,j) = 2.0*((value_type)std::rand()/(RAND_MAX)) - 1.0;
        TT.Value(j,j) = std::fabs(TT.Value(j,j));
      }
      if (verbose)
        TT.showMe(std::cout) << std::endl;

      DenseMatrixBaseDeviceType AA("AA"); 

      timer.reset();
      AA.mirror(TT);
      double t_mirror = timer.seconds();

      DenseMatrixBaseDeviceType BB("BB");
      BB.createConfTo(AA);

      timer.reset();
      DenseMatrixTools::copy(BB, AA);
      double t_copy = timer.seconds();

      // check
      DenseMatrixBaseHostType RR("RR");
      RR.createConfTo(BB);
      RR.mirror(BB);
      if (verbose)
        RR.showMe(std::cout) << std::endl;

      double err = 0.0;
      for (ordinal_type j=0;j<TT.NumCols();++j) 
        for (ordinal_type i=0;i<TT.NumRows();++i)
          err += std::fabs(TT.Value(i,j) - RR.Value(i,j));
      
      {
        const auto prec = std::cout.precision();
        std::cout.precision(4);
        
        std::cout << std::scientific
                  << "DenseMatrixBase:: dimension = " << m << " x " << m << ", "
                  << "Mirroring to device  = " << t_mirror << " [sec], " 
                  << "Elementwise copy on device = " << t_copy << " [sec], " 
                  << "Error = " << err 
                  << std::endl;
        
        std::cout.unsetf(std::ios::scientific);
        std::cout.precision(prec);
      }
    }

    return r_val;
  }
}

#endif
