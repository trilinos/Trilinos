#ifndef __TACHO_EXAMPLE_CRS_MATRIX_BASE_HPP__
#define __TACHO_EXAMPLE_CRS_MATRIX_BASE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixTools.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleCrsMatrixBase(const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace:: "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::   ";   HostSpaceType::print_configuration(std::cout, detail);

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>   CrsMatrixBaseHostType;
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,DeviceSpaceType> CrsMatrixBaseDeviceType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    
    // set user test matrix
    const ordinal_type m = 4;
    const size_type nnz = 12;
    size_type    ap_begin[m] = { 0, 3, 7, 10  };
    size_type    ap_end  [m] = { 3, 7, 10, 12 };
    ordinal_type aj[nnz] = { 0, 1, 3,
                             0, 1, 2, 3, 
                             1, 2, 3,
                             2, 3 };
    
    value_type   ax[nnz];
    for (size_type k=0;k<nnz;++k)
      ax[k] = 2.0*((value_type)std::rand()/(RAND_MAX)) - 1.0;      

    // wrap user data with kokkos view
    Kokkos::View<size_type*,   HostSpaceType,Kokkos::MemoryUnmanaged> ap_begin_ext(ap_begin, m);
    Kokkos::View<size_type*,   HostSpaceType,Kokkos::MemoryUnmanaged> ap_end_ext(ap_end, m);
    Kokkos::View<ordinal_type*,HostSpaceType,Kokkos::MemoryUnmanaged> aj_ext(aj, nnz);
    Kokkos::View<value_type*,  HostSpaceType,Kokkos::MemoryUnmanaged> ax_ext(ax, nnz);
    
    std::cout << "CrsMatrixBase:: test matrices "
              <<":: m = " << m << " , n = " << m << " , nnz = " << nnz << std::endl;

    {
      CrsMatrixBaseHostType AA("AA", m, m, nnz,
                               ap_begin_ext, ap_end_ext,
                               aj_ext, ax_ext);
      AA.showMe(std::cout) << std::endl;

      CrsMatrixBaseDeviceType BB("BB");
      timer.reset();
      BB.mirror(AA);
      double t_mirror = timer.seconds();
      
      //BB.showMe(std::cout) << std::endl;

      CrsMatrixBaseDeviceType CC("CC");
      CC.createConfTo(BB);

      timer.reset();
      CrsMatrixTools::copy(CC, BB);
      double t_copy = timer.seconds();

      //CC.showMe(std::cout) << std::endl;

      // check
      CrsMatrixBaseHostType RR("RR");
      RR.createConfTo(CC);
      RR.mirror(CC);

      RR.showMe(std::cout) << std::endl;

      double err = 0.0;
      for (ordinal_type k=0;k<nnz;++k) 
        err += std::fabs(AA.Value(k) - RR.Value(k));
      
      {
        const auto prec = std::cout.precision();
        std::cout.precision(4);
        
        std::cout << std::scientific
                  << "CrsMatrixBase:: dimension = " << m << " x " << m << ", " << " nnz = " << nnz << ", " 
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
