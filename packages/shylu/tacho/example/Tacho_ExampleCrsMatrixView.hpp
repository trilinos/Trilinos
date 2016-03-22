#ifndef __TACHO_EXAMPLE_CRS_MATRIX_VIEW_HPP__
#define __TACHO_EXAMPLE_CRS_MATRIX_VIEW_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixView.hpp"
#include "Tacho_CrsRowView.hpp"
#include "Tacho_CrsMatrixTools.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleCrsMatrixView(const bool verbose) {
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
    
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>   CrsMatrixBaseHostType;
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,DeviceSpaceType> CrsMatrixBaseDeviceType;
    typedef CrsRowView<CrsMatrixBaseDeviceType> CrsRowViewDeviceType;
    typedef CrsMatrixView<CrsMatrixBaseDeviceType> CrsMatrixViewDeviceType;

    int r_val = 0;

    //Kokkos::Impl::Timer timer;
    
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
      BB.mirror(AA);
      
      BB.showMe(std::cout) << std::endl;

      CrsMatrixViewDeviceType B(BB), C(BB, 0, 3, 0, 3), D(BB, 1, 2, 1, 2);
      Kokkos::View<CrsRowViewDeviceType*,DeviceSpaceType> rows("RowViews", 
                                                               B.NumRows() + C.NumRows() + D.NumRows());
      B.showMe(std::cout) << std::endl;
      C.showMe(std::cout) << std::endl;
      D.showMe(std::cout) << std::endl;

      ordinal_type offset = 0;

      B.setRowViewArray(Kokkos::subview(rows, range_type(offset, offset + B.NumRows()))); 
      offset += B.NumRows();
      B.showMe(std::cout) << std::endl;

      C.setRowViewArray(Kokkos::subview(rows, range_type(offset, offset + C.NumRows()))); 
      offset += C.NumRows();
      C.showMe(std::cout) << std::endl;

      D.setRowViewArray(Kokkos::subview(rows, range_type(offset, offset + D.NumRows()))); 
      offset += D.NumRows();
      D.showMe(std::cout) << std::endl;
    }

    return r_val;
  }
}

#endif
