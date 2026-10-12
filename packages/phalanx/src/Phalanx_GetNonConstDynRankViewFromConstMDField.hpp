// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_GET_NON_CONST_DYN_RANK_VIEW_FROM_CONST_MDFIELD_HPP
#define PHALANX_GET_NON_CONST_DYN_RANK_VIEW_FROM_CONST_MDFIELD_HPP

#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Teuchos_Assert.hpp"

namespace PHX {

  /** \brief Wrap a const MDField's data in a non-const unmanaged DynRankView.
    *
    * The device comes from the field, not from PHX::Device: an MDField's device
    * is whatever its template arguments say, falling back to PHX::Device only
    * when the caller did not name one.  The view aliases the field's data, so
    * naming any other device would type the view for a memory space the
    * pointer does not live in.
    *
    * A device argument has to appear at all.  Without one the view takes
    * Kokkos' default device, which only coincides with PHX::Device while the
    * memory space is the execution space's own -- so with
    * Phalanx_ENABLE_SHARED_SPACE the view was typed in CudaSpace over a
    * CudaUVMSpace pointer.  Handing that to something expecting a view on
    * PHX::Device, Intrepid2::CellTools for instance, is a cross-space
    * conversion, and nvcc reports it as a cicc segfault rather than an error.
    */
  template<typename Scalar,typename...Props>
  Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>
  getNonConstDynRankViewFromConstMDField(const PHX::MDField<const Scalar,Props...>& f) {

    using drv_type = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>;
    using nonconst_data_type = typename Sacado::ScalarType< typename drv_type::value_type >::type*;
    const int rank = f.rank();
    Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged> tmp;

#ifdef PHX_DEBUG
    TEUCHOS_ASSERT( (rank > 0) && (rank < 6) );
#endif

    if (Sacado::IsFad<Scalar>::value) {
      const int num_derivatives = Sacado::dimension_scalar(f.get_static_view());
      if (rank==1)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),num_derivatives);
      else if (rank==2)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),num_derivatives);
      else if (rank==3)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),num_derivatives);
      else if (rank==4)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3),num_derivatives);
      else if (rank==5)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3),f.extent(4),num_derivatives);
    }
    else {
      if (rank==1)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0));
      else if (rank==2)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1));
      else if (rank==3)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2));
      else if (rank==4)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3));
      else if (rank==5)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,typename PHX::MDField<const Scalar,Props...>::device_type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3),f.extent(4));
    }

    return tmp;
  }

} // namespace PHX

#endif
