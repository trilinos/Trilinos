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

  template<typename Scalar,typename...Props>
  Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>
  getNonConstDynRankViewFromConstMDField(const PHX::MDField<const Scalar,Props...>& f) {

    using drv_type = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>;
    using nonconst_data_type = typename Sacado::ScalarType< typename drv_type::value_type >::type*;
    const int rank = f.rank();
    Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged> tmp;

#ifdef PHX_DEBUG
    TEUCHOS_ASSERT( (rank > 0) && (rank < 6) );
#endif

    if (Sacado::IsFad<Scalar>::value) {
      const int num_derivatives = Kokkos::dimension_scalar(f.get_static_view());
      if (rank==1)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),num_derivatives);
      else if (rank==2)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),num_derivatives);
      else if (rank==3)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),num_derivatives);
      else if (rank==4)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3),num_derivatives);
      else if (rank==5)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3),f.extent(4),num_derivatives);
    }
    else {
      if (rank==1)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0));
      else if (rank==2)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1));
      else if (rank==3)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2));
      else if (rank==4)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3));
      else if (rank==5)
        tmp = Kokkos::DynRankView<Scalar,typename PHX::DevLayout<Scalar>::type,Kokkos::MemoryUnmanaged>(const_cast<nonconst_data_type>(f.get_static_view().data()),f.extent(0),f.extent(1),f.extent(2),f.extent(3),f.extent(4));
    }

    return tmp;
  }

} // namespace PHX

#endif
