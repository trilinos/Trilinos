// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_VIEW_FACTORY_HPP
#define PANZER_VIEW_FACTORY_HPP

#include "Kokkos_ViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

  //! Wrapper to simplify Panzer use of Sacado ViewFactory 
  template<typename InputArray,typename ... DimensionPack>
  Kokkos::DynRankView<typename InputArray::value_type,PHX::Device>
  createDynRankView(const InputArray& a, const std::string& name, const DimensionPack... dims)
  {
    using return_type = Kokkos::DynRankView<typename InputArray::value_type,PHX::Device>;
    return PHX::ViewFactory<InputArray>::template create_view<return_type>(a,name,dims...);
  }

}

#endif
