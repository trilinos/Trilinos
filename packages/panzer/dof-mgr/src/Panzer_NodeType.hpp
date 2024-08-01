// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_KOKKOS_NODE_TYPE_HPP
#define PANZER_KOKKOS_NODE_TYPE_HPP

//#include <Tpetra_KokkosCompat_DefaultNode.hpp>
//#include "Tpetra_Map.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace panzer {

  // New Tpetra default node type
  // typedef Tpetra::Map<>::node_type TpetraNodeType;
  typedef typename Tpetra::KokkosCompat::KokkosDeviceWrapperNode<PHX::Device> TpetraNodeType;

}

#endif
