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
  /**
   * \brief The Kokkos node type used to instantiate all Tpetra objects (Map, MultiVector, CrsGraph, etc.) throughout panzer.
   *
   * Tpetra classes are templated on a Kokkos node type in addition to
   * their scalar/ordinal types. This typedef pins that node type to a
   * `Tpetra::KokkosCompat::KokkosDeviceWrapperNode` wrapping panzer's
   * chosen Kokkos execution/memory space (`PHX::Device`, the same
   * device Phalanx evaluators run on), so that every Tpetra object
   * created by panzer executes on and shares the correct device. Use
   * this typedef wherever a Tpetra node type template argument is
   * needed, e.g. `Tpetra::Map<LocalOrdinal,GlobalOrdinal,TpetraNodeType>`.
   */
  typedef typename Tpetra::KokkosCompat::KokkosDeviceWrapperNode<PHX::Device> TpetraNodeType;

}

#endif
