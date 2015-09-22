#ifndef PANZER_KOKKOS_NODE_TYPE_HPP
#define PANZER_KOKKOS_NODE_TYPE_HPP

//#include "Kokkos_DefaultNode.hpp"
//#include "Tpetra_Map.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"

namespace panzer {

  // New Tpetra default node type
  // typedef Tpetra::Map<>::node_type TpetraNodeType;
  typedef typename Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device> TpetraNodeType;

}

#endif
