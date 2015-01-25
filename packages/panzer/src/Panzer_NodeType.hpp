#ifndef PANZER_KOKKOS_NODE_TYPE_HPP
#define PANZER_KOKKOS_NODE_TYPE_HPP

//#include "Kokkos_DefaultNode.hpp"
#include "Tpetra_Map.hpp"

namespace panzer {

  // Should really default off the phalanx node type to be constent
  // with assembly.  Can't do this until kokkos refactor branch merged
  // back in.

  // Old default node type
  //typedef KokkosClassic::DefaultNode::DefaultNodeType TpetraNodeType;

  // New Tpetra default node type
  typedef Tpetra::Map<>::node_type TpetraNodeType;

}

#endif
