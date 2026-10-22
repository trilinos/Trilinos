#ifndef KRINO_KRINO_KRINO_LIB_AKRI_TOPOLOGYDATA_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_TOPOLOGYDATA_HPP_
#include <stk_topology/topology.hpp>
#include <stk_topology/topology_detail/topology_data.hpp>

namespace krino {

template <stk::topology::topology_t TOPO>
struct TopologyData {
  static constexpr unsigned num_nodes() { return stk::topology_detail::topology_data<TOPO>::num_nodes; }
  static constexpr unsigned num_edges() { return stk::topology_detail::topology_data<TOPO>::num_edges; }
  static constexpr unsigned num_sides() { return (stk::topology_detail::topology_data<TOPO>::side_rank == stk::topology::FACE_RANK) ? stk::topology_detail::topology_data<TOPO>::num_faces : stk::topology_detail::topology_data<TOPO>::num_edges; }
  static constexpr unsigned spatial_dimension() { return stk::topology_detail::topology_data<TOPO>::spatial_dimension_vector[3] ? 3 : 2; }
};

// Topology of the interface (side) that separates the phases of a decomposed element of the
// given topology. For non-shell topologies this is the element's side topology; shells are a
// special case whose interface is a beam. Throws for an unrecognized shell topology.
stk::topology interface_topology(const stk::topology elementTopology);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_TOPOLOGYDATA_HPP_ */
