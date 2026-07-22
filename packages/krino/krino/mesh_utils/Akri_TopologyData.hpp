#ifndef KRINO_KRINO_KRINO_LIB_AKRI_TOPOLOGYDATA_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_TOPOLOGYDATA_HPP_
#include <stk_topology/topology_detail/topology_data.hpp>

namespace krino {

template <stk::topology::topology_t TOPO>
struct TopologyData {
  static constexpr unsigned num_nodes() { return stk::topology_detail::topology_data<TOPO>::num_nodes; }
  static constexpr unsigned num_edges() { return stk::topology_detail::topology_data<TOPO>::num_edges; }
  static constexpr unsigned num_sides() { return (stk::topology_detail::topology_data<TOPO>::side_rank == stk::topology::FACE_RANK) ? stk::topology_detail::topology_data<TOPO>::num_faces : stk::topology_detail::topology_data<TOPO>::num_edges; }
  static constexpr unsigned spatial_dimension() { return stk::topology_detail::topology_data<TOPO>::spatial_dimension_vector[3] ? 3 : 2; }
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_TOPOLOGYDATA_HPP_ */
