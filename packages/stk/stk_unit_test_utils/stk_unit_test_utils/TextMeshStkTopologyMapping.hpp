#ifndef TextMeshStkTopologyMapping_hpp
#define TextMeshStkTopologyMapping_hpp

#include "TextMeshDataTypes.hpp"

#include "stk_mesh/base/Types.hpp"                   // for EntityId, etc
#include "stk_topology/topology.hpp"                 // for topology, etc

#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector

namespace stk
{
namespace unit_test_util
{
struct StkTopologyMapEntry {
  using Ordinal = stk::mesh::ConnectivityOrdinal;
  using Permutation = unsigned;

  static constexpr Ordinal InvalidOrdinal = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
  static constexpr Permutation InvalidPermutation = stk::mesh::INVALID_PERMUTATION;

  stk::topology topology;

  StkTopologyMapEntry() : topology(stk::topology::INVALID_TOPOLOGY) {}

  StkTopologyMapEntry(stk::topology topoEntry) : topology(topoEntry) {}

  StkTopologyMapEntry(const StkTopologyMapEntry &topo) : topology(topo.topology) {}

  bool defined_on_spatial_dimension(const unsigned spatialDim) const
  {
    if (spatialDim > 3) {
      return false;
    }
    return topology.defined_on_spatial_dimension(spatialDim);
  }

  const std::string name() const { return topology.name(); }

  int num_nodes() const { return topology.num_nodes(); }

  bool operator==(const StkTopologyMapEntry &rhs) const { return topology == rhs.topology; }

  bool operator!=(const StkTopologyMapEntry &rhs) const { return !(*this == rhs); }

  int num_sides() const { return topology.num_sides(); }

  bool valid_side(unsigned side) const
  {
    unsigned numSides = topology.num_sides();
    if (side > 0 && side <= numSides) return true;
    return false;
  }

  std::string side_topology_name(unsigned side) const
  {
    if (!valid_side(side)) return "";

    unsigned ordinal = side - 1;
    stk::topology sideTopology = topology.side_topology(ordinal);
    return sideTopology.name();
  }

  StkTopologyMapEntry side_topology(unsigned side) const
  {
    if (!valid_side(side)) return StkTopologyMapEntry();

    unsigned ordinal = side - 1;
    stk::topology sideTopology = topology.side_topology(ordinal);
    return StkTopologyMapEntry(sideTopology);
  }

  unsigned side_topology_num_nodes(unsigned side) const
  {
    if (!valid_side(side)) return 0;

    unsigned ordinal = side - 1;
    stk::topology sideTopology = topology.side_topology(ordinal);
    return sideTopology.num_nodes();
  }

  std::vector<Ordinal> side_topology_node_indices(unsigned side) const
  {
    if (!valid_side(side)) return std::vector<Ordinal>();

    unsigned sideOrdinal = side - 1;
    stk::topology sideTopology = topology.side_topology(sideOrdinal);
    std::vector<Ordinal> elementNodeOrdinalVector(sideTopology.num_nodes());
    topology.side_node_ordinals(sideOrdinal, elementNodeOrdinalVector.data());
    return elementNodeOrdinalVector;
  }

  bool is_shell() const { return topology.is_shell(); }

  unsigned num_permutations() const { return topology.num_permutations(); }

  unsigned num_positive_permutations() const { return topology.num_positive_permutations(); }

  bool is_positive_polarity(Permutation permutation) const { return topology.is_positive_polarity(permutation); }

  bool valid_permutation(Permutation permutation) const { return permutation < num_permutations(); }

  bool fill_permutation_indices(Permutation permutation, std::vector<Ordinal>& nodeOrdinalVector) const
  {
    if (!valid_permutation(permutation)) return false;

    nodeOrdinalVector.resize(topology.num_nodes());
    topology.permutation_node_ordinals(permutation, nodeOrdinalVector.data());
    return true;
  }

  std::vector<Ordinal> permutation_indices(Permutation permutation) const
  {
    std::vector<Ordinal> nodeOrdinalVector;
    fill_permutation_indices(permutation, nodeOrdinalVector);
    return nodeOrdinalVector;
  }
};

inline std::ostream& operator<<(std::ostream& out, const StkTopologyMapEntry& t)
{
  return out << t.name();
}

class StkTopologyMapping : public text_mesh::TopologyMapping<StkTopologyMapEntry>
{
 public:
  StkTopologyMapEntry invalid_topology() const override { return StkTopologyMapEntry(); }

  void initialize_topology_map() override
  {
    m_nameToTopology = {
        {"NODE", StkTopologyMapEntry(stk::topology::NODE)},
        {"LINE_2", StkTopologyMapEntry(stk::topology::LINE_2)},
        {"LINE_3", StkTopologyMapEntry(stk::topology::LINE_3)},
        {"TRI_3", StkTopologyMapEntry(stk::topology::TRI_3)},
        {"TRIANGLE_3", StkTopologyMapEntry(stk::topology::TRI_3)},
        {"TRI_4", StkTopologyMapEntry(stk::topology::TRI_4)},
        {"TRIANGLE_4", StkTopologyMapEntry(stk::topology::TRI_4)},
        {"TRI_6", StkTopologyMapEntry(stk::topology::TRI_6)},
        {"TRIANGLE_6", StkTopologyMapEntry(stk::topology::TRI_6)},
        {"QUAD_4", StkTopologyMapEntry(stk::topology::QUAD_4)},
        {"QUADRILATERAL_4", StkTopologyMapEntry(stk::topology::QUAD_4)},
        {"QUAD_6", StkTopologyMapEntry(stk::topology::QUAD_6)},
        {"QUADRILATERAL_6", StkTopologyMapEntry(stk::topology::QUAD_6)},
        {"QUAD_8", StkTopologyMapEntry(stk::topology::QUAD_8)},
        {"QUADRILATERAL_8", StkTopologyMapEntry(stk::topology::QUAD_8)},
        {"QUAD_9", StkTopologyMapEntry(stk::topology::QUAD_9)},
        {"QUADRILATERAL_9", StkTopologyMapEntry(stk::topology::QUAD_9)},
        {"PARTICLE", StkTopologyMapEntry(stk::topology::PARTICLE)},
        {"LINE_2_1D", StkTopologyMapEntry(stk::topology::LINE_2_1D)},
        {"LINE_3_1D", StkTopologyMapEntry(stk::topology::LINE_3_1D)},
        {"BEAM_2", StkTopologyMapEntry(stk::topology::BEAM_2)},
        {"BEAM_3", StkTopologyMapEntry(stk::topology::BEAM_3)},
        {"SHELL_LINE_2", StkTopologyMapEntry(stk::topology::SHELL_LINE_2)},
        {"SHELL_LINE_3", StkTopologyMapEntry(stk::topology::SHELL_LINE_3)},
        {"SPRING_2", StkTopologyMapEntry(stk::topology::SPRING_2)},
        {"SPRING_3", StkTopologyMapEntry(stk::topology::SPRING_3)},
        {"TRI_3_2D", StkTopologyMapEntry(stk::topology::TRI_3_2D)},
        {"TRIANGLE_3_2D", StkTopologyMapEntry(stk::topology::TRI_3_2D)},
        {"TRI_4_2D", StkTopologyMapEntry(stk::topology::TRI_4_2D)},
        {"TRIANGLE_4_2D", StkTopologyMapEntry(stk::topology::TRI_4_2D)},
        {"TRI_6_2D", StkTopologyMapEntry(stk::topology::TRI_6_2D)},
        {"TRIANGLE_6_2D", StkTopologyMapEntry(stk::topology::TRI_6_2D)},
        {"QUAD_4_2D", StkTopologyMapEntry(stk::topology::QUAD_4_2D)},
        {"QUADRILATERAL_4_2D", StkTopologyMapEntry(stk::topology::QUAD_4_2D)},
        {"QUAD_8_2D", StkTopologyMapEntry(stk::topology::QUAD_8_2D)},
        {"QUADRILATERAL_8_2D", StkTopologyMapEntry(stk::topology::QUAD_8_2D)},
        {"QUAD_9_2D", StkTopologyMapEntry(stk::topology::QUAD_9_2D)},
        {"QUADRILATERAL_9_2D", StkTopologyMapEntry(stk::topology::QUAD_9_2D)},
        {"SHELL_TRI_3", StkTopologyMapEntry(stk::topology::SHELL_TRI_3)},
        {"SHELL_TRIANGLE_3", StkTopologyMapEntry(stk::topology::SHELL_TRI_3)},
        {"SHELL_TRI_4", StkTopologyMapEntry(stk::topology::SHELL_TRI_4)},
        {"SHELL_TRIANGLE_4", StkTopologyMapEntry(stk::topology::SHELL_TRI_4)},
        {"SHELL_TRI_6", StkTopologyMapEntry(stk::topology::SHELL_TRI_6)},
        {"SHELL_TRIANGLE_6", StkTopologyMapEntry(stk::topology::SHELL_TRI_6)},
        {"SHELL_QUAD_4", StkTopologyMapEntry(stk::topology::SHELL_QUAD_4)},
        {"SHELL_QUADRILATERAL_4", StkTopologyMapEntry(stk::topology::SHELL_QUAD_4)},
        {"SHELL_QUAD_8", StkTopologyMapEntry(stk::topology::SHELL_QUAD_8)},
        {"SHELL_QUADRILATERAL_8", StkTopologyMapEntry(stk::topology::SHELL_QUAD_8)},
        {"SHELL_QUAD_9", StkTopologyMapEntry(stk::topology::SHELL_QUAD_9)},
        {"SHELL_QUADRILATERAL_9", StkTopologyMapEntry(stk::topology::SHELL_QUAD_9)},
        {"TET_4", StkTopologyMapEntry(stk::topology::TET_4)},
        {"TETRAHEDRON_4", StkTopologyMapEntry(stk::topology::TET_4)},
        {"TET_8", StkTopologyMapEntry(stk::topology::TET_8)},
        {"TETRAHEDRON_8", StkTopologyMapEntry(stk::topology::TET_8)},
        {"TET_10", StkTopologyMapEntry(stk::topology::TET_10)},
        {"TETRAHEDRON_10", StkTopologyMapEntry(stk::topology::TET_10)},
        {"TET_11", StkTopologyMapEntry(stk::topology::TET_11)},
        {"TETRAHEDRON_11", StkTopologyMapEntry(stk::topology::TET_11)},
        {"PYRAMID_5", StkTopologyMapEntry(stk::topology::PYRAMID_5)},
        {"PYRAMID_13", StkTopologyMapEntry(stk::topology::PYRAMID_13)},
        {"PYRAMID_14", StkTopologyMapEntry(stk::topology::PYRAMID_14)},
        {"WEDGE_6", StkTopologyMapEntry(stk::topology::WEDGE_6)},
        {"WEDGE_12", StkTopologyMapEntry(stk::topology::WEDGE_12)},
        {"WEDGE_15", StkTopologyMapEntry(stk::topology::WEDGE_15)},
        {"WEDGE_18", StkTopologyMapEntry(stk::topology::WEDGE_18)},
        {"HEX_8", StkTopologyMapEntry(stk::topology::HEX_8)},
        {"HEXAHEDRON_8", StkTopologyMapEntry(stk::topology::HEX_8)},
        {"HEX_20", StkTopologyMapEntry(stk::topology::HEX_20)},
        {"HEXAHEDRON_20", StkTopologyMapEntry(stk::topology::HEX_20)},
        {"HEX_27", StkTopologyMapEntry(stk::topology::HEX_27)},
        {"HEXAHEDRON_27", StkTopologyMapEntry(stk::topology::HEX_27)}
    };
  }
};

namespace simple_fields {

struct STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
StkTopologyMapEntry : public stk::unit_test_util::StkTopologyMapEntry {};

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
StkTopologyMapping : public stk::unit_test_util::StkTopologyMapping {};

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk

#endif
