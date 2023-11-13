// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_CDFEM_Parent_Edge_h
#define Akri_CDFEM_Parent_Edge_h

#include <Akri_InterfaceID.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <functional>
#include <vector>
#include <map>

namespace stk { namespace mesh { class BulkData; } }

namespace krino {

class CDFEM_Parent_Edge {
public:
  CDFEM_Parent_Edge(const bool oneLSPerPhase,
    const std::vector<stk::mesh::Entity> & edgeNodes,
    const std::vector<double> & edgeNodePositions,
    const std::vector<std::vector<double> > & nodes_isovar)
  : my_edge_nodes(edgeNodes),
    my_edge_node_positions(edgeNodePositions)
  {
    STK_ThrowAssert(edgeNodePositions.size() == nodes_isovar.size());
    find_crossings(oneLSPerPhase, nodes_isovar);
  }

  CDFEM_Parent_Edge(const bool oneLSPerPhase,
      const std::vector<double> & edgeNodePositions,
      const std::vector<std::vector<double> > & nodes_isovar)
  : CDFEM_Parent_Edge(oneLSPerPhase, {}, edgeNodePositions, nodes_isovar) {}

  CDFEM_Parent_Edge(const bool oneLSPerPhase,
    const std::vector<stk::mesh::Entity> & edgeNodes,
    const std::vector<std::vector<double> > & nodes_isovar)
  : CDFEM_Parent_Edge(oneLSPerPhase, edgeNodes, {0.,1.}, nodes_isovar) {}

  CDFEM_Parent_Edge(const bool oneLSPerPhase, const std::vector<std::vector<double> > & nodes_isovar)
  : CDFEM_Parent_Edge(oneLSPerPhase, {0.,1.}, nodes_isovar) {}

  CDFEM_Parent_Edge(const std::vector<stk::mesh::Entity> & edgeNodes,
    const std::vector<double> & edgeNodePositions);
  CDFEM_Parent_Edge() {}

  // Must be larger than machine epsilon, but puts limit on smallest effective snap tolerance
  static double MinSize();

  bool valid() const { return !my_edge_node_positions.empty(); }

  void find_crossings(const bool oneLSPerPhase, const std::vector<std::vector<double> > & nodes_isovar);
  void collapse_small_segments_while_preserving_topology(const double snapTol);

  void adjust_crossing_locations_based_on_node_captured_domains(const bool oneLSPerPhase, const std::vector<int> & sortedParentNode0Domains, const std::vector<int> & sortedParentNode1Domains);

  bool have_crossing(const InterfaceID key) const { return my_crossings.find(key) != my_crossings.end(); }
  const CrossingMap & get_crossings() const { return my_crossings; }
  const CrossingMap & get_crossings_including_fake() const { return my_crossings_including_fake; }
  bool all_fake_crossings_are_really_fake() const;

  double get_crossing_position(const InterfaceID key) const {
    CrossingMap::const_iterator it = my_crossings.find(key);
    return it != my_crossings.end() ? it->second : -1.0;
  }
  // Crossing sign is defined as whether the interface is + or - for x > crossing point
  // For multiple LS problems we'll say that -1 corresponds to InterfaceID.first being lower for x > x_crossing_point
  // For uncrossed edges it will be the sign of both parent nodes
  int get_crossing_sign(const InterfaceID key) const {
    CrossingSignMap::const_iterator it = my_crossing_signs.find(key);
    return it->second;
  }
  std::tuple<double, int, bool> get_crossing_position_and_sign(const bool oneLSPerPhase, const InterfaceID key) const;
  unsigned get_num_nodes() const { return my_edge_node_positions.size(); }
  const std::vector<stk::mesh::Entity> & get_nodes() const { return my_edge_nodes; }
  bool have_any_crossings() const;
  std::pair<stk::mesh::Entity,stk::mesh::Entity> get_parent_nodes() const { return std::pair<stk::mesh::Entity,stk::mesh::Entity>(my_edge_nodes.front(), my_edge_nodes.back()); }
  int get_uncrossed_phase() const { return (edge_phases.size() == 1) ? (*edge_phases.begin()) : -1; }
  const std::set<int> & get_edge_phases() const { return edge_phases; }
  double get_edge_node_position(stk::mesh::Entity edgeNode) const;

  friend std::ostream & operator << (std::ostream &os, const CDFEM_Parent_Edge & edge);
  void debug_print_crossings() const;

private:
  void find_crossings_multiple_levelset(const std::vector<std::vector<double> > & nodes_isovar);
  std::pair<double, int> find_crossing_position_and_sign(const bool oneLSPerPhase, const InterfaceID key, const std::vector<std::vector<double> > & nodes_isovar) const;
  void find_crossings_including_fake_ones(const bool oneLSPerPhase, const std::vector<std::vector<double> > & nodes_isovar);
  double get_fake_crossing_position(const InterfaceID key) const {
    CrossingMap::const_iterator it = my_crossings_including_fake.find(key);
    return it != my_crossings_including_fake.end() ? it->second : -1.0;
  }
  int get_fake_crossing_sign(const InterfaceID key) const {
    CrossingSignMap::const_iterator it = my_crossing_signs_including_fake.find(key);
    return it->second;
  }
  void fixup_fake_crossing_locations_for_consistency();

  std::vector<stk::mesh::Entity> my_edge_nodes;
  std::vector<double> my_edge_node_positions;
  CrossingMap my_crossings;
  CrossingSignMap my_crossing_signs;
  CrossingMap my_crossings_including_fake;
  CrossingSignMap my_crossing_signs_including_fake;
  std::set<int> edge_phases;
};

inline std::ostream & operator << (std::ostream &os, const CDFEM_Parent_Edge & edge)
{
  const unsigned num_nodes = edge.get_num_nodes();
  os << "CDFEM parent edge has " << num_nodes << " nodes and phases { ";
  for (int phase : edge.get_edge_phases())
    os << phase << " ";
  os << "}\n  crossings: { ";
  const auto oldPrecision = os.precision();
  os.precision(16);
  for ( auto && crossing : edge.get_crossings() )
  {
    os << crossing.first << "@" << crossing.second << ", sign=" << edge.get_crossing_sign(crossing.first) << " ";
  }
  os << "}" << "\n  crossings including fake: { ";
  for ( auto && crossing : edge.get_crossings_including_fake() )
  {
    os << crossing.first << "@" << crossing.second << ", sign=" << edge.get_fake_crossing_sign(crossing.first) << " ";
  }
  os << "}" << "\n";
  os.precision(oldPrecision);
  return os;
}

} // namespace krino

#endif // Akri_CDFEM_Parent_Edge_h
