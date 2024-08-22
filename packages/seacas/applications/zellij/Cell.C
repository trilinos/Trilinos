// Copyright(C) 2021, 2022, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <numeric>

#include "Cell.h"

#include <Ioss_NodeBlock.h>
#include <Ioss_SmartAssert.h>
#include <algorithm>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <string>

//! \file
extern unsigned int debug_level;

template <> struct fmt::formatter<Loc> : formatter<std::string>
{
  // parse is inherited from formatter<std::string>.
  template <typename FormatContext> auto format(Loc l, FormatContext &ctx) const
  {
    std::string name = "unknown";
    switch (l) {
    case Loc::C: name = "Center"; break;
    case Loc::BL: name = "Bottom Left"; break;
    case Loc::B: name = "Bottom"; break;
    case Loc::BR: name = "Bottom Right"; break;
    case Loc::L: name = "Left"; break;
    case Loc::R: name = "Right"; break;
    case Loc::TL: name = "Top Left"; break;
    case Loc::T: name = "Top"; break;
    case Loc::TR: name = "Top Right"; break;
    }
    return formatter<std::string>::format(name, ctx);
  }
};

namespace {
  // Iterate over the interior nodes on the specified face.  Skips
  // the corner nodes on the I-J intersections.  Processes I-K and
  // J-K corners.
  template <typename INT>
  void process_face_nodes(const std::vector<INT> &node_map, std::vector<INT> &nodes,
                          std::vector<INT> &procs, const std::vector<int64_t> &face_nodes,
                          size_t KK, int rank)
  {
    for (size_t i = KK; i < face_nodes.size() - KK; i++) {
      auto idx = face_nodes[i];
      nodes.push_back(node_map[idx + 1]);
      procs.push_back(rank);
    }
  }

  // Iterate over the specified corner nodes.
  template <typename INT>
  void process_corner_nodes(const std::vector<INT> &node_map, std::vector<INT> &nodes,
                            std::vector<INT> &procs, const std::vector<int64_t> &face_nodes,
                            size_t KK, int rank, Loc location)
  {
    if (location == Loc::BL || location == Loc::TL) {
      for (size_t i = 0; i < KK; i++) {
        auto idx = face_nodes[i];
        nodes.push_back(node_map[idx + 1]);
        procs.push_back(rank);
      }
    }
    else {
      for (size_t i = face_nodes.size() - KK; i < face_nodes.size(); i++) {
        auto idx = face_nodes[i];
        nodes.push_back(node_map[idx + 1]);
        procs.push_back(rank);
      }
    }
  }

  // Return a vector (possibly empty) of the ranks of the cells that surround this
  // cell.  Includes the rank that this cell is on (cell_ranks[0]).
  std::vector<int> get_shared_ranks(const std::array<int, 9> &cell_ranks)
  {
    std::vector<int> ranks(9);
    std::copy(cell_ranks.begin(), cell_ranks.end(), ranks.begin());

    // Set all `-1` (non-neighbor) values to match center rank...
    for (auto &r : ranks) {
      if (r == -1) {
        r = cell_ranks[(int)Loc::C];
      }
    }
    Ioss::Utils::uniquify(ranks);
    return ranks;
  }
} // namespace

void Cell::initialize(size_t i, size_t j, std::shared_ptr<UnitCell> unit_cell)
{
  m_i        = i;
  m_j        = j;
  m_unitCell = unit_cell;

  // These are not necessarily the correct ranks, but at this point can determine
  // Whether this cell is surrounded by other cells, or is on the boundary.
  set_rank(Loc::C, 0);
  set_rank(Loc::L, m_i > 0 ? 0 : -1);
  set_rank(Loc::B, m_j > 0 ? 0 : -1);
  set_rank(Loc::BL, (m_i > 0 && m_j > 0) ? 0 : -1);
}

std::pair<double, double> Cell::get_coordinate_range(enum Axis axis) const
{
  if (axis == Axis::X) {
    return m_unitCell->minmax_x;
  }
  if (axis == Axis::Y) {
    return m_unitCell->minmax_y;
  }
  return std::make_pair(0.0, 0.0);
}

// Number of nodes that will be added to global node count when this cell is added to
// grid -- accounts for coincident nodes if cell has neighbor(s)
size_t Cell::added_node_count(enum Mode mode, bool equivalence_nodes) const
{
  // If no neighbors (to -I, -J), then all nodes would be added...
  auto count = m_unitCell->m_region->get_property("node_count").get_int();

  if (equivalence_nodes) {
    if (mode == Mode::GLOBAL) {
      if (has_neighbor_i()) {
        count -= (m_unitCell->cell_JJ * m_unitCell->cell_KK);
      }

      if (has_neighbor_j()) {
        count -= (m_unitCell->cell_II * m_unitCell->cell_KK);
      }

      if (has_neighbor_i() && has_neighbor_j()) {
        count += m_unitCell->cell_KK;
      }
    }
    else if (mode == Mode::PROCESSOR) {
      if (has_neighbor_i() && !processor_boundary(Loc::L)) {
        count -= (m_unitCell->cell_JJ * m_unitCell->cell_KK);
      }

      if (has_neighbor_j() && !processor_boundary(Loc::B)) {
        count -= (m_unitCell->cell_II * m_unitCell->cell_KK);
      }

      if (has_neighbor_i() && has_neighbor_j() && !processor_boundary(Loc::L) &&
          !processor_boundary(Loc::B)) {
        count += m_unitCell->cell_KK;
      }

      // Now the "corner case" ;-) If there is a processor boundary below, but the cell to the BL is
      // on the same rank as this cell, then we have already counted the IJ-line nodes, so need to
      // subtract that count...
      if (processor_boundary(Loc::B) && processor_boundary(Loc::L) &&
          rank(Loc::BL) == rank(Loc::C)) {
        count -= m_unitCell->cell_KK;
      }

      // Now the other "corner case"
      if (processor_boundary(Loc::B) && rank(Loc::BR) == rank(Loc::C)) {
        count -= m_unitCell->cell_KK;
      }
    }
  }
  return count;
}

std::array<int, 9> Cell::categorize_processor_boundary_nodes(int the_rank) const
{
  // Create a "unit cell" to categorize processor boundary nodes...
  std::array<int, 9> bnd_nodes{0};

  // Bottom...
  if (rank(Loc::B) == the_rank) {
    bnd_nodes[(int)Loc::B] = 1;
    if ((rank(Loc::BL) != rank(Loc::C)) && (rank(Loc::L) != rank(Loc::C))) {
      bnd_nodes[(int)Loc::BL] = 1;
    }
    if (rank(Loc::BR) != rank(Loc::C)) {
      bnd_nodes[(int)Loc::BR] = 1;
    }
  }

  // Left
  if (rank(Loc::L) == the_rank) {
    bnd_nodes[(int)Loc::L]  = 1;
    bnd_nodes[(int)Loc::TL] = 1;
    if ((rank(Loc::BL) != rank(Loc::C)) && (rank(Loc::B) != rank(Loc::C))) {
      bnd_nodes[(int)Loc::BL] = 1;
    }
  }

  // Top
  if (rank(Loc::T) == the_rank) {
    bnd_nodes[(int)Loc::T]  = 1;
    bnd_nodes[(int)Loc::TR] = 1;
    if (rank(Loc::L) != rank(Loc::C)) {
      bnd_nodes[(int)Loc::TL] = 1;
    }
  }

  // Right
  if (rank(Loc::R) == the_rank) {
    bnd_nodes[(int)Loc::R]  = 1;
    bnd_nodes[(int)Loc::TR] = 1;
    if ((rank(Loc::BR) != rank(Loc::C)) && (rank(Loc::B) != rank(Loc::C))) {
      bnd_nodes[(int)Loc::BR] = 1;
    }
  }

  // Bottom Left
  if (rank(Loc::BL) == the_rank) {
    // If left and bottom *don't* match the_rank, then need to add this node
    if ((rank(Loc::L) != the_rank) && (rank(Loc::B) != the_rank) &&
        (rank(Loc::L) != rank(Loc::C)) && (rank(Loc::B) != rank(Loc::C))) {
      bnd_nodes[(int)Loc::BL] = 1;
    }
  }

  // Bottom Right
  if (rank(Loc::BR) == the_rank) {
    // If bottom *doesn't* match the_rank, then need to add this node
    if ((rank(Loc::B) != the_rank) && (rank(Loc::B) != rank(Loc::C))) {
      bnd_nodes[(int)Loc::BR] = 1;
    }
  }

  // Top Left
  if (rank(Loc::TL) == the_rank) {
    // If left *doesn't* match the_rank, then need to add this node
    if ((rank(Loc::L) != the_rank) && (rank(Loc::L) != rank(Loc::C))) {
      bnd_nodes[(int)Loc::TL] = 1;
    }
  }

  // Top Right
  if (rank(Loc::TR) == the_rank) {
    bnd_nodes[(int)Loc::TR] = 1;
  }

  return bnd_nodes;
}

size_t Cell::processor_boundary_node_count() const
{
  // Get list of ranks that this cell shares nodes with...
  auto ranks = get_shared_ranks(m_ranks);
  if (ranks.size() == 1) {
    // `ranks` contains center node, so if size == 1, does not touch
    // any other processor.
    return 0;
  }

  // Iterate `ranks` and for each rank, "color" the `bnd_nodes` that that rank touches...
  // Skip center.
  size_t b_count = 0;
  for (auto the_rank : ranks) {
    if (the_rank == rank(Loc::C)) {
      continue;
    }

    // a "unit cell" categorizing processor boundary nodes...
    // Size is 9.  Value is '1' if nodes at this location are shared with rank `rank`
    auto bnd_nodes = categorize_processor_boundary_nodes(the_rank);

    // Now count how many nodes we have added...
    // Edges (B, T, L, R) without corners
    b_count += (m_unitCell->cell_II - 2) * (bnd_nodes[(int)Loc::B] + bnd_nodes[(int)Loc::T]);
    b_count += (m_unitCell->cell_JJ - 2) * (bnd_nodes[(int)Loc::L] + bnd_nodes[(int)Loc::R]);

    // Now the corners (BL, BR, TL, TR)
    b_count += bnd_nodes[(int)Loc::BL] + bnd_nodes[(int)Loc::BR] + bnd_nodes[(int)Loc::TL] +
               bnd_nodes[(int)Loc::TR];
  }

  // The counts above only account for a single KK plane.  Now multiply by `m_unitCell->KK` to get
  // total count.
  b_count *= m_unitCell->cell_KK;
  m_communicationNodeCount = b_count;
  return b_count;
}

template void Cell::populate_node_communication_map(const std::vector<int> &node_map,
                                                    std::vector<int>       &nodes,
                                                    std::vector<int>       &procs) const;
template void Cell::populate_node_communication_map(const std::vector<int64_t> &node_map,
                                                    std::vector<int64_t>       &nodes,
                                                    std::vector<int64_t>       &procs) const;

template <typename INT>
void Cell::populate_node_communication_map(const std::vector<INT> &node_map,
                                           std::vector<INT> &nodes, std::vector<INT> &procs) const
{
  if (m_communicationNodeCount == 0) {
    return;
  }

  nodes.reserve(m_communicationNodeCount);
  procs.reserve(m_communicationNodeCount);

  // Get list of ranks that this cell shares nodes with...
  auto ranks = get_shared_ranks(m_ranks);
  SMART_ASSERT(ranks.size() > 1);

  auto KK = m_unitCell->cell_KK;

  for (auto shared_rank : ranks) {
    if (shared_rank == rank(Loc::C)) {
      continue;
    }

    // a "unit cell" categorizing processor boundary nodes...
    // Size is 9.  Value is '1' if nodes at this location are shared with rank `rank`
    auto bnd_nodes = categorize_processor_boundary_nodes(shared_rank);

    // Handle Edges, but skip nodes on corners.  They are handled later.
    if (bnd_nodes[(int)Loc::B] == 1) {
      process_face_nodes(node_map, nodes, procs, m_unitCell->min_J_face, KK, shared_rank);
    }

    if (bnd_nodes[(int)Loc::T] == 1) {
      process_face_nodes(node_map, nodes, procs, m_unitCell->max_J_face, KK, shared_rank);
    }

    if (bnd_nodes[(int)Loc::L] == 1) {
      process_face_nodes(node_map, nodes, procs, m_unitCell->min_I_face, KK, shared_rank);
    }

    if (bnd_nodes[(int)Loc::R] == 1) {
      process_face_nodes(node_map, nodes, procs, m_unitCell->max_I_face, KK, shared_rank);
    }

    // Now the corners...
    if (bnd_nodes[(int)Loc::BL] == 1) {
      process_corner_nodes(node_map, nodes, procs, m_unitCell->min_J_face, KK, shared_rank,
                           Loc::BL);
    }

    if (bnd_nodes[(int)Loc::BR] == 1) {
      process_corner_nodes(node_map, nodes, procs, m_unitCell->min_J_face, KK, shared_rank,
                           Loc::BR);
    }

    if (bnd_nodes[(int)Loc::TL] == 1) {
      process_corner_nodes(node_map, nodes, procs, m_unitCell->max_J_face, KK, shared_rank,
                           Loc::TL);
    }

    if (bnd_nodes[(int)Loc::TR] == 1) {
      process_corner_nodes(node_map, nodes, procs, m_unitCell->max_J_face, KK, shared_rank,
                           Loc::TR);
    }
  }
  SMART_ASSERT(nodes.size() == procs.size())(nodes.size())(procs.size());
  SMART_ASSERT(nodes.size() == m_communicationNodeCount)(nodes.size())(m_communicationNodeCount);
}

std::vector<int> Cell::categorize_nodes(enum Mode mode) const
{
  auto nodes = m_unitCell->categorize_nodes(has_neighbor_i(), has_neighbor_j());
  if (mode == Mode::PROCESSOR) {
    // If there is a processor boundary to the left, then need to change categorization of
    // all nodes on the left to '0'
    if (processor_boundary(Loc::L)) {
      const auto &min_I_face = m_unitCell->min_I_face;
      for (const auto &node : min_I_face) {
        nodes[node] -= 1;
      }
    }
    if (processor_boundary(Loc::B)) {
      const auto &min_J_face = m_unitCell->min_J_face;
      for (const auto &node : min_J_face) {
        nodes[node] -= 2;
      }
    }

    // Now the "corner case" ;-) If there is a processor boundary below, but the cell to the BL is
    // on the same rank as this cell, then we have already counted the IJ-line nodes, so need to
    // categorize those nodes as already accounted for in a previous map...
    if (processor_boundary(Loc::B) && processor_boundary(Loc::L) && rank(Loc::BL) == rank(Loc::C)) {
      // Want KK() nodes -- First KK of min_i and of min_j.  But since they match, can "unzero"
      // min_i[0..KK)
      for (size_t i = 0; i < m_unitCell->cell_KK; i++) {
        nodes[m_unitCell->min_I_face[i]] = -1;
      }
    }
    // Now the other "corner case"
    if (processor_boundary(Loc::B) && rank(Loc::BR) == rank(Loc::C)) {
      // Want KK() nodes -- First KK of max_i and Last KK of min_j.  But since they match, can
      // "unzero" max_i[0..KK)
      for (size_t i = 0; i < m_unitCell->cell_KK; i++) {
        nodes[m_unitCell->max_I_face[i]] = -1;
      }
    }
  }
  return nodes;
}

template std::vector<int64_t> Cell::generate_node_map(Mode, bool, int64_t) const;
template std::vector<int>     Cell::generate_node_map(Mode, bool, int) const;

template <typename INT>
std::vector<INT> Cell::generate_node_map(Mode mode, bool equivalence_nodes, INT /*dummy*/) const
{
  // Size is node_count + 1 to handle the 1-based connectivity values.
  size_t           cell_node_count = m_unitCell->m_region->get_property("node_count").get_int();
  std::vector<INT> map(cell_node_count + 1);

  INT offset = mode == Mode::PROCESSOR ? m_localNodeIdOffset : m_globalNodeIdOffset;

  if (!equivalence_nodes || !(has_neighbor_i() || has_neighbor_j())) {
    std::iota(map.begin(), map.end(), offset);
  }
  else if (has_neighbor_i() || has_neighbor_j()) {
    // At least one neighboring cell and the nodes are being equivalenced
    // Generate map for the "non-neighbored" nodes (not contiguous with a neighbor cell)
    auto categorized_nodes = categorize_nodes(mode);
    SMART_ASSERT(categorized_nodes.size() == cell_node_count)
    (categorized_nodes.size())(cell_node_count);
    offset++; // To deal with 1-based node numbers.
    for (size_t n = 0; n < cell_node_count; n++) {
      if (categorized_nodes[n] == 0) {
        map[n + 1] = offset++;
      }
    }
  }

  if (equivalence_nodes && has_neighbor_i() &&
      (mode == Mode::GLOBAL || (mode == Mode::PROCESSOR && rank(Loc::C) == rank(Loc::L)))) {
    // Get the neighbor cell...
    // iterate my unit cell's min_I_face() nodes to get index into map
    // At this index, set value to this cells min_I_nodes() node
    // which was created by the neighbor when he was processed...
    SMART_ASSERT(min_I_nodes.size() == m_unitCell->min_I_face.size())
    (m_i)(m_j)(min_I_nodes.size())(m_unitCell->min_I_face.size());

    for (size_t i = 0; i < m_unitCell->min_I_face.size(); i++) {
      auto idx = m_unitCell->min_I_face[i] + 1;
      auto val = min_I_nodes[i];
      map[idx] = (INT)val;
    }
  }

  if (equivalence_nodes && has_neighbor_j() &&
      (mode == Mode::GLOBAL || (mode == Mode::PROCESSOR && rank(Loc::C) == rank(Loc::B)))) {
    SMART_ASSERT(min_J_nodes.size() == m_unitCell->min_J_face.size())
    (m_i)(m_j)(min_J_nodes.size())(m_unitCell->min_J_face.size());

    for (size_t i = 0; i < m_unitCell->min_J_face.size(); i++) {
      auto idx = m_unitCell->min_J_face[i] + 1;
      auto val = min_J_nodes[i];
      map[idx] = (INT)val;
    }
  }

  if (mode == Mode::PROCESSOR) {
    // Now the "corner case" ;-) If there is a processor boundary below, but the cell to the BL is
    // on the same rank as this cell, then we have already counted the IJ-line nodes, so need to
    // categorize those nodes as already accounted for in a previous map...
    if (processor_boundary(Loc::B) && processor_boundary(Loc::L) && rank(Loc::BL) == rank(Loc::C)) {
      // Want KK() nodes -- First KK of min_i and of min_j.  But since they match, can "unzero"
      // min_i[0..KK)
      auto KK = m_unitCell->cell_KK;
      for (size_t i = 0; i < KK; i++) {
        auto idx = m_unitCell->min_J_face[i] + 1;
        auto val = min_J_nodes[i];
        map[idx] = (INT)val;
      }
    }
    // Now the other "corner case"
    if (processor_boundary(Loc::B) && rank(Loc::BR) == rank(Loc::C)) {
      // Want KK() nodes -- First KK of max_i and Last KK of min_j.  But since they match, can
      // "unzero" max_i[0..KK)
      auto KK       = m_unitCell->cell_KK;
      auto j_offset = min_J_nodes.size() - KK;
      for (size_t i = 0; i < KK; i++) {
        auto idx = m_unitCell->min_J_face[j_offset + i] + 1;
        auto val = min_J_nodes[j_offset + i];
        map[idx] = (INT)val;
      }
    }
  }
  // Can now clean out the `min_I_nodes` and `min_J_nodes` lists since the data will no longer be
  // needed.
  Ioss::Utils::clear(min_I_nodes);
  Ioss::Utils::clear(min_J_nodes);

  return map;
}

template void Cell::populate_neighbor(Loc location, const std::vector<int64_t> &map,
                                      const Cell &neighbor) const;
template void Cell::populate_neighbor(Loc location, const std::vector<int> &map,
                                      const Cell &neighbor) const;

template <typename INT>
void Cell::populate_neighbor(Loc location, const std::vector<INT> &map, const Cell &neighbor) const
{
  switch (location) {
  case Loc::L:
    neighbor.min_I_nodes.resize(m_unitCell->max_I_face.size());
    for (size_t i = 0; i < m_unitCell->max_I_face.size(); i++) {
      auto idx                = m_unitCell->max_I_face[i] + 1;
      auto val                = map[idx];
      neighbor.min_I_nodes[i] = val;
    }
    if (debug_level & 8) {
      fmt::print("\nCell {} {}\n", neighbor.m_i, neighbor.m_j);
      fmt::print("min_I_nodes: {}\n", fmt::join(neighbor.min_I_nodes, " "));
    }
    break;

  case Loc::B:
    neighbor.min_J_nodes.resize(m_unitCell->max_J_face.size());
    for (size_t i = 0; i < m_unitCell->max_J_face.size(); i++) {
      auto idx                = m_unitCell->max_J_face[i] + 1;
      auto val                = map[idx];
      neighbor.min_J_nodes[i] = val;
    }
    if (debug_level & 8) {
      fmt::print("min_J_nodes: {}\n", fmt::join(neighbor.min_J_nodes, " "));
    }
    break;

  case Loc::BR: {
    neighbor.min_J_nodes.resize(m_unitCell->max_J_face.size());
    auto KK       = m_unitCell->cell_KK;
    auto j_offset = neighbor.min_J_nodes.size() - KK;
    for (size_t i = 0; i < KK; i++) {
      auto idx                           = m_unitCell->max_J_face[i] + 1;
      auto val                           = map[idx];
      neighbor.min_J_nodes[j_offset + i] = val;
    }
  } break;

  case Loc::BL: {
    neighbor.min_J_nodes.resize(m_unitCell->max_J_face.size());
    auto KK       = m_unitCell->cell_KK;
    auto j_offset = neighbor.min_J_nodes.size() - KK;
    for (size_t i = 0; i < KK; i++) {
      auto idx                = m_unitCell->max_J_face[j_offset + i] + 1;
      auto val                = map[idx];
      neighbor.min_J_nodes[i] = val;
    }
  } break;

  default:
    fmt::print(stderr, "\nINTERNAL ERROR: Unhandled direction in populate_neighbor(): {}\n",
               location);
    exit(EXIT_FAILURE);
  }
}
