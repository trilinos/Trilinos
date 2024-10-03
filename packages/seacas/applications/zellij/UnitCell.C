// Copyright(C) 2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "UnitCell.h"
#include <vector>

#include "Ioss_ElementBlock.h"
#include "Ioss_FaceGenerator.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_Region.h"
#include "Ioss_SmartAssert.h"
#include "Ioss_Sort.h"
#include "fmt/format.h"
#include "fmt/ranges.h"

//! \file

extern unsigned int debug_level;

namespace {
  bool on_boundary(Flg ijk, const Ioss::Face &face, std::vector<int> &categorized_nodes)
  {
    int result = (unsigned int)ijk & categorized_nodes[face.connectivity_[0] - 1] &
                 categorized_nodes[face.connectivity_[1] - 1] &
                 categorized_nodes[face.connectivity_[2] - 1] &
                 categorized_nodes[face.connectivity_[3] - 1];
    return result != 0;
  }

  Ioss::ElementBlock *get_element_block(Ioss::ElementBlock *block, int64_t elem_id,
                                        std::shared_ptr<Ioss::Region> &region)
  {
    auto *new_block = block;
    if (block == nullptr || !block->contains(elem_id)) {
      new_block = region->get_element_block(elem_id);
      assert(new_block != nullptr);
    }
    return new_block;
  }

  bool approx_equal(double A, double B)
  {
    static double maxRelDiff = 1000.0 * std::numeric_limits<double>::epsilon();
    static double maxDiff    = 100.0 * maxRelDiff;

    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = std::abs(A - B);
    if (diff <= maxDiff) {
      return true;
    }

    A              = std::abs(A);
    B              = std::abs(B);
    double largest = std::max(B, A);

    return (diff <= largest * maxRelDiff);
  }

  void gather_face_nodes(std::vector<double> &coord, std::pair<double, double> &minmax,
                         std::vector<int64_t> &min_face, std::vector<int64_t> &max_face)
  {
    for (size_t i = 0; i < coord.size(); i++) {
      if (approx_equal(coord[i], minmax.first)) {
        min_face.push_back(i);
      }
      if (approx_equal(coord[i], minmax.second)) {
        max_face.push_back(i);
      }
    }
    SMART_ASSERT(min_face.size() == max_face.size())(min_face.size())(max_face.size());
  }

  template <typename INT>
  void sort_face_nodes(std::vector<INT> &face_nodes, const std::vector<double> &coord_j,
                       const std::vector<double> &coord_i)
  {
    // NOTE: The `!approx_equal(coord_j[a], coord_j[b])` portion of the less than comparison
    // is meant to handle values close to but one or both are not exactly zero.
    Ioss::sort(face_nodes.begin(), face_nodes.end(), [&coord_j](size_t a, size_t b) {
      return !approx_equal(coord_j[a], coord_j[b]) && float(coord_j[a]) < float(coord_j[b]);
    });
    std::stable_sort(face_nodes.begin(), face_nodes.end(), [&coord_i](size_t a, size_t b) {
      return !approx_equal(coord_i[a], coord_i[b]) && float(coord_i[a]) < float(coord_i[b]);
    });
  }
} // namespace

UnitCell::UnitCell(std::shared_ptr<Ioss::Region> region) : m_region(region)
{
  std::vector<double> coord_x;
  std::vector<double> coord_y;
  std::vector<double> coord_z;

  auto *nb = region->get_node_blocks()[0];
  nb->get_field_data("mesh_model_coordinates_x", coord_x);
  nb->get_field_data("mesh_model_coordinates_y", coord_y);
  nb->get_field_data("mesh_model_coordinates_z", coord_z);

  const auto x_min_max_it = std::minmax_element(coord_x.begin(), coord_x.end());
  const auto y_min_max_it = std::minmax_element(coord_y.begin(), coord_y.end());

  minmax_x = std::make_pair(*x_min_max_it.first, *x_min_max_it.second);
  minmax_y = std::make_pair(*y_min_max_it.first, *y_min_max_it.second);

  if (debug_level & 4) {
    fmt::print("Min / Max X = {} ... {}\n", minmax_x.first, minmax_x.second);
    fmt::print("Min / Max Y = {} ... {}\n", minmax_y.first, minmax_y.second);
  }
  // Now iterate all nodes and categorize if on a face -- minx, maxx, miny, maxy,
  gather_face_nodes(coord_x, minmax_x, min_I_face, max_I_face);
  gather_face_nodes(coord_y, minmax_y, min_J_face, max_J_face);
  if (debug_level & 4) {
    fmt::print("Nodes on X Face = {}\n", min_I_face.size());
    fmt::print("Nodes on Y Face = {}\n", min_J_face.size());
  }

  sort_face_nodes(min_I_face, coord_z, coord_y);
  sort_face_nodes(max_I_face, coord_z, coord_y);
  sort_face_nodes(min_J_face, coord_z, coord_x);
  sort_face_nodes(max_J_face, coord_z, coord_x);
  if (debug_level & 4) {
    // Output each set of nodes --
    fmt::print("\nSORTED:\n");
    fmt::print("\tMin/Max I Face:\n");
    for (size_t i = 0; i < min_I_face.size(); i++) {
      auto min_I = min_I_face[i];
      auto max_I = max_I_face[i];
      fmt::print("\t\t{:10}: {:12.4e} {:12.4e} {:12.4e}\t{:10}: {:12.4e} {:12.4e} {:12.4e}\n",
                 min_I, coord_x[min_I], coord_y[min_I], coord_z[min_I], max_I, coord_x[max_I],
                 coord_y[max_I], coord_z[max_I]);
    }
    fmt::print("\tMin/Max J Face:\n");
    for (size_t i = 0; i < min_J_face.size(); i++) {
      auto min_J = min_J_face[i];
      auto max_J = max_J_face[i];
      fmt::print("\t\t{:10}: {:12.4e} {:12.4e} {:12.4e}\t{:10}: {:12.4e} {:12.4e} {:12.4e}\n",
                 min_J, coord_x[min_J], coord_y[min_J], coord_z[min_J], max_J, coord_x[max_J],
                 coord_y[max_J], coord_z[max_J]);
    }
  }

#ifndef NDEBUG
  for (size_t i = 0; i < min_I_face.size(); i++) {
    auto minI = min_I_face[i];
    auto maxI = max_I_face[i];
    SMART_ASSERT((size_t)minI < coord_y.size());
    SMART_ASSERT((size_t)maxI < coord_y.size());
    SMART_ASSERT((size_t)minI < coord_z.size());
    SMART_ASSERT((size_t)maxI < coord_z.size());
    SMART_ASSERT(approx_equal(coord_y[minI], coord_y[maxI]))(coord_y[minI])(coord_y[maxI]);
    SMART_ASSERT(approx_equal(coord_z[minI], coord_z[maxI]))(coord_z[minI])(coord_z[maxI]);
  }
  for (size_t i = 0; i < min_J_face.size(); i++) {
    auto minJ = min_J_face[i];
    auto maxJ = max_J_face[i];
    SMART_ASSERT((size_t)minJ < coord_x.size());
    SMART_ASSERT((size_t)maxJ < coord_x.size());
    SMART_ASSERT((size_t)minJ < coord_z.size());
    SMART_ASSERT((size_t)maxJ < coord_z.size());
    SMART_ASSERT(approx_equal(coord_x[minJ], coord_x[maxJ]))(coord_x[minJ])(coord_z[maxJ]);
    SMART_ASSERT(approx_equal(coord_z[minJ], coord_z[maxJ]))(coord_z[minJ])(coord_z[maxJ]);
  }
#endif

  // Determine 'K' -- size of minI_minJ corner list.
  cell_KK = 0;
  for (; cell_KK < min_I_face.size(); cell_KK++) {
    if (min_I_face[cell_KK] != min_J_face[cell_KK]) {
      break;
    }

#ifndef NDEBUG
    // Checking that we get the same on the maxI_maxJ corner...
    size_t i_x = max_I_face.size() - 1 - cell_KK;
    size_t i_y = max_J_face.size() - 1 - cell_KK;
    SMART_ASSERT(max_I_face[i_x] == max_J_face[i_y])
    (cell_KK)(i_x)(i_y)(max_I_face[i_x])(max_J_face[i_y]);
#endif
  }

  SMART_ASSERT(min_I_face.size() % cell_KK == 0)(min_I_face.size())(cell_KK);
  SMART_ASSERT(min_J_face.size() % cell_KK == 0)(min_J_face.size())(cell_KK);

  cell_II = min_J_face.size() / cell_KK;
  cell_JJ = min_I_face.size() / cell_KK;

  if (debug_level & 4) {
    fmt::print("\nUnitCell {}:\n", m_region->name());
    fmt::print("\tThe minI face contains {} nodes\n", min_I_face.size());
    fmt::print("\tThe minJ face contains {} nodes\n", min_J_face.size());
    fmt::print("\tThe calculated cell shape is {} x {} x {}\n", cell_II, cell_JJ, cell_KK);
  }
}

std::vector<int> UnitCell::categorize_nodes(bool neighbor_i, bool neighbor_j, bool all_faces) const
{
  // Create a vector of `node_count` length which has the following values:
  // 0: Node that is not shared with any neighbors.
  // 1: Node on min_I face
  // 2: Node on min_J face
  // 3: Node on min_I-min_J line

  auto             node_count = m_region->get_property("node_count").get_int();
  std::vector<int> node_category(node_count);

  if (neighbor_i || all_faces) {
    for (auto node : min_I_face) {
      node_category[node] = (int)Flg::MIN_I;
    }
  }

  if (neighbor_j || all_faces) {
    for (auto node : min_J_face) {
      node_category[node] += (int)Flg::MIN_J;
    }
  }

  if (all_faces) {
    for (auto node : max_I_face) {
      node_category[node] += (int)Flg::MAX_I;
    }

    for (auto node : max_J_face) {
      node_category[node] += (int)Flg::MAX_J;
    }
  }
  return node_category;
}

void UnitCell::categorize_z_nodes(std::vector<int> &categorized_nodes)
{
  std::vector<double> coord_z;

  auto *nb = m_region->get_node_blocks()[0];
  nb->get_field_data("mesh_model_coordinates_z", coord_z);

  const auto z_min_max_it = std::minmax_element(coord_z.begin(), coord_z.end());
  auto       minmax_z     = std::make_pair(*z_min_max_it.first, *z_min_max_it.second);

  std::vector<int64_t> min_K_face;
  std::vector<int64_t> max_K_face;
  gather_face_nodes(coord_z, minmax_z, min_K_face, max_K_face);

  for (auto node : min_K_face) {
    categorized_nodes[node] += (int)Flg::MIN_K;
  }
  for (auto node : max_K_face) {
    categorized_nodes[node] += (int)Flg::MAX_K;
  }
}

void UnitCell::generate_boundary_faces(unsigned int which_faces)
{
  Ioss::FaceGenerator face_generator(*m_region);
  bool                block_by_block = false;
  bool                local_ids      = true;
  if (m_region->get_database()->int_byte_size_api() == 4) {
    face_generator.generate_faces((int)0, block_by_block, local_ids);
  }
  else {
    face_generator.generate_faces((int64_t)0, block_by_block, local_ids);
  }

  // Get vector of all boundary faces which need to be categorized into
  // which face of the unit cell they are on
  auto categorized_nodes = categorize_nodes(false, false, true);
  categorize_z_nodes(categorized_nodes);

  if (debug_level & 128) {
    fmt::print("Node Category: {}\n", fmt::join(categorized_nodes, " "));
  }

  std::array<enum Flg, 6> boundary_flag{Flg::MIN_I, Flg::MAX_I, Flg::MIN_J,
                                        Flg::MAX_J, Flg::MIN_K, Flg::MAX_K};
  auto                   &faces = face_generator.faces("ALL");
  Ioss::ElementBlock     *block = nullptr;
  for (const auto &face : faces) {
    if (face.element_count() == 1) {
      block             = get_element_block(block, face.element[0] / 10, m_region);
      auto block_offset = block->get_offset();
      for (int i = 0; i < 6; i++) {
        if ((which_faces & (unsigned)boundary_flag[i]) &&
            on_boundary(boundary_flag[i], face, categorized_nodes)) {

          // Complexity:  The `face.element[0]` is 10 * local_element_id + face.
          // We need to convert that to the offset within the element block
          // For example, if we have two element blocks with 100 elements each, then
          // local_element_id 128 would be the 28 element in the second element block
          // and we would need to store '28 * 10 + face' instead of '128 * 10 + face'
          // so when we process this element later, we can correctly calculate the
          // local_element_id in the output file...
          auto unit_local_element_id = face.element[0] / 10;
          auto face_ordinal          = face.element[0] % 10;
          auto unit_block_location   = (unit_local_element_id - block_offset) * 10 + face_ordinal;

          if (debug_level & 128) {
            fmt::print("Element {}, Side {} is on boundary {}, block {}\n",
                       unit_block_location / 10, unit_block_location % 10, i, block->name());
          }
          boundary_blocks[i].m_faces[block->name()].push_back(unit_block_location);
          break;
        }
      }
    }
  }

  if (which_faces & (unsigned)Flg::MIN_I) {
    SMART_ASSERT(boundary_blocks[(int)Bnd::MIN_I].size() == (cell_JJ - 1) * (cell_KK - 1))
    (boundary_blocks[(int)Bnd::MIN_I].size())(cell_JJ - 1)(cell_KK - 1);
  }
  if (which_faces & (unsigned)Flg::MAX_I) {
    SMART_ASSERT(boundary_blocks[(int)Bnd::MAX_I].size() == (cell_JJ - 1) * (cell_KK - 1))
    (boundary_blocks[(int)Bnd::MAX_I].size())(cell_JJ - 1)(cell_KK - 1);
  }

  if (which_faces & (unsigned)Flg::MIN_J) {
    SMART_ASSERT(boundary_blocks[(int)Bnd::MIN_J].size() == (cell_II - 1) * (cell_KK - 1))
    (boundary_blocks[(int)Bnd::MIN_J].size())(cell_II - 1)(cell_KK - 1);
  }
  if (which_faces & (unsigned)Flg::MAX_J) {
    SMART_ASSERT(boundary_blocks[(int)Bnd::MAX_J].size() == (cell_II - 1) * (cell_KK - 1))
    (boundary_blocks[(int)Bnd::MAX_J].size())(cell_II - 1)(cell_KK - 1);
  }
}
