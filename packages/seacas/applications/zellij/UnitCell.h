// Copyright(C) 2021, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

#include "Ioss_Region.h"

//! \file

enum class Bnd { MIN_I = 0, MAX_I = 1, MIN_J = 2, MAX_J = 3, MIN_K = 4, MAX_K = 5 };
enum class Flg { MIN_I = 1, MAX_I = 4, MIN_J = 2, MAX_J = 8, MIN_K = 16, MAX_K = 32 };

struct GeneratedSideBlock
{
  //! List of faces (10*element_offset + face) in this SideBlock
  //! indexed by element block name...
  //!
  //! Note that the `element_offset` is the 1-based index of the
  //! element within the element block that it is a member of
  GeneratedSideBlock()                           = default;
  GeneratedSideBlock(const GeneratedSideBlock &) = delete;

  std::map<std::string, std::vector<int64_t>> m_faces;

  size_t size() const
  {
    size_t count = 0;
    for (const auto &faces : m_faces) {
      count += faces.second.size();
    }
    return count;
  }
};

class UnitCell
{
public:
  explicit UnitCell(std::shared_ptr<Ioss::Region> region);
  UnitCell(const UnitCell &) = delete;

  //! Create a vector of `node_count` length which has the following values:
  //! * 0: Node that is not shared with any neighbors.
  //! * 1: Node on `min_I` face
  //! * 2: Node on `min_J` face
  //! * 3: Node on `min_I-min_J` line
  //! If `all_faces` is true, then also categorize the max_I (4) and max_J (8) faces and don't
  //! consider neighbors (want all boundaries marked).
  std::vector<int> categorize_nodes(bool neighbor_i, bool neighbor_j, bool all_faces = false) const;

  std::shared_ptr<Ioss::Region> m_region{nullptr};

  //!@{ The local node ids of the nodes that are on each face of this
  //! unit cell.  The I/J face nodes are in a structured configuration
  //! and these vectors are sorted such that the `min_I_face` and
  //! `max_I_face` nodes are at the same "parametric" location.  In
  //! other words, the nodes in the `max_I_face` of one unit cell will
  //! line up with the nodes in the `min_I_face` of the neighboring
  //! cell.  They are also ordered such that the first `cell_KK` nodes
  //! in the `min_I_face` and `min_J_face` lists should be the same
  //! and are the nodes at the intersection of those two faces.
  std::vector<int64_t> min_I_face{};
  std::vector<int64_t> max_I_face{};
  std::vector<int64_t> min_J_face{};
  std::vector<int64_t> max_J_face{};
  //!@}

  ///@{
  //! A pair containing the
  //! minimum and maximum coordinate extent in the `x` and `y`
  //! directions.
  std::pair<double, double> minmax_x{};
  std::pair<double, double> minmax_y{};
  ///@}

  void generate_boundary_faces(unsigned int which_faces);

  //! Used by `generate_boundary_faces()` to categorize nodes on the +/- Z faces of unit cell.
  void categorize_z_nodes(std::vector<int> &categorized_nodes);

  std::array<GeneratedSideBlock, 6> boundary_blocks{};

  ///@{
  //! The outer boundary of a UnitCell has a
  //! structured-configuration of the boundary faces on the non-K
  //! faces.  The `cell_II`, `cell_JJ`, and `cell_KK` variables give
  //! the dimensions of this structured face mesh for the `min_I`,
  //! `max_I`, `min_J`, and `max_J` faces.  Note that the `min_K` and
  //! `max_K` faces are *NOT* structured, but they will have `cell_II`
  //! and `cell_JJ` nodes on the boundary of the faces.
  size_t cell_II{};
  size_t cell_JJ{};
  size_t cell_KK{};
  ///@}
};

using UnitCellMap = std::map<std::string, std::shared_ptr<UnitCell>>;
