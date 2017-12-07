// Copyright(C) 1999-2010 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Ioss_BoundingBox.h"  // for AxisAlignedBoundingBox
#include "Ioss_FieldManager.h" // for FieldManager
#include <Ioss_DatabaseIO.h>   // for DatabaseIO
#include <Ioss_Field.h>        // for Field, etc
#include <Ioss_Property.h>     // for Property
#include <Ioss_Region.h>
#include <Ioss_StructuredBlock.h>
#include <cstddef> // for size_t
#include <numeric>
#include <string> // for string
#include <vector> // for vector

namespace {
  const std::string SCALAR() { return std::string("scalar"); }
  const std::string VECTOR_2D() { return std::string("vector_2d"); }
  const std::string VECTOR_3D() { return std::string("vector_3d"); }

  int sign(int value) { return value < 0 ? -1 : 1; }

  int del(int v1, int v2) { return static_cast<int>(std::abs(v1) == std::abs(v2)); }

} // namespace

namespace Ioss {
  class Field;

  /** \brief Create a structured block.
   *
   *  \param[in] io_database The database associated with the region containing the structured
   * block.
   *  \param[in] my_name The structured block's name.
   *  \param[in] index_dim The dimensionality of the block -- 1D, 2D, 3D
   *  \param[in] ni The number of intervals in the (i) direction.
   *  \param[in] nj The number of intervals in the (j) direction. Zero if 1D
   *  \param[in] nk The number of intervals in the (k) direction. Zero if 2D
   */
  StructuredBlock::StructuredBlock(DatabaseIO *io_database, const std::string &my_name,
                                   int index_dim, int ni, int nj, int nk, int off_i, int off_j,
                                   int off_k)
      : EntityBlock(io_database, my_name, "Hex8", ni * (nj > 0 ? nj : 1) * (nk > 0 ? nk : 1)),
        m_ni(ni), m_nj(nj), m_nk(nk), m_offsetI(off_i), m_offsetJ(off_j), m_offsetK(off_k),
        m_niGlobal(m_ni), m_njGlobal(m_nj), m_nkGlobal(m_nk), m_nodeOffset(0), m_cellOffset(0),
        m_nodeGlobalOffset(0), m_cellGlobalOffset(0),
        m_nodeBlock(io_database, my_name + "_nodes", (m_ni + 1) * (m_nj + 1) * (m_nk + 1),
                    index_dim)
  {
    add_properties_and_fields(index_dim);
  }

  StructuredBlock::StructuredBlock(DatabaseIO *io_database, const std::string &my_name,
                                   int index_dim, Ioss::IJK_t &ordinal, Ioss::IJK_t &offset,
                                   Ioss::IJK_t &global_ordinal)
      : EntityBlock(io_database, my_name, "Hex8", ordinal[0] * ordinal[1] * ordinal[2]),
        m_ni(ordinal[0]), m_nj(ordinal[1]), m_nk(ordinal[2]), m_offsetI(offset[0]),
        m_offsetJ(offset[1]), m_offsetK(offset[2]), m_niGlobal(global_ordinal[0]),
        m_njGlobal(global_ordinal[1]), m_nkGlobal(global_ordinal[2]), m_nodeOffset(0),
        m_cellOffset(0), m_nodeGlobalOffset(0), m_cellGlobalOffset(0),
        m_nodeBlock(io_database, my_name + "_nodes", (m_ni + 1) * (m_nj + 1) * (m_nk + 1),
                    index_dim)
  {
    add_properties_and_fields(index_dim);
  }

  void StructuredBlock::add_properties_and_fields(int index_dim)
  {
    assert(index_dim == 1 || index_dim == 2 || index_dim == 3);

    int64_t cell_count        = 0;
    int64_t node_count        = 0;
    int64_t global_cell_count = 0;
    int64_t global_node_count = 0;

    if (index_dim == 1) {
      cell_count = m_ni;
      node_count = cell_count == 0 ? 0 : (m_ni + 1);

      global_cell_count = m_niGlobal;
      global_node_count = global_cell_count == 0 ? 0 : (m_niGlobal + 1);
    }
    else if (index_dim == 2) {
      cell_count = static_cast<int64_t>(m_ni) * m_nj;
      node_count = cell_count == 0 ? 0 : static_cast<int64_t>(m_ni + 1) * (m_nj + 1);

      global_cell_count = static_cast<int64_t>(m_niGlobal) * m_njGlobal;
      global_node_count =
          global_cell_count == 0 ? 0 : static_cast<int64_t>(m_niGlobal + 1) * (m_njGlobal + 1);
    }
    else if (index_dim == 3) {
      cell_count = static_cast<int64_t>(m_ni) * m_nj * m_nk;
      node_count = cell_count == 0 ? 0 : static_cast<int64_t>(m_ni + 1) * (m_nj + 1) * (m_nk + 1);

      global_cell_count = static_cast<int64_t>(m_niGlobal) * m_njGlobal * m_nkGlobal;
      global_node_count = global_cell_count == 0 ? 0
                                                 : static_cast<int64_t>(m_niGlobal + 1) *
                                                       (m_njGlobal + 1) * (m_nkGlobal + 1);
    }

    properties.add(Property("component_degree", index_dim));
    properties.add(Property("node_count", node_count));
    properties.add(Property("cell_count", cell_count));
    properties.add(Property("global_node_count", global_node_count));
    properties.add(Property("global_cell_count", global_cell_count));

    properties.add(Property("ni", m_ni));
    properties.add(Property("nj", m_nj));
    properties.add(Property("nk", m_nk));

    properties.add(Property("ni_global", m_niGlobal));
    properties.add(Property("nj_global", m_njGlobal));
    properties.add(Property("nk_global", m_nkGlobal));

    properties.add(Property("offset_i", m_offsetI));
    properties.add(Property("offset_j", m_offsetJ));
    properties.add(Property("offset_k", m_offsetK));

    std::string vector_name;
    if (index_dim == 1) {
      vector_name = SCALAR();
    }
    else if (index_dim == 2) {
      vector_name = VECTOR_2D();
    }
    else if (index_dim == 3) {
      vector_name = VECTOR_3D();
    }
    fields.add(
        Ioss::Field("cell_ids", Ioss::Field::INTEGER, SCALAR(), Ioss::Field::MESH, cell_count));

    fields.add(Ioss::Field("cell_node_ids", Ioss::Field::INTEGER, SCALAR(), Ioss::Field::MESH,
                           node_count));

    fields.add(Ioss::Field("mesh_model_coordinates", Ioss::Field::REAL, vector_name,
                           Ioss::Field::MESH, node_count));

    // Permit access 1-coordinate at a time
    fields.add(Ioss::Field("mesh_model_coordinates_x", Ioss::Field::REAL, SCALAR(),
                           Ioss::Field::MESH, node_count));
    if (index_dim > 1) {
      fields.add(Ioss::Field("mesh_model_coordinates_y", Ioss::Field::REAL, SCALAR(),
                             Ioss::Field::MESH, node_count));
    }

    if (index_dim > 2) {
      fields.add(Ioss::Field("mesh_model_coordinates_z", Ioss::Field::REAL, SCALAR(),
                             Ioss::Field::MESH, node_count));
    }
  }

  StructuredBlock::~StructuredBlock() = default;

  StructuredBlock *StructuredBlock::clone(DatabaseIO *database) const
  {
    int index_dim = properties.get("component_degree").get_int();

    IJK_t ijk{{m_ni, m_nj, m_nk}};
    IJK_t offset{{m_offsetI, m_offsetJ, m_offsetK}};
    IJK_t ijk_glob{{m_niGlobal, m_njGlobal, m_nkGlobal}};

    auto block = new StructuredBlock(database, name(), index_dim, ijk, offset, ijk_glob);

    block->m_zoneConnectivity   = m_zoneConnectivity;
    block->m_boundaryConditions = m_boundaryConditions;

    return block;
  }

  Property StructuredBlock::get_implicit_property(const std::string &my_name) const
  {
    return GroupingEntity::get_implicit_property(my_name);
  }

  int64_t StructuredBlock::internal_get_field_data(const Field &field, void *data,
                                                   size_t data_size) const
  {
    return get_database()->get_field(this, field, data, data_size);
  }

  int64_t StructuredBlock::internal_put_field_data(const Field &field, void *data,
                                                   size_t data_size) const
  {
    return get_database()->put_field(this, field, data, data_size);
  }

  AxisAlignedBoundingBox StructuredBlock::get_bounding_box() const
  {
    return get_database()->get_bounding_box(this);
  }

  std::ostream &operator<<(std::ostream &os, const ZoneConnectivity &zgc)
  {
    std::array<std::string, 7> tf{{"-k", "-j", "-i", " ", "i", "j", "k"}};

    // 0 -3 -k
    // 1 -2 -j
    // 2 -1 -i
    // 3
    // 4  1  i
    // 5  2  j
    // 6  3  k
    std::string transform = "[i..";
    transform += tf[zgc.m_transform[0] + 3];
    transform += " j..";
    transform += tf[zgc.m_transform[1] + 3];
    transform += " k..";
    transform += tf[zgc.m_transform[2] + 3];
    transform += "] ";

    os << "\t\t" << zgc.m_donorName << "[P" << zgc.m_donorProcessor << "]:\tDZ " << zgc.m_donorZone
       << "\tName '" << zgc.m_connectionName << "' shares " << zgc.get_shared_node_count()
       << " nodes. (Owned = " << (zgc.owns_shared_nodes() ? "true" : "false") << ")."
       << "\n\t\t\t\tRange: [" << zgc.m_rangeBeg[0] << ".." << zgc.m_rangeEnd[0] << ", "
       << zgc.m_rangeBeg[1] << ".." << zgc.m_rangeEnd[1] << ", " << zgc.m_rangeBeg[2] << ".."
       << zgc.m_rangeEnd[2] << "]\tDonor Range: [" << zgc.m_donorRangeBeg[0] << ".."
       << zgc.m_donorRangeEnd[0] << ", " << zgc.m_donorRangeBeg[1] << ".." << zgc.m_donorRangeEnd[1]
       << ", " << zgc.m_donorRangeBeg[2] << ".." << zgc.m_donorRangeEnd[2] << "]";
    return os;
  }

  std::vector<int> ZoneConnectivity::get_range(int ordinal) const
  {
    // Return the integer values for the specified range for the specified ordinal (1,2,3) ->
    // (i,j,k)
    ordinal--;
    int size  = std::abs(m_rangeBeg[ordinal] - m_rangeEnd[ordinal]) + 1;
    int delta = sign(m_rangeEnd[ordinal] - m_rangeBeg[ordinal]);
    assert(delta == 1 || delta == -1);

    std::vector<int> range(size);
    for (int i = 0; i < size; i++) {
      range[i] = m_rangeBeg[ordinal] + i * delta;
    }
    return range;
  }

  std::array<int, 9> ZoneConnectivity::transform_matrix() const
  {
    std::array<int, 9> t_matrix{};
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        t_matrix[3 * i + j] = sign(m_transform[j]) * del(m_transform[j], i + 1);
      }
    }
    return t_matrix;
  }

  Ioss::IJK_t ZoneConnectivity::transform(const Ioss::IJK_t &index_1) const
  {
    auto t_matrix = transform_matrix();

    Ioss::IJK_t diff{};
    Ioss::IJK_t donor{};

    diff[0] = index_1[0] - m_rangeBeg[0];
    diff[1] = index_1[1] - m_rangeBeg[1];
    diff[2] = index_1[2] - m_rangeBeg[2];

    donor[0] =
        t_matrix[0] * diff[0] + t_matrix[1] * diff[1] + t_matrix[2] * diff[2] + m_donorRangeBeg[0];
    donor[1] =
        t_matrix[3] * diff[0] + t_matrix[4] * diff[1] + t_matrix[5] * diff[2] + m_donorRangeBeg[1];
    donor[2] =
        t_matrix[6] * diff[0] + t_matrix[7] * diff[1] + t_matrix[8] * diff[2] + m_donorRangeBeg[2];

    assert(std::fabs(donor[0] - m_donorRangeBeg[0]) <=
           std::fabs(m_donorRangeBeg[0] - m_donorRangeEnd[0]));
    assert(std::fabs(donor[1] - m_donorRangeBeg[1]) <=
           std::fabs(m_donorRangeBeg[1] - m_donorRangeEnd[1]));
    assert(std::fabs(donor[2] - m_donorRangeBeg[2]) <=
           std::fabs(m_donorRangeBeg[2] - m_donorRangeEnd[2]));
    return donor;
  }

  // ----------------------------------------------------------------------------

  Ioss::IJK_t ZoneConnectivity::inverse_transform(const Ioss::IJK_t &index_1) const
  {
    auto t_matrix = transform_matrix();

    Ioss::IJK_t diff{};
    Ioss::IJK_t index{};

    diff[0] = index_1[0] - m_donorRangeBeg[0];
    diff[1] = index_1[1] - m_donorRangeBeg[1];
    diff[2] = index_1[2] - m_donorRangeBeg[2];

    index[0] =
        t_matrix[0] * diff[0] + t_matrix[3] * diff[1] + t_matrix[6] * diff[2] + m_rangeBeg[0];
    index[1] =
        t_matrix[1] * diff[0] + t_matrix[4] * diff[1] + t_matrix[7] * diff[2] + m_rangeBeg[1];
    index[2] =
        t_matrix[2] * diff[0] + t_matrix[5] * diff[1] + t_matrix[8] * diff[2] + m_rangeBeg[2];

    return index;
  }

  int BoundaryCondition::which_parent_face() const
  {
    // Determine which "face" of the parent block this BC is applied to.
    // min X, max X, min Y, max Y, min Z, max Z -- -1, 1, -2, 2, -3, 3
    if (m_rangeBeg[0] == m_rangeEnd[0]) {
      return (m_rangeBeg[0] == 1) ? -1 : 1;
    }
    if (m_rangeBeg[1] == m_rangeEnd[1]) {
      return (m_rangeBeg[1] == 1) ? -2 : 2;
    }
    return (m_rangeBeg[2] == 1) ? -3 : 3;
  }

  std::ostream &operator<<(std::ostream &os, const BoundaryCondition &bc)
  {
    os << "\t\tBC Name '" << bc.m_bcName << "' owns " << bc.get_face_count() << " faces."
       << "\n\t\t\t\tRange: [" << bc.m_rangeBeg[0] << ".." << bc.m_rangeEnd[0] << ", "
       << bc.m_rangeBeg[1] << ".." << bc.m_rangeEnd[1] << ", " << bc.m_rangeBeg[2] << ".."
       << bc.m_rangeEnd[2] << "]";
    return os;
  }

} // namespace Ioss
