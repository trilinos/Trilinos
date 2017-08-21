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

#ifndef IOSS_Ioss_StructuredBlock_h
#define IOSS_Ioss_StructuredBlock_h

#include <Ioss_BoundingBox.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_EntityBlock.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_Property.h>
#include <array>
#include <cassert>
#include <string>
namespace Ioss {
  class Region;

  struct ZoneConnectivity
  {
    ZoneConnectivity(const std::string name, int owner_zone, const std::string donor_name,
                     int donor_zone, const Ioss::IJK_t p_transform, const Ioss::IJK_t range_beg,
                     const Ioss::IJK_t range_end, const Ioss::IJK_t donor_beg,
                     const Ioss::IJK_t donor_end, bool owns_nodes, bool intra_block = false)
        : m_connectionName(std::move(name)), m_donorName(std::move(donor_name)),
          m_transform(std::move(p_transform)), m_rangeBeg(std::move(range_beg)),
          m_rangeEnd(std::move(range_end)), m_donorRangeBeg(std::move(donor_beg)),
          m_donorRangeEnd(std::move(donor_end)), m_ownerZone(owner_zone), m_donorZone(donor_zone),
          m_ownerProcessor(-1), m_donorProcessor(-1), m_sameRange(false),
          m_ownsSharedNodes(owns_nodes), m_intraBlock(intra_block), m_isActive(true)
    {
      if (!m_intraBlock) {
        m_ownerRange[0] = m_rangeBeg[0];
        m_ownerRange[1] = m_rangeBeg[1];
        m_ownerRange[2] = m_rangeBeg[2];
        m_ownerRange[3] = m_rangeEnd[0];
        m_ownerRange[4] = m_rangeEnd[1];
        m_ownerRange[5] = m_rangeEnd[2];

        m_donorRange[0] = m_donorRangeBeg[0];
        m_donorRange[1] = m_donorRangeBeg[1];
        m_donorRange[2] = m_donorRangeBeg[2];
        m_donorRange[3] = m_donorRangeEnd[0];
        m_donorRange[4] = m_donorRangeEnd[1];
        m_donorRange[5] = m_donorRangeEnd[2];
      }
    }

    ZoneConnectivity(const ZoneConnectivity &copy_from) = default;

    // Return number of nodes in the connection shared with the donor zone.
    size_t get_shared_node_count() const
    {
      size_t snc = 1;
      for (int i = 0; i < 3; i++) {
        snc *= (std::abs(m_rangeEnd[i] - m_rangeBeg[i]) + 1);
      }
      return snc;
    }

    bool owns_shared_nodes() const { return m_ownsSharedNodes; }

    std::array<int, 9> transform_matrix() const;
    Ioss::IJK_t transform(const Ioss::IJK_t &index_1) const;
    Ioss::IJK_t inverse_transform(const Ioss::IJK_t &index_1) const;

    std::vector<int> get_range(int ordinal) const;

    // The "original" owner and donor range -- that is, they have not been subsetted
    // due to block decompositions in a parallel run.  These should be the same on
    // all processors...  Primarily used to make parallel collective output easier...
    std::array<int, 6> m_ownerRange{};
    std::array<int, 6> m_donorRange{};

    std::string m_connectionName;
    std::string m_donorName;
    Ioss::IJK_t m_transform;
    Ioss::IJK_t m_rangeBeg;
    Ioss::IJK_t m_rangeEnd;
    Ioss::IJK_t m_donorRangeBeg;
    Ioss::IJK_t m_donorRangeEnd;

    friend std::ostream &operator<<(std::ostream &os, const ZoneConnectivity &zgc);

    // NOTE: Shared nodes are "owned" by the zone with the lowest zone id.
    int  m_ownerZone;      // "id" of zone that owns this connection
    int  m_donorZone;      // "id" of zone that is donor of this connection
    int  m_ownerProcessor; // processor that owns the owner zone
    int  m_donorProcessor; // processor that owns the donor zone
    bool m_sameRange; // True if owner and donor range should always match...(special use during
                      // decomp)
    bool m_ownsSharedNodes;
    bool m_intraBlock; // True if this zc is created due to processor decompositions in a parallel
                       // run.
    bool m_isActive;   // True if non-zero range...
  };

  struct BoundaryCondition
  {
    BoundaryCondition(const std::string name, const Ioss::IJK_t range_beg,
                      const Ioss::IJK_t range_end)
        : m_bcName(std::move(name)), m_rangeBeg(std::move(range_beg)),
          m_rangeEnd(std::move(range_end))
    {
      m_ownerRange[0] = m_rangeBeg[0];
      m_ownerRange[1] = m_rangeBeg[1];
      m_ownerRange[2] = m_rangeBeg[2];
      m_ownerRange[3] = m_rangeEnd[0];
      m_ownerRange[4] = m_rangeEnd[1];
      m_ownerRange[5] = m_rangeEnd[2];

#ifndef NDEBUG
      int same_count = (m_rangeBeg[0] == m_rangeEnd[0] ? 1 : 0) +
                       (m_rangeBeg[1] == m_rangeEnd[1] ? 1 : 0) +
                       (m_rangeBeg[2] == m_rangeEnd[2] ? 1 : 0);
      assert(same_count == 1 || (same_count == 3 && m_rangeBeg[0] == 0));
#endif
    }

    BoundaryCondition(const BoundaryCondition &copy_from) = default;

    // Determine which "face" of the parent block this BC is applied to.
    int which_parent_face() const;

    // Return number of cell faces in the BC
    size_t get_face_count() const
    {
      if (m_rangeBeg[0] == 0 || m_rangeEnd[0] == 0 || m_rangeBeg[1] == 0 || m_rangeEnd[1] == 0 ||
          m_rangeBeg[2] == 0 || m_rangeEnd[2] == 0) {
        return 0;
      }

      size_t cell_count = 1;
      for (int i = 0; i < 3; i++) {
        auto diff = std::abs(m_rangeEnd[i] - m_rangeBeg[i]);
        cell_count *= ((diff == 0) ? 1 : diff);
      }
      return cell_count;
    }

    bool is_active() const
    {
      return (m_rangeBeg[0] != 0 || m_rangeEnd[0] != 0 || m_rangeBeg[1] != 0 ||
              m_rangeEnd[1] != 0 || m_rangeBeg[2] != 0 || m_rangeEnd[2] != 0);
    }

    std::string m_bcName;

    // The "original" owner range -- that is, is has not been subsetted
    // due to block decompositions in a parallel run.  It should be the same on
    // all processors...  Primarily used to make parallel collective output easier...
    std::array<int, 6> m_ownerRange{};

    // These are potentially subsetted due to parallel decompositions...
    Ioss::IJK_t m_rangeBeg;
    Ioss::IJK_t m_rangeEnd;

    friend std::ostream &operator<<(std::ostream &os, const BoundaryCondition &bc);
  };

  class DatabaseIO;

  /** \brief A structured zone -- i,j,k
   */
  class StructuredBlock : public EntityBlock
  {
  public:
    StructuredBlock(DatabaseIO *io_database, const std::string &my_name, int index_dim, int ni,
                    int nj = 0, int nk = 0, int off_i = 0, int off_j = 0, int off_k = 0);
    StructuredBlock(DatabaseIO *io_database, const std::string &my_name, int index_dim,
                    Ioss::IJK_t &ordinal, Ioss::IJK_t &offset, Ioss::IJK_t &global_ordinal);

    StructuredBlock *clone(DatabaseIO *database) const;

    ~StructuredBlock() override;

    std::string type_string() const override { return "StructuredBlock"; }
    std::string short_type_string() const override { return "structuredblock"; }
    EntityType  type() const override { return STRUCTUREDBLOCK; }

    const Ioss::NodeBlock &get_node_block() const { return m_nodeBlock; }

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string &my_name) const override;

    AxisAlignedBoundingBox get_bounding_box() const;

    /** \brief Set the 'offset' for the block.
     *
     *  The 'offset' is used to map a cell or node location within a
     *  structured block to the model implicit cell or node location
     *  on a single processor.  zero-based.
     *
     *  The 'global' offsets do the same except for they apply over
     *  the entire model on all processors. zero-based.
     *
     *  For example, the file descriptor (1-based) of
     *  the 37th cell in the 4th block is calculated by:
     *
     *  file_descriptor = offset of block 4 + 37
     *
     *  This can also be used to determine which structured block
     *  a cell with a file_descriptor maps into. An particular
     *  structured block contains all cells in the range:
     *
     *  offset < file_descriptor <= offset+number_cells_per_block
     *
     *  Note that for nodes, the nodeOffset does not take into account
     *  the nodes that are shared between blocks.
     */
    void set_node_offset(size_t offset) { m_nodeOffset = offset; }
    void set_cell_offset(size_t offset) { m_cellOffset = offset; }
    void set_node_global_offset(size_t offset) { m_nodeGlobalOffset = offset; }
    void set_cell_global_offset(size_t offset) { m_cellGlobalOffset = offset; }

    size_t get_node_offset() const { return m_nodeOffset; }
    size_t get_cell_offset() const { return m_cellOffset; }
    size_t get_node_global_offset() const { return m_nodeGlobalOffset; }
    size_t get_cell_global_offset() const { return m_cellGlobalOffset; }

    // Get the local (relative to this block on this processor) node
    // id at the specified i,j,k location (1 <= i,j,k <= ni+1,nj+1,nk+1).  1-based.
    size_t get_block_local_node_id(int ii, int jj, int kk) const
    {
      auto i = ii - m_offsetI;
      auto j = jj - m_offsetJ;
      auto k = kk - m_offsetK;
      assert(i > 0 && i <= m_ni + 1 && j > 0 && j <= m_nj + 1 && k > 0 && k <= m_nk + 1);
      return static_cast<size_t>(k - 1) * (m_ni + 1) * (m_nj + 1) +
             static_cast<size_t>(j - 1) * (m_ni + 1) + i;
    }

    // Get the local (relative to this block on this processor) cell
    // id at the specified i,j,k location (1 <= i,j,k <= ni,nj,nk).  1-based.
    size_t get_block_local_cell_id(int ii, int jj, int kk) const
    {
      auto i = ii - m_offsetI;
      auto j = jj - m_offsetJ;
      auto k = kk - m_offsetK;
      assert(i > 0 && i <= m_ni + 1 && j > 0 && j <= m_nj + 1 && k > 0 && k <= m_nk + 1);
      return static_cast<size_t>(k - 1) * m_ni * m_nj + static_cast<size_t>(j - 1) * m_ni + i;
    }

    // Get the global (over all processors) cell
    // id at the specified i,j,k location (1 <= i,j,k <= ni,nj,nk).  1-based.
    size_t get_global_cell_id(int i, int j, int k) const
    {
      return m_cellGlobalOffset + static_cast<size_t>(k - 1) * m_niGlobal * m_njGlobal +
             static_cast<size_t>(j - 1) * m_niGlobal + i;
    }

    // Get the global (over all processors) node
    // offset at the specified i,j,k location (1 <= i,j,k <= ni,nj,nk).  0-based, does not account
    // for shared nodes.
    size_t get_global_node_offset(int i, int j, int k) const
    {
      return m_nodeGlobalOffset + static_cast<size_t>(k - 1) * (m_niGlobal + 1) * (m_njGlobal + 1) +
             static_cast<size_t>(j - 1) * (m_niGlobal + 1) + i - 1;
    }

    // Get the local (relative to this block on this processor) node id at the specified
    // i,j,k location (1 <= i,j,k <= ni+1,nj+1,nk+1).  0-based.
    size_t get_block_local_node_offset(int ii, int jj, int kk) const
    {
      auto i = ii - m_offsetI;
      auto j = jj - m_offsetJ;
      auto k = kk - m_offsetK;
      assert(i > 0 && i <= m_ni + 1 && j > 0 && j <= m_nj + 1 && k > 0 && k <= m_nk + 1);
      return static_cast<size_t>(k - 1) * (m_ni + 1) * (m_nj + 1) +
             static_cast<size_t>(j - 1) * (m_ni + 1) + i - 1;
    }

    // Get the local (on this processor) cell-node offset at the specified
    // i,j,k location (1 <= i,j,k <= ni+1,nj+1,nk+1).  0-based.
    size_t get_local_node_offset(int i, int j, int k) const

    {
      return get_block_local_node_offset(i, j, k) + m_nodeOffset;
    }

    // Get the global node id at the specified
    // i,j,k location (1 <= i,j,k <= ni+1,nj+1,nk+1).  1-based.
    size_t get_global_node_id(int i, int j, int k) const
    {
      return get_global_node_offset(i, j, k) + 1;
    }

    std::vector<int> get_cell_node_ids(bool add_offset) const
    {
      size_t           node_count = get_property("node_count").get_int();
      std::vector<int> ids(node_count);
      get_cell_node_ids(ids.data(), add_offset);
      return ids;
    }

    template <typename INT> size_t get_cell_node_ids(INT *idata, bool add_offset) const
    {
      // Fill 'idata' with the cell node ids which are the
      // 1-based location of each node in this zone
      // The location is based on the "model" (all processors) zone.
      // If this is a parallel decomposed model, then
      // this block may be a subset of the "model" zone
      //
      // if 'add_offset' is true, then add the m_cellGlobalOffset
      // which changes the location to be the location in the
      // entire "mesh" instead of within a "zone" (all processors)

      size_t index  = 0;
      size_t offset = add_offset ? m_nodeGlobalOffset : 0;

      if (m_nk == 0 && m_nj == 0 && m_ni == 0) {
        return index;
      }

      for (int kk = 0; kk < m_nk + 1; kk++) {
        size_t k = m_offsetK + kk;
        for (int jj = 0; jj < m_nj + 1; jj++) {
          size_t j = m_offsetJ + jj;
          for (int ii = 0; ii < m_ni + 1; ii++) {
            size_t i = m_offsetI + ii;

            size_t ind = k * (m_niGlobal + 1) * (m_njGlobal + 1) + j * (m_niGlobal + 1) + i;

            idata[index++] = ind + offset + 1;
          }
        }
      }

      for (auto idx_id : m_globalIdMap) {
        idata[idx_id.first] = idx_id.second;
      }

      return index;
    }

    template <typename INT> size_t get_cell_ids(INT *idata, bool add_offset) const
    {
      // Fill 'idata' with the cell ids which are the
      // 1-based location of each cell in this zone
      // The location is based on the "model" zone.
      // If this is a parallel decomposed model, then
      // this block may be a subset of the "model" zone
      //
      // if 'add_offset' is true, then add the m_cellGlobalOffset
      // which changes the location to be the location in the
      // entire "mesh" instead of within a "zone"

      size_t index  = 0;
      size_t offset = add_offset ? m_cellGlobalOffset : 0;

      if (m_nk == 0 && m_nj == 0 && m_ni == 0) {
        return index;
      }

      for (int kk = 0; kk < m_nk; kk++) {
        size_t k = m_offsetK + kk;
        for (int jj = 0; jj < m_nj; jj++) {
          size_t j = m_offsetJ + jj;
          for (int ii = 0; ii < m_ni; ii++) {
            size_t i = m_offsetI + ii;

            size_t ind = k * m_niGlobal * m_njGlobal + j * m_niGlobal + i;

            idata[index++] = ind + offset + 1;
          }
        }
      }
      return index;
    }

    bool contains(size_t global_offset) const
    {
      return (global_offset >= m_nodeOffset &&
              global_offset < m_nodeOffset + get_property("node_count").get_int());
    }

  protected:
    int64_t internal_get_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

    int64_t internal_put_field_data(const Field &field, void *data,
                                    size_t data_size) const override;

  private:
    void add_properties_and_fields(int index_dim);

    int m_ni;
    int m_nj;
    int m_nk;

    int m_offsetI; // Valid 'i' ordinal runs from m_offsetI+1 to m_offsetI+m_ni
    int m_offsetJ;
    int m_offsetK;

    int m_niGlobal; // The ni,nj,nk of the master block this is a subset of.
    int m_njGlobal;
    int m_nkGlobal;

    size_t m_nodeOffset;
    size_t m_cellOffset;

    size_t m_nodeGlobalOffset;
    size_t m_cellGlobalOffset;

    Ioss::NodeBlock m_nodeBlock;

  public:
    std::vector<ZoneConnectivity>  m_zoneConnectivity;
    std::vector<BoundaryCondition> m_boundaryConditions;
    std::vector<size_t>            m_blockLocalNodeIndex;
    std::vector<std::pair<size_t, size_t>> m_globalIdMap;
    std::vector<std::pair<size_t, size_t>> m_sharedNode;
  };
} // namespace Ioss
#endif
