// Copyright(C) 1999-2017 National Technology & Engineering Solutions
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

#ifndef IOSS_Iocgns_ParallelDatabaseIO_h
#define IOSS_Iocgns_ParallelDatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>    // for DatabaseUsage
#include <Ioss_DatabaseIO.h> // for DatabaseIO
#include <Ioss_IOFactory.h>  // for IOFactory
#include <Ioss_Map.h>        // for Map
#include <Ioss_State.h>      // for State
#include <iostream>          // for ostream
#include <memory>
#include <stddef.h> // for size_t
#include <stdint.h> // for int64_t
#include <string>   // for string

#include <cgns/Iocgns_DecompositionData.h>

#include <cgnslib.h>

namespace Ioss {
  class CommSet;
  class EdgeBlock;
  class EdgeSet;
  class ElementBlock;
  class ElementSet;
  class ElementTopology;
  class FaceBlock;
  class FaceSet;
  class Field;
  class GroupingEntity;
  class NodeBlock;
  class NodeSet;
  class Region;
  class SideBlock;
  class SideSet;
  class EntityBlock;
} // namespace Ioss

namespace Iocgns {

  using CGNSIntVector = std::vector<cgsize_t>;

  class ParallelDatabaseIO : public Ioss::DatabaseIO
  {
  public:
    enum class entity_type { NODE, ELEM };

    ParallelDatabaseIO(Ioss::Region *region, const std::string &filename,
                       Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                       const Ioss::PropertyManager &props);

    ~ParallelDatabaseIO();

    // Check capabilities of input/output database...  Returns an
    // unsigned int with the supported Ioss::EntityTypes or'ed
    // together. If "return_value & Ioss::EntityType" is set, then the
    // database supports that type (e.g. return_value & Ioss::FACESET)
    unsigned entity_field_support() const override;

    int64_t node_global_to_local__(int64_t global, bool must_exist) const override;
    int64_t element_global_to_local__(int64_t global) const override;

    void release_memory__() override;

    void openDatabase__() const override;
    void closeDatabase__() const override;
    bool node_major() const override { return false; }

    bool needs_shared_node_information() const override { return false; }

    // This isn't quite true since a CGNS library with cgsize_t == 64-bits can read
    // a file with 32-bit ints. However,...
    int int_byte_size_db() const override { return CG_SIZEOF_SIZE; }

    bool begin__(Ioss::State state) override;
    bool end__(Ioss::State state) override;

    bool begin_state__(Ioss::Region *region, int state, double time) override;
    bool end_state__(Ioss::Region *region, int state, double time) override;

    // Metadata-related functions.
    void read_meta_data__() override;
    void write_meta_data();
    void write_results_meta_data();

  private:
    void    handle_structured_blocks();
    void    handle_unstructured_blocks();
    size_t  finalize_structured_blocks();
    int64_t handle_node_ids(void *ids, int64_t num_to_get) const;
    void    finalize_database() override;
    void    get_step_times__() override;
    void    write_adjacency_data();

    int64_t get_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override;
    int64_t get_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    int64_t put_field_internal(const Ioss::Region *reg, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceBlock *nb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                               void *data, size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::EdgeSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::FaceSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::ElementSet *ns, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;
    int64_t put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field, void *data,
                               size_t data_size) const override;

    // ID Mapping functions.
    const Ioss::Map &get_map(entity_type type) const;
    const Ioss::Map &get_map(Ioss::Map &entity_map, int64_t entityCount, int64_t file_offset,
                             int64_t file_count, entity_type type) const;

    int64_t handle_element_ids(const Ioss::ElementBlock *eb, void *ids, size_t num_to_get,
                               size_t offset, size_t count) const;

    // Bulk Data
    void resolve_zone_shared_nodes(const CGNSIntVector &nodes, CGNSIntVector &connectivity_map,
                                   size_t &owned_node_count, size_t &owned_node_offset) const;

    std::vector<int64_t> get_processor_zone_node_offset() const;

    mutable int   cgnsFilePtr{-1};
    mutable int   cgnsSerFilePtr{-1};
    CG_ZoneType_t m_zoneType{CG_ZoneTypeNull};

    mutable std::unique_ptr<DecompositionDataBase> decomp;

    int m_currentVertexSolutionIndex     = 0;
    int m_currentCellCenterSolutionIndex = 0;

    mutable std::vector<size_t> m_zoneOffset; // Offset for local zone/block element ids to global.

    mutable std::vector<size_t>
                                       m_bcOffset; // The BC Section element offsets in unstructured output.
    mutable std::vector<double>        m_timesteps; // Should be able to get this from region?
    std::map<std::string, int>         m_zoneNameMap;
    mutable std::map<int, Ioss::Map *> m_globalToBlockLocalNodeMap;
    mutable CGNSIntVector
        m_elemGlobalImplicitMap; // Position of this element in the global-implicit ordering
  };
} // namespace Iocgns
#endif
