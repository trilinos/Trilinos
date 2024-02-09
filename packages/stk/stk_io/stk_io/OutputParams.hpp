// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
//

#ifndef stk_io_OutputParams_hpp
#define stk_io_OutputParams_hpp

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_DBUsage.h>                // for DatabaseUsage
#include <Ioss_Field.h>                  // for Field, Field::RoleType, etc
#include <stddef.h>                      // for size_t
#include <stk_mesh/base/Types.hpp>       // for EntityId, EntityRank
#include <stk_topology/topology.hpp>     // for topology
#include <string>                        // for string, operator<, etc
#include <utility>                       // for pair
#include <vector>                        // for vector
#include <unordered_map>
#include "Ioss_EntityType.h"             // for EntityType
#include "Ioss_GroupingEntity.h"
#include "Ioss_Region.h"
#include "Ioss_Utils.h"
#include "Ioss_DatabaseIO.h"
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
#include "stk_mesh/base/FieldBase.hpp"  // for FieldState
#include "stk_mesh/base/Part.hpp"        // for Part
#include "stk_mesh/base/MetaData.hpp"
#include "stk_io/MeshField.hpp"
#include "stk_io/FieldAndName.hpp"

#include <stk_mesh/base/GetEntities.hpp>            // for count_selected_en...
#include <stk_util/parallel/ParallelReduce.hpp>     // for all_reduce_sum

namespace Ioss { class ElementTopology; }
namespace Ioss { class EntityBlock; }
namespace Ioss { class GroupingEntity; }
namespace Ioss { class Region; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class FieldRestriction; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class Selector; } }

namespace stk {
namespace io {

bool is_part_element_block_io_part(const stk::mesh::Part &part);
bool is_part_assembly_io_part(const stk::mesh::Part &part);

struct OutputParams
{
public:
    OutputParams(Ioss::Region & region, const mesh::BulkData &bulk) :
        m_ioRegion(&region),
        m_bulkData(bulk)
    {
        initialize_output_selectors();
        initialize_block_sizes();
    }

    OutputParams(const mesh::BulkData &bulk) :
        m_ioRegion(nullptr),
        m_bulkData(bulk)
    {
        initialize_output_selectors();
        initialize_block_sizes();
    }

    void set_io_region(Ioss::Region* region) {
      m_ioRegion = region;
    }
    Ioss::Region *io_region_ptr() const {
        return m_ioRegion;
    }
    Ioss::Region &io_region() const {
        STK_ThrowRequireMsg(m_ioRegion != nullptr, "Region is null"); return *m_ioRegion;
    }
    const mesh::BulkData &bulk_data() const {
        return m_bulkData;
    }

    const stk::mesh::Selector *get_subset_selector() const {
        return m_subsetSelector;
    }
    void set_subset_selector(const stk::mesh::Selector *subset_selector) {
        m_subsetSelector = subset_selector;
    }
    bool has_subset_selector() const {
        return nullptr != m_subsetSelector;
    }

    const stk::mesh::Selector *get_shared_selector() const {
        return m_sharedSelector;
    }
    void set_shared_selector(const stk::mesh::Selector *shared_selector) {
        m_sharedSelector = shared_selector;
    }
    bool has_shared_selector() const {
        return nullptr != m_sharedSelector;
    }

    const stk::mesh::Selector *get_output_selector(stk::topology::rank_t rank) const {
        return is_valid_rank(rank) ? m_outputSelector[rank] : nullptr;
    }
    void set_output_selector(stk::topology::rank_t rank, const stk::mesh::Selector *output_selector) {
        if(is_valid_rank(rank)) {
            m_outputSelector[rank] = output_selector;
        }
    }
    bool has_output_selector(stk::topology::rank_t rank) const {
        return is_valid_rank(rank) ? (nullptr != m_outputSelector[rank]) : false;
    }

    const stk::mesh::Selector *get_skin_mesh_selector() const {
        return m_skinMeshSelector;
    }
    void set_skin_mesh_selector(const stk::mesh::Selector *skin_mesh_selector) {
        m_skinMeshSelector = skin_mesh_selector;
    }
    bool has_skin_mesh_selector() const {
        return nullptr != m_skinMeshSelector;
    }

    bool get_sort_stk_parts_by_name() const {
        return m_sortStkPartsByName;
    }
    void set_sort_stk_parts_by_name(const bool sortStkPartsByName) {
        m_sortStkPartsByName = sortStkPartsByName;
    }

    bool get_use_nodeset_for_block_node_fields() const {
        return m_useNodesetForBlockNodeFields;
    }
    void set_use_nodeset_for_block_node_fields(const bool useNodesetForBlockNodeFields) {
        m_useNodesetForBlockNodeFields = useNodesetForBlockNodeFields;
    }

    bool get_use_nodeset_for_sideset_node_fields() const {
        return m_useNodesetForSidesetNodeFields;
    }
    void set_use_nodeset_for_sideset_node_fields(const bool useNodesetForSidesetNodeFields) {
        m_useNodesetForSidesetNodeFields = useNodesetForSidesetNodeFields;
    }

    bool check_field_existence_when_creating_nodesets() const {
        return m_checkFieldExistenceWhenCreatingNodesets;
    }
    void check_field_existence_when_creating_nodesets(const bool checkFieldExistenceWhenCreatingNodesets) {
        m_checkFieldExistenceWhenCreatingNodesets = checkFieldExistenceWhenCreatingNodesets;
    }

    bool get_use_part_id_for_output() const {
        return m_usePartIdForOutput;
    }
    void set_use_part_id_for_output(const bool usePartIdForOutput) {
        m_usePartIdForOutput = usePartIdForOutput;
    }

    bool get_has_ghosting() const {
        return m_hasGhosting;
    }
    void set_has_ghosting(const bool hasGhosting) {
        m_hasGhosting = hasGhosting;
    }

    bool get_has_adaptivity() const {
        return m_hasAdaptivity;
    }
    void set_has_adaptivity(const bool hasAdaptivity) {
        m_hasAdaptivity = hasAdaptivity;
    }

    bool get_is_restart() const {
        return m_isRestart;
    }
    void set_is_restart(const bool restart) {
        m_isRestart = restart;
    }

    bool get_enable_edge_io() const {
        return m_enableEdgeIO;
    }
    void set_enable_edge_io(const bool enableEdgeIO) {
        m_enableEdgeIO = enableEdgeIO;
    }

    const std::vector<stk::io::FieldAndName>& get_additional_attribute_fields() const {
        return m_additionalAttributeFields;
    }
    void set_additional_attribute_fields(const std::vector<stk::io::FieldAndName>& additionalAttributeFields) {
        m_additionalAttributeFields = additionalAttributeFields;
    }

    void set_filter_empty_entity_blocks(const bool filterEmptyEntityBlocks) {
      m_filterEmptyEntityBlocks = filterEmptyEntityBlocks;
    }

    bool get_filter_empty_entity_blocks() const {
      return m_filterEmptyEntityBlocks;
    }

    void set_filter_empty_assembly_entity_blocks(const bool filterEmptyAssemblyEntityBlocks) {
      m_filterEmptyAssemblyEntityBlocks = filterEmptyAssemblyEntityBlocks;
    }

    bool get_filter_empty_assembly_entity_blocks() const {
      return m_filterEmptyAssemblyEntityBlocks;
    }

    const std::unordered_map<unsigned, size_t>& get_block_sizes() const {
      return m_blockSizes;
    }
private:
    OutputParams();
    OutputParams(const OutputParams &);

    void initialize_output_selectors()
    {
        for (unsigned rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEM_RANK; rank++) {
            m_outputSelector[rank] = nullptr;
        }
    }

    void initialize_block_sizes()
    {
      const stk::mesh::MetaData & meta = m_bulkData.mesh_meta_data();
      const mesh::PartVector & parts = meta.get_parts();
      stk::mesh::ConstPartVector elementParts;
      elementParts.reserve(parts.size());

      bool hasAssembly = false;

      for (const stk::mesh::Part * part : parts) {
        if (is_part_element_block_io_part(*part)) {
          elementParts.push_back(part);
        }

        if (is_part_assembly_io_part(*part)) {
          hasAssembly = true;
        }
      }

      if(!hasAssembly) return;

      size_t length = elementParts.size();
      std::vector<size_t> localBlockSizes(length, 0);
      std::vector<size_t> globalBlockSizes(length, 0);

      for(size_t i=0; i<length; ++i) {
        const stk::mesh::Part *part = elementParts[i];
        localBlockSizes[i] = stk::mesh::count_entities(m_bulkData, stk::topology::ELEMENT_RANK, *part);
      }

      stk::all_reduce_sum( m_bulkData.parallel(), localBlockSizes.data(), globalBlockSizes.data(), length);

      for(size_t i=0; i<length; ++i) {
        const stk::mesh::Part *part = elementParts[i];
        m_blockSizes[part->mesh_meta_data_ordinal()] = globalBlockSizes[i];
      }
    }

    bool is_valid_rank(stk::topology::rank_t rank) const {return ((rank >= stk::topology::NODE_RANK) && (rank <= stk::topology::ELEM_RANK)); }

    Ioss::Region * m_ioRegion = nullptr;
    const mesh::BulkData &m_bulkData;
    const stk::mesh::Selector *m_subsetSelector = nullptr;
    const stk::mesh::Selector *m_sharedSelector = nullptr;
    const stk::mesh::Selector *m_outputSelector[stk::topology::ELEM_RANK+1];
    const stk::mesh::Selector *m_skinMeshSelector = nullptr;
    bool m_sortStkPartsByName = false;
    bool m_useNodesetForBlockNodeFields = true;
    bool m_useNodesetForSidesetNodeFields = true;
    bool m_checkFieldExistenceWhenCreatingNodesets = true;
    bool m_usePartIdForOutput = true;
    bool m_hasGhosting = false;
    bool m_hasAdaptivity = false;
    bool m_isRestart = false;
    bool m_enableEdgeIO = false;
    std::vector<stk::io::FieldAndName> m_additionalAttributeFields;
    std::unordered_map<unsigned, size_t> m_blockSizes;

    bool m_filterEmptyEntityBlocks = false;
    bool m_filterEmptyAssemblyEntityBlocks = false;

};

}//namespace io
}//namespace stk
#endif
