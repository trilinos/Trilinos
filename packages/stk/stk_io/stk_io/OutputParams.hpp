// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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
#include "Ioss_EntityType.h"             // for EntityType
#include "Ioss_GroupingEntity.h"
#include "Ioss_Region.h"
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
#include "stk_mesh/base/FieldBase.hpp"  // for FieldState
#include "stk_mesh/base/Part.hpp"        // for Part
#include "stk_mesh/base/MetaData.hpp"

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


struct OutputParams
{
public:
    OutputParams(Ioss::Region & region, const mesh::BulkData &bulk) :
        m_ioRegion(region),
        m_bulkData(bulk) {}

    Ioss::Region &io_region() const {return m_ioRegion;}
    const mesh::BulkData &bulk_data() const {return m_bulkData;}

    const stk::mesh::Selector *get_subset_selector() const {return m_subsetSelector;}
    void set_subset_selector(const stk::mesh::Selector *subset_selector) {m_subsetSelector = subset_selector;}
    bool has_subset_selector() const {return nullptr != m_subsetSelector;}

    const stk::mesh::Selector *get_shared_selector() const {return m_sharedSelector;}
    void set_shared_selector(const stk::mesh::Selector *shared_selector) {m_sharedSelector = shared_selector;}
    bool has_shared_selector() const {return nullptr != m_sharedSelector;}

    const stk::mesh::Selector *get_output_selector() const {return m_outputSelector;}
    void set_output_selector(const stk::mesh::Selector *output_selector) {m_outputSelector = output_selector;}
    bool has_output_selector() const {return nullptr != m_outputSelector;}

    bool get_sort_stk_parts_by_name() const {return m_sortStkPartsByName;}
    void set_sort_stk_parts_by_name(const bool sortStkPartsByName) {m_sortStkPartsByName = sortStkPartsByName;}

    bool get_use_nodeset_for_block_node_fields() const {return m_useNodesetForBlockNodeFields;}
    void set_use_nodeset_for_block_node_fields(const bool useNodesetForBlockNodeFields) {m_useNodesetForBlockNodeFields = useNodesetForBlockNodeFields;}

    bool get_use_nodeset_for_sideset_node_fields() const {return m_useNodesetForSidesetNodeFields;}
    void set_use_nodeset_for_sideset_node_fields(const bool useNodesetForSidesetNodeFields) {m_useNodesetForSidesetNodeFields = useNodesetForSidesetNodeFields;}

    bool check_field_existence_when_creating_nodesets() const {return m_checkFieldExistenceWhenCreatingNodesets;}
    void check_field_existence_when_creating_nodesets(const bool checkFieldExistenceWhenCreatingNodesets) {m_checkFieldExistenceWhenCreatingNodesets = checkFieldExistenceWhenCreatingNodesets;}

    bool get_use_part_id_for_output() const {return m_usePartIdForOutput;}
    void set_use_part_id_for_output(const bool usePartIdForOutput) {m_usePartIdForOutput = usePartIdForOutput;}

private:
    Ioss::Region & m_ioRegion;
    const mesh::BulkData &m_bulkData;
    const stk::mesh::Selector *m_subsetSelector = nullptr;
    const stk::mesh::Selector *m_sharedSelector = nullptr;
    const stk::mesh::Selector *m_outputSelector = nullptr;
    bool m_sortStkPartsByName = false;
    bool m_useNodesetForBlockNodeFields = true;
    bool m_useNodesetForSidesetNodeFields = true;
    bool m_checkFieldExistenceWhenCreatingNodesets = true;
    bool m_usePartIdForOutput = true;

private:
    OutputParams();
    OutputParams(const OutputParams &);
};

}//namespace io
}//namespace stk
#endif
