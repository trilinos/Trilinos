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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_io/IOHelpers.hpp>
#include <cassert>                             // for assert
#include <cstdint>                             // for int64_t
#include <iostream>                            // for operator<<, basic_ostream
#include <map>                                 // for _Rb_tree_iterator
#include <memory>                              // for allocator_traits<>::va...
#include <stdexcept>                           // for runtime_error
#include <stk_io/IossBridge.hpp>               // for field_data_from_ioss
#include <stk_mesh/base/BulkData.hpp>          // for BulkData
#include <stk_mesh/base/FEMHelpers.hpp>        // for get_side_entity_for_el...
#include <stk_mesh/base/Field.hpp>             // for Field
#include <stk_mesh/base/GetEntities.hpp>       // for get_entities
#include <stk_mesh/base/MetaData.hpp>          // for MetaData
#include <stk_util/util/ReportHandler.hpp>     // for ThrowErrorMsgIf, Throw...
#include <utility>                             // for pair
#include "Ioss_DatabaseIO.h"                   // for DatabaseIO
#include "Ioss_ElementBlock.h"                 // for ElementBlock
#include "Ioss_ElementTopology.h"              // for ElementTopology
#include "Ioss_EntityType.h"                   // for NODEBLOCK, SIDESET
#include "Ioss_Field.h"                        // for Field, Field::ATTRIBUTE
#include "Ioss_GroupingEntity.h"               // for GroupingEntity
#include "Ioss_IOFactory.h"                    // for NameList
#include "Ioss_NodeBlock.h"                    // for NodeBlock
#include "Ioss_NodeSet.h"                      // for NodeSet
#include "Ioss_Property.h"                     // for Property
#include "Ioss_Region.h"                       // for Region, ElementBlockCo...
#include "Ioss_SideBlock.h"                    // for SideBlock
#include "Ioss_SideSet.h"                      // for SideSet
#include "Ioss_State.h"                        // for STATE_TRANSIENT
#include "Ioss_VariableType.h"                 // for VariableType
#include "ProcessSetsOrBlocks.hpp"             // for get_part_for_grouping_...
#include "StkIoUtils.hpp"                      // for part_primary_entity_rank
#include "stk_io/OutputParams.hpp"             // for OutputParams
#include "stk_mesh/base/Bucket.hpp"            // for Bucket
#include "stk_mesh/base/Entity.hpp"            // for Entity
#include "stk_mesh/base/EntityKey.hpp"         // for EntityKey, operator<<
#include "stk_mesh/base/FieldBase.hpp"         // for FieldBase, field_bytes...
#include "stk_mesh/base/FieldState.hpp"        // for FieldState
#include "stk_mesh/base/Part.hpp"              // for Part
#include "stk_mesh/base/Selector.hpp"          // for Selector, operator&
#include "stk_mesh/base/Types.hpp"             // for EntityRank, PartVector
#include "stk_topology/topology.hpp"           // for topology, topology::NO...
#include "stk_util/parallel/CommSparse.hpp"    // for CommSparse, pack_and_c...
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/util/ParameterList.hpp"     // for Par...

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk {
namespace io {

void write_defined_global_any_fields(std::shared_ptr<Ioss::Region> region,
                                     std::vector<stk::io::GlobalAnyVariable> &global_any_fields)
{
    for (size_t i=0; i < global_any_fields.size(); i++) {
        const std::string &name = global_any_fields[i].m_name;
        const std::any *value = global_any_fields[i].m_value;
        stk::util::ParameterType::Type type = global_any_fields[i].m_type;
        impl::write_parameter(region, name, *value, type);
    }
}

size_t get_entities(const stk::mesh::Part &part,
                    const stk::mesh::BulkData &bulk,
                    std::vector<stk::mesh::Entity> &entities,
                    const stk::mesh::Selector *subset_selector)
{
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    stk::mesh::EntityRank type = stk::io::part_primary_entity_rank(part);

    stk::mesh::Selector own = meta.locally_owned_part();
    stk::mesh::Selector selector = part & own;
    if (subset_selector) selector &= *subset_selector;

    const bool sortById = true;
    stk::mesh::get_entities(bulk, type, selector, entities, sortById);
    return entities.size();
}

void pack_distribution_factor(const stk::mesh::BulkData& bulk,
                              CommBuffer& buf,
                              const stk::mesh::FieldBase* distFact,
                              stk::mesh::Entity side)
{
    stk::mesh::EntityKey key = bulk.entity_key(side);
    buf.pack<stk::mesh::EntityKey>(key);

    const stk::mesh::Bucket & bucket = bulk.bucket(side);
    const unsigned size = field_bytes_per_entity(*distFact, bucket);
    if (size) {
        unsigned char * const ptr = reinterpret_cast<unsigned char *>(stk::mesh::field_data(*distFact , side));
        buf.pack<unsigned char>(ptr , size);
    }
}

void unpack_distribution_factor(const stk::mesh::BulkData& bulk,
                                CommBuffer& buf,
                                const stk::mesh::FieldBase* distFact)
{
    stk::mesh::EntityKey key;
    buf.unpack<stk::mesh::EntityKey>(key);

    stk::mesh::Entity side = bulk.get_entity(key);
    STK_ThrowRequire(bulk.is_valid(side));
    const stk::mesh::Bucket & bucket = bulk.bucket(side);
    const unsigned size = field_bytes_per_entity( *distFact, bucket );
    if (size)
    {
        unsigned char * ptr = reinterpret_cast<unsigned char *>( stk::mesh::field_data(*distFact , side));
        buf.unpack<unsigned char>(ptr , size);
    }
}

void communicate_shared_side_entity_fields(const stk::mesh::BulkData& bulk,
                                           const stk::mesh::FieldBase* distFact,
                                           std::vector<stk::mesh::Entity>& sides)
{
    stk::CommSparse comm(bulk.parallel());

    const bool anythingToUnpack =
        stk::pack_and_communicate(comm, [&comm, &bulk, &distFact, &sides]() {
           std::vector<int> sharedProcs;
            for (stk::mesh::Entity side : sides) {
                STK_ThrowRequireMsg(bulk.is_valid(side),"communicate_shared_side_entity_fields, invalid side");
                sharedProcs.clear();
                bulk.comm_shared_procs(side, sharedProcs);
                for (int sharedProc : sharedProcs)
                {
                    CommBuffer & buffer = comm.send_buffer(sharedProc);
                    pack_distribution_factor(bulk, buffer, distFact, side);
                }
            }
        });

    if (anythingToUnpack) {
        for (int p = 0 ; p < bulk.parallel_size(); ++p) {
            CommBuffer & buf = comm.recv_buffer(p);
            while (buf.remaining()) {
                unpack_distribution_factor(bulk, buf, distFact);
            }
        }
    }
}

template <typename INT>
void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
    assert(sset->type() == Ioss::SIDESET);

    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    Ioss::Region *region = sset->get_database()->get_region();
    const std::string universalAlias = region->get_alias("universal_sideset", Ioss::SIDESET);
    if (sset->name() == universalAlias)
        return;

    size_t block_count = sset->block_count();
    for (size_t i=0; i < block_count; i++) {
        Ioss::SideBlock *block = sset->get_block(i);
        if (stk::io::include_entity(block)) {

            stk::mesh::Part * sb_part = get_part_for_grouping_entity(*region, meta, block);
            if (sb_part == nullptr){
                sb_part = get_part_for_grouping_entity(*region, meta, sset);
            }
            if (sb_part == nullptr){
                continue;
            }

            // Get topology of the sides being defined to see if they
            // are 'faces' or 'edges'.  This is needed since for shell-type
            // elements, (and actually all elements) a sideset can specify either a face or an edge...
            // For a quad shell, sides 1,2 are faces and 3,4,5,6 are edges.

            // NOTE: This assumes that the sides within a side block are homogenous.  If the side_set
            //       is not split into homogenous side_blocks, then the topology will not necessarily
            //       be the same and this could fail (a sideset of mixed edges and faces)
            int par_dimen = block->topology()->parametric_dimension();

            stk::mesh::EntityRank side_rank = par_dimen == 1 ? stk::topology::EDGE_RANK : stk::topology::FACE_RANK;

            // Would be nice to do:
            //    std::vector<stk::mesh::Entity> sides ;
            //    get_entities(*sb_part, bulk, sides, nullptr);
            // But, we need the entities in the exact same order and count  as they appear on the
            // mesh file.  The get_entities can give them in different order and will ignore duplicated
            // "sides" that occur in some exodus files.

            std::vector<INT> elemSidePairs;
            block->get_field_data("element_side", elemSidePairs);
            size_t side_count = elemSidePairs.size()/2;

            // NOTE: Needed for shell side topology offset translation to IOSS
            std::int64_t sideSetOffset = Ioss::Utils::get_side_offset(block);

            std::vector<stk::mesh::Entity> sides;
            sides.reserve(side_count);

            for (size_t is = 0; is < side_count; ++is) {
                INT elemId = elemSidePairs[is*2];
                INT elemSide = elemSidePairs[is*2+1] + sideSetOffset - 1;
                stk::mesh::Entity const elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elemId);
                if (!bulk.is_valid(elem)) { continue; }
                stk::topology elemTopo = bulk.bucket(elem).topology();
                if(par_dimen == 0)
                {
                  stk::topology faceTopo = elemTopo.sub_topology(elemTopo.side_rank(), elemSide);

                  Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(faceTopo.name(), false);
                  par_dimen = ioss_topo->parametric_dimension();
                  side_rank = par_dimen == 1 ? stk::topology::EDGE_RANK : stk::topology::FACE_RANK;
                }
                STKIORequire(par_dimen == 1 || par_dimen == 2);

                if (par_dimen == 1) {
                  // conversion from face ordinal to edge ordinal for shells
                  if (elemTopo.is_shell()) {
                    elemSide -= elemTopo.num_faces();
                  }
                }
                stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_id_side_pair_of_rank(bulk, elemId, elemSide, side_rank);

                if (bulk.is_valid(side)) {
                    sides.push_back(side);
                }
                else {
                    std::ostringstream os;
                    os<<"P"<<bulk.parallel_rank()<<" STK IO Warning, process_surface_entity_df: side "
                      << elemSide << " on element " << elemId
                      <<" in sideset: " << sset->name()
                      <<" does not have a corresponding face entity in the current mesh database"
                      <<" (possibly referenced by a reduced restart file)."<<std::endl;
                    std::cerr<<os.str();
                }
            }

            const stk::mesh::FieldBase *df_field = stk::io::get_distribution_factor_field(*sb_part);
            if (df_field != nullptr) {
                stk::io::field_data_from_ioss(bulk, df_field, sides, block, "distribution_factors");
                communicate_shared_side_entity_fields(bulk, df_field, sides);
            }

            // Add all attributes as fields.
            // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
            // named components of the 'attribute' field, so add them instead.
            Ioss::NameList names;
            block->field_describe(Ioss::Field::ATTRIBUTE, &names);
            for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
                if(is_skipped_attribute_field(*I, names.size()))
                    continue;
                stk::mesh::FieldBase *field = meta.get_field(side_rank, *I);
                if (field)
                    stk::io::field_data_from_ioss(bulk, field, sides, block, *I);
            }
        }
    }
}

void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
    if (stk::io::db_api_int_size(sset) == 4) {
        process_surface_entity_df<int>(sset, bulk);
    } else {
        process_surface_entity_df<int64_t>(sset, bulk);
    }
}


template <typename INT>
void process_node_coords_and_attributes(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
    // This must be called after the "process_element_blocks" call
    // since there may be nodes that exist in the database that are
    // not part of the analysis mesh due to subsetting of the element
    // blocks.

    // Currently, all nodes found in the finite element mesh are defined
    // as nodes in the stk_mesh database. If some of the element blocks
    // are omitted, then there will be disconnected nodes defined.
    // However, if we only define nodes that are connected to elements,
    // then we risk missing "free" nodes that the user may want to have
    // existing in the model.
    const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
    assert(node_blocks.size() == 1);

    Ioss::NodeBlock *nb = node_blocks[0];

    size_t node_count = nb->get_property("entity_count").get_int();

    std::vector<INT> ids;
    nb->get_field_data("ids", ids);

    std::vector<stk::mesh::Entity> nodes;
    nodes.reserve(node_count);
    for (size_t i=0; i < ids.size(); i++) {
        stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, ids[i]);
        nodes.push_back(node);
    }

    // Temporary (2013/04/02) kluge for Salinas porting to stk-based mesh.
    // Salinas uses the "implicit id" which is the ordering of the nodes
    // in the "undecomposed" or "serial" mesh as user-visible ids
    // instead of the "global" ids. If there exists a stk-field with the
    // name "implicit_node_ids", then populate the field with the correct
    // data.
    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    stk::mesh::FieldBase *implicit_node_id_field = meta.get_field(stk::topology::NODE_RANK, "implicit_node_ids");
    if (implicit_node_id_field) {
        stk::io::field_data_from_ioss(bulk, implicit_node_id_field, nodes, nb, "implicit_ids");
    }


    stk::mesh::FieldBase const* coord_field = meta.coordinate_field();
    stk::io::field_data_from_ioss(bulk, coord_field, nodes, nb, "mesh_model_coordinates");

    // Add all attributes as fields.
    // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
    // named components of the 'attribute' field, so add them instead.
    Ioss::NameList names;
    nb->field_describe(Ioss::Field::ATTRIBUTE, &names);
    for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(is_skipped_attribute_field(*I, names.size()))
            continue;
        stk::mesh::FieldBase *field = meta.get_field(stk::topology::NODE_RANK, *I);
        if (field)
            stk::io::field_data_from_ioss(bulk, field, nodes, nb, *I);
    }
}

template void process_node_coords_and_attributes<int>(Ioss::Region &region, stk::mesh::BulkData &bulk);
template void process_node_coords_and_attributes<int64_t>(Ioss::Region &region, stk::mesh::BulkData &bulk);


template <typename INT>
void process_elem_attributes_and_implicit_ids(Ioss::Region &region, stk::mesh::BulkData &bulk, const bool shouldAutoLoadAttributes)
{
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
            it != elem_blocks.end(); ++it) {
        Ioss::ElementBlock *entity = *it;

        if (stk::io::include_entity(entity)) {
            stk::mesh::Part* part = get_part_for_grouping_entity(region, meta, entity);
            if (part == nullptr)
            {
                continue;
            }

            // See if we need to get the list of elements...
            bool elements_needed = false;
            Ioss::NameList names;
            entity->field_describe(Ioss::Field::ATTRIBUTE, &names);

            stk::mesh::FieldBase *implicit_elem_id_field = meta.get_field(stk::topology::ELEMENT_RANK, "implicit_element_ids");
            if (implicit_elem_id_field) {
                elements_needed = true;
            } else {
                for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
                    if(is_skipped_attribute_field(*I, names.size()))
                        continue;
                    stk::mesh::FieldBase *field = meta.get_field(stk::topology::ELEMENT_RANK, *I);
                    if (field) {
                        elements_needed = true;
                        break;
                    }
                }
            }

            if (!elements_needed)
                continue;

            std::vector<INT> elem_ids ;
            entity->get_field_data("ids", elem_ids);
            size_t element_count = elem_ids.size();

            std::vector<stk::mesh::Entity> elements;
            elements.reserve(element_count);

            for(size_t i=0; i<element_count; ++i) {
                stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elem_ids[i]);
                elements.push_back(elem);
            }

            // Temporary (2013/04/17) kluge for Salinas porting to stk-based mesh.
            // Salinas uses the "implicit id" which is the ordering of the nodes
            // in the "undecomposed" or "serial" mesh as user-visible ids
            // instead of the "global" ids. If there exists a stk-field with the
            // name "implicit_element_ids", then populate the field with the correct
            // data.
            if (implicit_elem_id_field) {
                stk::io::field_data_from_ioss(bulk, implicit_elem_id_field, elements, entity, "implicit_ids");
            }

            // Add all element attributes as fields.
            // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
            // named components of the 'attribute' field, so add them instead.
            if (shouldAutoLoadAttributes)
            {
                for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
                    if(is_skipped_attribute_field(*I, names.size()))
                        continue;
                    stk::mesh::FieldBase *field = meta.get_field(stk::topology::ELEMENT_RANK, *I);
                    if (field)
                        stk::io::field_data_from_ioss(bulk, field, elements, entity, *I);
                }
            }
        }
    }
}

template void process_elem_attributes_and_implicit_ids<int>(Ioss::Region &region, stk::mesh::BulkData &bulk, const bool shouldAutoLoadAttributes);
template void process_elem_attributes_and_implicit_ids<int64_t>(Ioss::Region &region, stk::mesh::BulkData &bulk, const bool shouldAutoLoadAttributes);

stk::mesh::Part* get_part_for_nodeset(const stk::mesh::MetaData &meta, Ioss::NodeSet *entity, NodesetMap &nodesetMap)
{
    const std::string & name = entity->name();
    stk::mesh::Part* part = meta.get_part(name);

    if(nullptr == part) {
        NodesetMap::iterator it = nodesetMap.find(entity);
        if(it != nodesetMap.end())  {
            part = it->second;
        }
    }

    if(nullptr == part) {
        Ioss::Region *region = entity->get_database()->get_region();
        part = get_part_for_grouping_entity(*region, meta, entity);
    }

    return part;
}

// ========================================================================
template <typename INT>
void process_nodesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
    // Should only process nodes that have already been defined via the element
    // blocks connectivity lists.
    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
    const stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    NodesetMap nodesetMap;

    populate_hidden_nodesets(region, meta, nodesetMap);

    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
            it != node_sets.end(); ++it) {
        Ioss::NodeSet *entity = *it;

        if (stk::io::include_entity(entity)) {
            stk::mesh::Part* part = get_part_for_nodeset(meta, entity, nodesetMap);
            stk::mesh::PartVector add_parts;

            if(part != nullptr) {
                add_parts.push_back(part);
            }
            else {
                continue;
            }

            std::vector<INT> node_ids ;
            size_t node_count = entity->get_field_data("ids", node_ids);

            std::vector<stk::mesh::Entity> nodes(node_count);
            for(size_t i=0; i<node_count; ++i) {
                nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, node_ids[i] );
                if (!bulk.is_valid(nodes[i]) ) {
                    nodes[i] = stk::mesh::Entity();
                }
            }

            if(nullptr != part) {
                std::string distributionFactorsPerNodesetFieldName = "distribution_factors_" + part->name();

                stk::mesh::Field<double> *df_field_per_nodeset =
                        meta.get_field<double>(stk::topology::NODE_RANK, distributionFactorsPerNodesetFieldName);

                if (df_field_per_nodeset != nullptr) {
                    stk::io::field_data_from_ioss(bulk, df_field_per_nodeset, nodes, entity, "distribution_factors");
                }
                else {
                    const stk::mesh::FieldBase *dfField = stk::io::get_distribution_factor_field(*part);
                    if (nullptr != dfField)
                    {
                        stk::io::field_data_from_ioss(bulk, dfField, nodes, entity, "distribution_factors");
                    }
                }
            }

            // Add all attributes as fields.
            // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
            // named components of the 'attribute' field, so add them instead.
            Ioss::NameList names;
            entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
            for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
                if(is_skipped_attribute_field(*I, names.size()))
                    continue;
                stk::mesh::FieldBase *field = meta.get_field(stk::topology::NODE_RANK, *I);
                if (field)
                    stk::io::field_data_from_ioss(bulk, field, nodes, entity, *I);
            }
        }
    }
}

template void process_nodesets_df<int>(Ioss::Region &region, stk::mesh::BulkData &bulk);
template void process_nodesets_df<int64_t>(Ioss::Region &region, stk::mesh::BulkData &bulk);

void process_sidesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
    const Ioss::SideSetContainer& side_sets = region.get_sidesets();

    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {
        Ioss::SideSet *entity = *it;

        if (stk::io::include_entity(entity)) {
            process_surface_entity_df(entity, bulk);
        }
    }
}

bool is_processed_field(stk::mesh::Part &part,
                        stk::mesh::EntityRank part_type,
                        Ioss::GroupingEntity *io_entity,
                        const stk::io::FieldAndName &namedField)
{
    const stk::mesh::FieldBase *f = namedField.field();
    const std::string& field_name = namedField.db_name();
    // There is ugliness here to deal with the output of fields on the nodes of a part,
    // either on the nodeblock or on a nodeset part created from the part nodes.

    bool isCorrectFieldRank = (f->entity_rank() == part_type);
    bool forceNodalOutput   = (namedField.m_forceNodeblockOutput && io_entity->type() == Ioss::NODEBLOCK);
    bool isFieldOnPart      = stk::io::is_field_on_part(f, part_type, part);
    bool isFieldOnEntity    = ( io_entity->type() != Ioss::NODEBLOCK && io_entity->field_exists(field_name));

    bool isProcessed = isCorrectFieldRank && (forceNodalOutput || isFieldOnPart || isFieldOnEntity);

    return isProcessed;
}

void put_field_data(OutputParams &params,
                    stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    const std::vector<stk::io::FieldAndName> &namedFields,
                    Ioss::Field::RoleType filter_role,
                    const stk::mesh::FieldState *state)
{
    bool hasFieldsToProcess = false;
    for (const stk::io::FieldAndName& namedField: namedFields) {
        hasFieldsToProcess |= is_processed_field(part, part_type, io_entity, namedField);
    }

    if (hasFieldsToProcess)
    {
        std::vector<stk::mesh::Entity> entities;
        impl::fill_output_entities(io_entity, &part, part_type, params, entities);

        for (const stk::io::FieldAndName& namedField: namedFields) {
            const stk::mesh::FieldBase *f = namedField.field();
            const std::string& field_name = namedField.db_name();

            if (is_processed_field(part, part_type, io_entity, namedField)) {
                if(state != nullptr)
                {
                    f = f->field_state(*state);
                    if(f == nullptr) {
                        f = namedField.field();
                    }
                }

                stk::io::field_data_to_ioss(params.bulk_data(), f, entities, io_entity, field_name, filter_role);
            }
        }
    }
}

void put_field_data(OutputParams &params,
                    stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    const std::vector<stk::io::FieldAndName> &namedFields,
                    const stk::mesh::FieldState *state)
{
    put_field_data(params, part, part_type, io_entity, namedFields, Ioss::Field::Field::TRANSIENT, state);
}

namespace impl {

template <typename DataType>
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           DataType globalVarData)
{
    STK_ThrowErrorMsgIf (output_region == nullptr,
                     "There is no Output mesh region associated with this Mesh Data.");
    STK_ThrowErrorMsgIf (output_region->get_state() != Ioss::STATE_TRANSIENT,
                     "The output region " << output_region->name() <<
                     " is not in the correct state for outputting data at this time.");
    STK_ThrowErrorMsgIf (!output_region->field_exists(globalVarName),
                     "The field named '" << globalVarName << "' does not exist.");
    output_region->put_field_data(globalVarName, &globalVarData, sizeof(DataType));
}

template
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           int globalVarData);
template
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           int64_t globalVarData);
template
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           double globalVarData);

template <typename DataType>
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           std::vector<DataType> &globalVarData)
{
    STK_ThrowErrorMsgIf (output_region == nullptr,
                     "There is no Output mesh region associated with this Mesh Data.");
    STK_ThrowErrorMsgIf (output_region->get_state() != Ioss::STATE_TRANSIENT,
                     "The output region " << output_region->name() <<
                     " is not in the correct state for outputting data at this time.");
    STK_ThrowErrorMsgIf (!output_region->field_exists(globalVarName),
                     "The field named '" << globalVarName << "' does not exist "
                     "on output region "  << output_region->name());
    size_t comp_count = output_region->get_fieldref(globalVarName).raw_storage()->component_count();
    STK_ThrowErrorMsgIf (comp_count != globalVarData.size(),
                     "On output region "  << output_region->name() <<
                     ", the field named '" << globalVarName << "' was registered with size "
                     << comp_count
                     << " but the output size is " << globalVarData.size());

    output_region->put_field_data(globalVarName, globalVarData);
}

template
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           std::vector<int> &globalVarData);
template
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           std::vector<int64_t> &globalVarData);
template
void write_global(std::shared_ptr<Ioss::Region> output_region, const std::string &globalVarName,
                           std::vector<double> &globalVarData);

bool has_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName)
{
  STK_ThrowErrorMsgIf (input_region == nullptr,
                   "There is no Input mesh region associated with this Mesh Data.");

  return input_region->field_exists(globalVarName);
}

template <typename DataType>
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          DataType &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found)
{
    STK_ThrowErrorMsgIf (input_region == nullptr,
                     "There is no Input mesh region associated with this Mesh Data.");

    if (input_region->field_exists(globalVarName)) {
        input_region->get_fieldref(globalVarName).check_type(iossType);
        input_region->get_field_data(globalVarName, &globalVarData, sizeof(DataType));
        return true;
    }
    else {
        if (abort_if_not_found) {
            std::ostringstream msg;
            msg << "ERROR: The field named '" << globalVarName << "' does not exist "
                    << "on input region "  << input_region->name();
            throw std::runtime_error( msg.str() );
        }
        return false;
    }
}

template
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          int &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);
template
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          int64_t &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);
template
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          double &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);

template <typename DataType>
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          std::vector<DataType> &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found)
{
    STK_ThrowErrorMsgIf (input_region == nullptr,
                     "There is no Input mesh region associated with this Mesh Data.");

    if (input_region->field_exists(globalVarName)) {
        input_region->get_fieldref(globalVarName).check_type(iossType);
        input_region->get_field_data(globalVarName, globalVarData);
        return true;
    }
    else {
        if (abort_if_not_found) {
            std::ostringstream msg;
            msg << "ERROR: The field named '" << globalVarName << "' does not exist "
                    << "on input region "  << input_region->name();
            throw std::runtime_error( msg.str() );
        }
        return false;
    }
}

template
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          std::vector<int> &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);
template
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          std::vector<int64_t> &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);
template
bool read_global(std::shared_ptr<Ioss::Region> input_region, const std::string &globalVarName,
                          std::vector<double> &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);


void write_parameter(std::shared_ptr<Ioss::Region> output_region,
                              const std::string &name, const stk::util::Parameter &param)
{
  write_parameter(output_region, name, param.value, param.type);
}

void write_parameter(std::shared_ptr<Ioss::Region> output_region,
                              const std::string &name, const std::any &any_value,
                              stk::util::ParameterType::Type type)
{
    try {
        switch(type)
        {
            case stk::util::ParameterType::INTEGER: {
                int value = std::any_cast<int>(any_value);
                write_global(output_region, name, value);
                break;
            }

            case stk::util::ParameterType::INT64: {
                int64_t value = std::any_cast<int64_t>(any_value);
                write_global(output_region, name, value);
                break;
            }

            case stk::util::ParameterType::DOUBLE: {
                double value = std::any_cast<double>(any_value);
                write_global(output_region, name, value);
                break;
            }

            case stk::util::ParameterType::DOUBLEVECTOR: {
                std::vector<double> vec = std::any_cast<std::vector<double> >(any_value);
                write_global(output_region, name, vec);
                break;
            }

            case stk::util::ParameterType::INTEGERVECTOR: {
                std::vector<int> vec = std::any_cast<std::vector<int> >(any_value);
                write_global(output_region, name, vec);
                break;
            }

            case stk::util::ParameterType::INT64VECTOR: {
                std::vector<int64_t> vec = std::any_cast<std::vector<int64_t> >(any_value);
                write_global(output_region, name, vec);
                break;
            }

            default: {
                std::cerr << "WARNING: '" << name << "' is not a supported type. It's value cannot be output." << std::endl;
                break;
            }
        }
    }
    catch(...) {
        std::cerr << "ERROR: Actual type of parameter named '" << name
                << "' does not match the declared type. Something went wrong."
                << " Maybe you need to call add_global_ref() instead of add_global() or vice-versa.";
        throw;
    }
}

bool read_parameter(std::shared_ptr<Ioss::Region> input_region,
                             const std::string &globalVarName,
                             stk::util::Parameter& param,
                             bool abort_if_not_found)
{
  return read_parameter(input_region, globalVarName, param.value, param.type, abort_if_not_found);
}

bool read_parameter(std::shared_ptr<Ioss::Region> input_region,
                             const std::string &globalVarName,
                             std::any &any_value, stk::util::ParameterType::Type type,
                             bool abort_if_not_found)
{
    bool success = false;
    switch(type)
    {
        case stk::util::ParameterType::INTEGER: {
            int value = 0;
            success = read_global(input_region, globalVarName, value, Ioss::Field::INTEGER,
                                           abort_if_not_found);
            any_value = value;
            break;
        }

        case stk::util::ParameterType::INT64: {
            int64_t value = 0;
            success = read_global(input_region, globalVarName, value, Ioss::Field::INT64,
                                           abort_if_not_found);
            any_value = value;
            break;
        }

        case stk::util::ParameterType::DOUBLE: {
            double value = 0;
            success = read_global(input_region, globalVarName, value, Ioss::Field::REAL,
                                           abort_if_not_found);
            any_value = value;
            break;
        }

        case stk::util::ParameterType::DOUBLEVECTOR: {
            std::vector<double> vec;
            success = read_global(input_region, globalVarName, vec, Ioss::Field::REAL,
                                           abort_if_not_found);
            any_value = vec;
            break;
        }

        case stk::util::ParameterType::INTEGERVECTOR: {
            std::vector<int> vec;
            success = read_global(input_region, globalVarName, vec, Ioss::Field::INTEGER,
                                           abort_if_not_found);
            any_value = vec;
            break;
        }

        case stk::util::ParameterType::INT64VECTOR: {
            std::vector<int64_t> vec;
            success = read_global(input_region, globalVarName, vec, Ioss::Field::INT64,
                                           abort_if_not_found);
            any_value = vec;
            break;
        }

        default: {
            std::cerr << "WARNING: '" << globalVarName << "' is not a supported type. It's value cannot be input." << std::endl;
            break;
        }
    }
    return success;
}

void add_global(std::shared_ptr<Ioss::Region> region,
                         const std::string &globalVarName,
                         const std::string &storage,
                         Ioss::Field::BasicType dataType,
                         int copies,
                         Ioss::Field::RoleType role)
{
    STK_ThrowErrorMsgIf(region->field_exists(globalVarName),
                     "On region named " << region->name() <<
                     " Attempt to add global variable '" << globalVarName << "' twice.");

    region->field_add(Ioss::Field(globalVarName, dataType, storage, copies, role, 1));
}

void add_global(std::shared_ptr<Ioss::Region> region,
                         const std::string &globalVarName,
                         int component_count,
                         Ioss::Field::BasicType dataType,
                         int copies,
                         Ioss::Field::RoleType role)
{
    if (component_count == 1) {
        add_global(region, globalVarName, "scalar", dataType, copies, role);
    } else {
        std::ostringstream type;
        type << "Real[" << component_count << "]";
        add_global(region, globalVarName, type.str(), dataType, copies, role);
    }
}

void fill_output_entities(Ioss::GroupingEntity *io_entity,
                                   stk::mesh::Part *part,
                                   stk::mesh::EntityRank part_type,
                                   OutputParams &params,
                                   std::vector<stk::mesh::Entity> &entities)
{
    if(io_entity->type() == Ioss::SIDEBLOCK)
    {
        // Temporary Kluge to handle sideblocks which contain internally generated sides
        // where the "ids" field on the io_entity doesn't work to get the correct side...
        // NOTE: Could use this method for all entity types, but then need to correctly
        // specify whether shared entities are included/excluded (See IossBridge version).

        Ioss::SideBlock *block = dynamic_cast<Ioss::SideBlock *>(io_entity);
        if(block==nullptr)
        {
            std::ostringstream msg;
            msg << " INTERNAL ERROR: grouping entity failed to a SideBlock. Contact sierra-help@sandia.gov for assistance.\n";
            throw std::runtime_error(msg.str());
        }

        std::vector<int> elem_side_ids;

        fill_data_for_side_block( params, *io_entity, part,
                                  block->parent_element_topology(),
                                  elem_side_ids, entities);

        size_t num_sides = entities.size();
        if(num_sides != static_cast<size_t>(io_entity->get_property("entity_count").get_int()))
        {
            std::ostringstream msg;
            msg << " INTERNAL_ERROR: Number of sides on part " << part->name() << " (" << num_sides
                    << ") does not match number of sides in the associated Ioss SideBlock named "
                    << io_entity->name() << " (" << io_entity->get_property("entity_count").get_int()
                    << ").";
            throw std::runtime_error(msg.str());
        }
    } else {
        stk::io::get_output_entity_list(io_entity, part_type, params, entities);
    }
}

}

} // namespace io
} // namespace stk
