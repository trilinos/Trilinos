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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_util/environment/Env.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ionit_Initializer.h>                       // for Initializer
#include <assert.h>                                  // for assert
#include <stdlib.h>                                  // for exit, etc
#include <string.h>                                  // for memcpy
#include <cstdint>                                   // for int64_t
#include <iostream>                                  // for operator<<, etc
#include <iterator>
#include <limits>                                    // for numeric_limits
#include <map>
#include <stdexcept>                                 // for runtime_error
#include <stk_io/IOHelpers.hpp>
#include <stk_io/OutputFile.hpp>                     // for OutputFile
#include <stk_io/IossBridge.hpp>                     // for FieldAndName, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData, etc
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>                   // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>                // for MetaData, etc
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>    // for ThrowErrorMsgIf, etc
#include <utility>                                   // for pair, make_pair
#include "Ioss_CodeTypes.h"                          // for NameList
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"                         // for DatabaseIO
#include "Ioss_ElementBlock.h"                       // for ElementBlock
#include "Ioss_ElementTopology.h"                    // for ElementTopology
#include "Ioss_EntityType.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"                     // for GroupingEntity
#include "Ioss_IOFactory.h"                          // for IOFactory
#include "Ioss_NodeBlock.h"                          // for NodeBlock
#include "Ioss_NodeSet.h"                            // for NodeSet
#include "Ioss_ParallelUtils.h"                      // for ParallelUtils
#include "Ioss_Property.h"                           // for Property
#include "Ioss_PropertyManager.h"                    // for PropertyManager
#include "Ioss_Region.h"                             // for Region, etc
#include "Ioss_SideBlock.h"                          // for SideBlock
#include "Ioss_SideSet.h"                            // for SideSet
#include "Ioss_State.h"
#include "Ioss_VariableType.h"                       // for VariableType
#include "ProcessSetsOrBlocks.hpp"
#include "SidesetTranslator.hpp"
#include "StkIoUtils.hpp"
#include "Teuchos_RCP.hpp"                           // for RCP::operator->, etc
#include "boost/any.hpp"                             // for any_cast, any
#include "stk_io/DatabasePurpose.hpp"                // for DatabasePurpose, etc
#include "stk_io/MeshField.hpp"                      // for MeshField, etc
#include "stk_io/SidesetUpdater.hpp"
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for FieldBase
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/FieldState.hpp"              // for FieldState
#include "stk_mesh/base/Part.hpp"                    // for Part
#include "stk_mesh/base/Selector.hpp"                // for Selector, etc
#include "stk_mesh/base/Types.hpp"                   // for FieldVector, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/parallel/Parallel.hpp"            // for ParallelMachine, etc
#include "stk_util/util/ParameterList.hpp"           // for Type, etc
#include "stk_util/diag/StringUtil.hpp"           // for Type, etc
#include "stk_util/util/string_case_compare.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace io {

impl::OutputFile::~OutputFile()
{
    if(m_region->get_state() == Ioss::STATE_TRANSIENT)
    {
        m_region->end_mode(Ioss::STATE_TRANSIENT);
    }
    stk::io::delete_selector_property(*m_region);
    delete m_multiStateSuffixes;
}

Ioss::DatabaseIO *impl::OutputFile::get_output_database()
{
    if(m_region.get() == nullptr) {
        return nullptr;
    }

    return m_region->get_database();
}

bool impl::OutputFile::set_multistate_suffixes(std::vector<std::string>& multiStateSuffixes)
{
    if(nullptr != m_multiStateSuffixes) {
        delete m_multiStateSuffixes;
        m_multiStateSuffixes = nullptr;
    }

    m_multiStateSuffixes = new std::vector<std::string>(multiStateSuffixes);
    return true;
}

void impl::OutputFile::setup_output_params(OutputParams &params) const
{
    bool sort_stk_parts_by_name = false;

    if(m_region->property_exists("sort_stk_parts")) {
        sort_stk_parts_by_name = (m_region->get_property("sort_stk_parts").get_int() != 0);
    }
    params.set_subset_selector(m_subsetSelector.get());
    params.set_shared_selector(m_sharedSelector.get());
    params.set_output_selector(stk::topology::NODE_RANK, m_outputSelector[stk::topology::NODE_RANK].get());
    params.set_output_selector(stk::topology::EDGE_RANK, m_outputSelector[stk::topology::EDGE_RANK].get());
    params.set_output_selector(stk::topology::FACE_RANK, m_outputSelector[stk::topology::FACE_RANK].get());
    params.set_output_selector(stk::topology::ELEM_RANK, m_outputSelector[stk::topology::ELEM_RANK].get());
    params.set_sort_stk_parts_by_name(sort_stk_parts_by_name);
    params.set_use_nodeset_for_block_node_fields(m_useNodesetForBlockNodesFields);
    params.set_use_nodeset_for_sideset_node_fields(m_useNodesetForSidesetNodesFields);
    params.check_field_existence_when_creating_nodesets(m_checkFieldExistenceWhenCreatingNodesets);
    params.set_use_part_id_for_output(m_usePartIdForOutput);
    params.set_has_ghosting(m_hasGhosting);
    params.set_has_adaptivity(m_hasAdaptivity);
    params.set_is_skin_mesh(m_isSkinMesh);
    params.set_additional_attribute_fields(m_additionalAttributeFields);
}

void impl::OutputFile::set_input_region(const Ioss::Region *input_region)
{
    ThrowErrorMsgIf (m_region.get() == input_region,
                     "Attempting to set the input region to the output region");

    m_inputRegion = input_region;
}

void impl::OutputFile::write_output_mesh(const stk::mesh::BulkData& bulk_data,
                                         const std::vector<std::vector<int>> &attributeOrdering)
{
    if ( m_meshDefined == false )
    {
        m_meshDefined = true;

        // If using hdf5 as the underlying file type for exodus/netcdf,
        // it is more picky about overwriting an existing file -- if the
        // file is open, then it will abort; it will only overwrite an existing
        // file if it is not open.  Since overwriting restart files (input/output)
        // is a common usecase, we need to check at this point whether there are
        // any existing input files with the same name as the file we are attempting
        // to create here. However, due to symbolic links and other junk, it is often
        // difficult to determine that the files are the same, so..., If m_input_region
        // refers to a file, just close it since we should be done with it at this time...
        if (m_inputRegion) {
            m_inputRegion->get_database()->closeDatabase();
        }

        // used in stk_adapt/stk_percept
        stk::io::OutputParams params(*m_region, bulk_data);
        setup_output_params(params);

        stk::io::define_output_db(params, attributeOrdering, m_inputRegion);

        if(!m_appendingToMesh)
            stk::io::write_output_db(params);

        //Attempt to avoid putting state change into the interface.  We'll see . . .
        m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }
}

template<typename INT>
void internal_fill_output_entities_for_sideblock(stk::io::OutputParams& params, Ioss::GroupingEntity *ge, stk::mesh::Part *part, stk::mesh::EntityVector& sides)
{
    std::vector<INT> elem_side_ids;
    Ioss::SideBlock *block = dynamic_cast<Ioss::SideBlock*>(ge);
    ThrowRequireMsg(nullptr != block, "Input GroupingEntity is not a sideblock");
    const Ioss::ElementTopology *parent_topology = block->parent_element_topology();
    fill_data_for_side_block(params, *ge, part, parent_topology, elem_side_ids, sides);
}

std::vector<stk::mesh::Entity> impl::OutputFile::get_output_entities(const stk::mesh::BulkData& bulk_data, const std::string &name)
{
    stk::mesh::Part *part = nullptr;
    const stk::mesh::MetaData& meta = bulk_data.mesh_meta_data();

    std::vector<stk::mesh::Entity> entities;

    stk::io::OutputParams params(*m_region, bulk_data);
    setup_output_params(params);

    Ioss::GroupingEntity *ge = m_region->get_entity(name);
    ThrowErrorMsgIf (ge == nullptr,
                     "Could not find grouping entity with name: " + name);

    Ioss::EntityType type = ge->type();
    stk::mesh::EntityRank part_type = stk::mesh::InvalidEntityRank;
    if(type == Ioss::NODEBLOCK) {
        part_type = stk::topology::NODE_RANK;
    } else if(type == Ioss::NODESET) {
        part_type = stk::topology::NODE_RANK;
    } else if(type == Ioss::ELEMENTBLOCK) {
        part_type = stk::topology::ELEMENT_RANK;
    } else if(type == Ioss::SIDESET) {
        part = meta.get_part(name);
        ThrowRequireMsg(nullptr != part, "Could not find a sideset with name: " + name);
        part_type = part->primary_entity_rank();

        Ioss::SideBlock *sb = m_region->get_sideblock(name);
        if(nullptr != sb) {
            ge = sb;
            type = ge->type();
        }
    } else if(type == Ioss::SIDEBLOCK) {
        part = meta.get_part(name);
        ThrowRequireMsg(nullptr != part, "Could not find a sideblock with name: " + name);
        part_type = part->primary_entity_rank();
    }

    if(type == Ioss::SIDEBLOCK) {
        Ioss::Region &io_region = params.io_region();
        bool ints64bit = db_api_int_size(&io_region) == 8;
        if (ints64bit) {
            internal_fill_output_entities_for_sideblock<int64_t>(params, ge, part, entities);
        } else {
            internal_fill_output_entities_for_sideblock<int>(params, ge, part, entities);
        }
    } else {
        get_output_entity_list(ge, part_type, params, entities);
    }

    return entities;
}

void impl::OutputFile::flush_output() const
{
    m_region->get_database()->flush_database();
}

void impl::OutputFile::add_attribute_field(stk::mesh::FieldBase &field, const OutputVariableParams &var)
{
    const std::string &alternate_name = var.name();

    ThrowErrorMsgIf (m_fieldsDefined,
                     "Attempting to add attribute fields after fields have already been written to the database.");
    ThrowErrorMsgIf (alternate_name.empty(),
                     "Attempting to add attribute field " << field.name() << " with no name.");

    stk::io::FieldAndName *existingEntry = nullptr;

    bool fieldAlreadyExists=false;
    for (size_t i=0;i<m_additionalAttributeFields.size();i++) {
        if ( &field          == m_additionalAttributeFields[i].field() &&
                alternate_name == m_additionalAttributeFields[i].db_name() ) {
            existingEntry = &m_additionalAttributeFields[i];
            m_additionalAttributeFields[i].set_db_name(alternate_name);
            fieldAlreadyExists = true;
            break;
        }
    }

    if (!fieldAlreadyExists) {
        stk::io::FieldAndName namedField(&field, alternate_name, field.entity_rank());
        namedField.set_output_params(var);

        if (m_dbPurpose == stk::io::WRITE_RESTART) {
            namedField.set_use_alias(false);
        }

        m_additionalAttributeFields.push_back(namedField);
    } else {
        const std::vector<std::string>& entities = var.get_subset_entities();
        for(const std::string &entity : entities) {
            existingEntry->add_subset_entity(entity);
        }
    }
}

void impl::OutputFile::add_field(stk::mesh::FieldBase &field, const OutputVariableParams &var, stk::mesh::EntityRank var_type)
{
    const std::string &alternate_name = var.name();

    ThrowErrorMsgIf (m_fieldsDefined,
                     "Attempting to add fields after fields have already been written to the database.");
    ThrowErrorMsgIf (alternate_name.empty(),
                     "Attempting to output results field " << field.name() << " with no name.");

    bool fieldAlreadyExists=false;
    for (size_t i=0;i<m_namedFields.size();i++) {
        if ( &field == m_namedFields[i].field() && alternate_name == m_namedFields[i].db_name() ) {
            m_namedFields[i].set_db_name(alternate_name);
            fieldAlreadyExists = true;
            break;
        }
    }

    if (!fieldAlreadyExists) {
        if (m_dbPurpose == stk::io::WRITE_RESTART) {
            int state_count = field.number_of_states();
            STKIORequire(state_count < 7);

            if(state_count > 1) {
                std::vector<std::string>* multiStateSuffixes = state_count > 2 ? m_multiStateSuffixes : nullptr;

                int num_states_to_write = std::max(state_count-1, 1);
                for(int state=0; state < num_states_to_write; state++) {
                    stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
                    stk::mesh::FieldBase *statedField = field.field_state(state_identifier);
                    std::string field_name_with_suffix = stk::io::get_stated_field_name(alternate_name, state_identifier, multiStateSuffixes);
                    stk::io::FieldAndName namedField(statedField, field_name_with_suffix, var_type);
                    namedField.set_use_alias(false);
                    namedField.set_output_params(var);
                    m_namedFields.push_back(namedField);
                    stk::io::set_field_role(*statedField, Ioss::Field::TRANSIENT);
                }
            } else {
                stk::io::FieldAndName namedField(&field, alternate_name, var_type);
                namedField.set_use_alias(false);
                namedField.set_output_params(var);
                m_namedFields.push_back(namedField);
                stk::io::set_field_role(field, Ioss::Field::TRANSIENT);
            }

        } else {
            stk::io::FieldAndName namedField(&field, alternate_name, var_type);
            namedField.set_output_params(var);
            m_namedFields.push_back(namedField);
            stk::io::set_field_role(field, Ioss::Field::TRANSIENT);
        }
    }
}

void impl::OutputFile::add_user_data(const std::vector<std::string>& partNames, const std::string &alternate_name, stk::io::DataLocation loc)
{
    ThrowErrorMsgIf (m_fieldsDefined,
                     "Attempting to add fields after fields have already been written to the database.");
    ThrowErrorMsgIf (alternate_name.empty(),
                     "Attempting to output results field with no name.");

    bool fieldAlreadyExists=false;
    for (size_t i=0;i<m_userData.size();i++) {
        if ( alternate_name == m_userData[i].db_name() ) {
            m_userData[i].set_db_name(alternate_name);
            fieldAlreadyExists = true;
            break;
        }
    }

    if (!fieldAlreadyExists) {
        stk::io::UserDataAndName namedData(partNames, alternate_name, loc);
        m_userData.push_back(namedData);
    }

}

void impl::OutputFile::add_global_ref(const std::string &name, const boost::any *value, stk::util::ParameterType::Type type)
{
    ThrowErrorMsgIf (m_fieldsDefined,
                     "On region named " << m_region->name() <<
                     " Attempting to add global variable after data has already been written to the database.");
    std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, *value);
    internal_add_global(m_region, name, parameter_type.first, parameter_type.second);
    m_globalAnyFields.emplace_back(name, value, type);
}

bool impl::OutputFile::has_global(const std::string &globalVarName) const
{
    return m_region->field_exists(globalVarName);
}

void impl::OutputFile::add_global(const std::string &name, const boost::any &value, stk::util::ParameterType::Type type)
{
    ThrowErrorMsgIf (m_fieldsDefined,
                     "On region named " << m_region->name() <<
                     " Attempting to add global variable after data has already been written to the database.");
    std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, value);
    m_anyGlobalVariablesDefined = true;  // This output file has at least 1 global variable.
    internal_add_global(m_region, name, parameter_type.first, parameter_type.second);
}

void impl::OutputFile::add_global(const std::string &globalVarName, Ioss::Field::BasicType dataType)
{
    ThrowErrorMsgIf (m_fieldsDefined,
                     "On region named " << m_region->name() <<
                     " Attempting to add global variable after data has already been written to the database.");
    m_anyGlobalVariablesDefined = true;  // This output file has at least 1 global variable.
    internal_add_global(m_region, globalVarName, "scalar", dataType);
}

void impl::OutputFile::add_global(const std::string &globalVarName, int component_count, Ioss::Field::BasicType dataType)
{
    ThrowErrorMsgIf (m_fieldsDefined,
                     "On region named " << m_region->name() <<
                     " Attempting to add global variable after data has already been written to the database.");
    m_anyGlobalVariablesDefined = true;  // This output file has at least 1 global variable.
    internal_add_global(m_region, globalVarName, component_count, dataType);
}

void impl::OutputFile::add_global(const std::string &globalVarName, const std::string &storage, Ioss::Field::BasicType dataType)
{
    ThrowErrorMsgIf (m_fieldsDefined,
                     "On region named " << m_region->name() <<
                     " Attempting to add global variable after data has already been written to the database.");
    m_anyGlobalVariablesDefined = true;  // This output file has at least 1 global variable.
    internal_add_global(m_region, globalVarName, storage, dataType);
}

void impl::OutputFile::write_global(const std::string &globalVarName,
                                    const boost::any &value, stk::util::ParameterType::Type type)
{
    internal_write_parameter(m_region, globalVarName, value, type);
}

void impl::OutputFile::write_global(const std::string &globalVarName, std::vector<double>& globalVarData)
{
    internal_write_global(m_region, globalVarName, globalVarData);
}

void impl::OutputFile::write_global(const std::string &globalVarName, std::vector<int>& globalVarData)
{
    internal_write_global(m_region, globalVarName, globalVarData);
}

void impl::OutputFile::write_global(const std::string &globalVarName, int globalVarData)
{
    internal_write_global(m_region, globalVarName, globalVarData);
}

void impl::OutputFile::write_global(const std::string &globalVarName, double globalVarData)
{
    internal_write_global(m_region, globalVarName, globalVarData);
}

void impl::OutputFile::setup_output_file(const std::string &filename, stk::ParallelMachine communicator,
                                         Ioss::PropertyManager &property_manager, char const* type,
                                         bool openFileImmediately)
{
    ThrowErrorMsgIf (filename.empty(),
                     "No filename was specified for the output file creation.");
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(type, filename,
                                                    Ioss::WRITE_RESULTS,
                                                    communicator,
                                                    property_manager);

    if (dbo == nullptr || (openFileImmediately && !dbo->ok())) {
        std::cerr << "ERROR: Could not open output database '" << filename << "' of type 'exodus'\n";
        std::exit(EXIT_FAILURE);
    }

    // There is no "label" for the output region; just use the filename for now so
    // Ioss messages will specify a unique or identifiable instance.
    m_region = Teuchos::rcp(new Ioss::Region(dbo, filename));

    if(property_manager.exists("APPEND_OUTPUT") && property_manager.get("APPEND_OUTPUT").get_int() == Ioss::DB_APPEND)
        m_appendingToMesh = true;
}

void impl::OutputFile::begin_output_step(double time, const stk::mesh::BulkData& bulk_data,
                                         const std::vector<std::vector<int>> &attributeOrdering)
{
    if (!m_fieldsDefined) {
        define_output_fields(bulk_data, attributeOrdering);
    }

    //Attempt to avoid putting state change into the interface.  We'll see . . .
    Ioss::State currentState = m_region->get_state();
    if(currentState == Ioss::STATE_DEFINE_TRANSIENT) {
        m_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }

    if(m_region->get_state() != Ioss::STATE_TRANSIENT)
    {
        m_region->begin_mode(Ioss::STATE_TRANSIENT);
    }
    m_currentOutputStep = m_region->add_state(time);
    m_region->begin_state(m_currentOutputStep);
}

// ========================================================================
// Iterate over all fields defined in the stk mesh data structure.
// If the field has the io_attribute set, then define that field
// on the corresponding io entity on the output mesh database.
// The database field will have the same name as the stk field.
//
// To export the data to the database, call
// process_output_request().
void impl::OutputFile::define_output_fields(const stk::mesh::BulkData& bulk_data,
                                            const std::vector<std::vector<int>> &attributeOrdering)
{
    if(m_fieldsDefined) {
        return;
    }

    write_output_mesh(bulk_data, attributeOrdering);

    Ioss::Region *region = m_region.get();

    // Sort the fields in m_named_fields based on the field ordinal.
    // This is needed so all processors will process the fields in the same order
    // Not guaranteed that the application has added the fields in order, so make
    // the guarantee here...
    std::sort(m_namedFields.begin(), m_namedFields.end(), fieldOrdinalSort);

    //          std::cerr << "In define_output_fields" << std::endl;

    const stk::mesh::MetaData &meta_data = bulk_data.mesh_meta_data();
    // Special processing for nodeblock (all nodes in model)...
    stk::io::ioss_add_fields(meta_data.universal_part(), stk::topology::NODE_RANK,
                             region->get_node_blocks()[0], m_namedFields);

    const stk::mesh::PartVector &all_parts = meta_data.get_parts();
    for(auto part : all_parts) {
        bool isIoPart = stk::io::is_part_io_part(*part);

        // Check whether this part should be output to database.
        if(isIoPart) {
            stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
            // Get Ioss::GroupingEntity corresponding to this part...
            std::string partName = getPartName(*part);
            Ioss::GroupingEntity *entity = region->get_entity(partName);

            if(entity!= nullptr)
            {
                if (entity->type() == Ioss::SIDESET)
                {
                    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
                    int block_count = sset->block_count();

                    for(int i = 0; i < block_count; i++)
                    {
                        Ioss::SideBlock *fb = sset->get_block(i);
                        stk::io::ioss_add_fields(*part, rank, fb, m_namedFields);
                    }
                }
                else {
                    stk::io::ioss_add_fields(*part, rank, entity, m_namedFields);
                }
                // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
                // (should probably do edges and faces also...)
                // Get Ioss::GroupingEntity corresponding to the nodes on this part...
                if(rank != stk::topology::NODE_RANK) {
                    Ioss::GroupingEntity *node_entity = nullptr;

                    if (rank == meta_data.side_rank() || rank == stk::topology::ELEM_RANK) {
                        bool use_nodeset = (rank == meta_data.side_rank()) ? m_useNodesetForSidesetNodesFields
                                : m_useNodesetForBlockNodesFields;

                        if (use_nodeset) {
                            std::string nodes_name = partName + s_entity_nodes_suffix;
                            node_entity = region->get_entity(nodes_name);
                        }
                    }

                    if(node_entity != nullptr) {
                        stk::io::ioss_add_fields(*part, stk::topology::NODE_RANK, node_entity, m_namedFields);
                    }
                }
            }
        }
    }
    m_fieldsDefined = true;
}

int impl::OutputFile::process_output_request(double time, const stk::mesh::BulkData& bulk_data,
                                             const std::vector<std::vector<int>> &attributeOrdering)
{
    ThrowErrorMsgIf(m_anyGlobalVariablesDefined,
                    "The output database " << m_region->name() << " has defined global variables, "
                    "but is calling the process_output_request() function which does not output global "
                    "variables.  Call begin_output_step() instead.");

    begin_output_step(time, bulk_data, attributeOrdering);
    write_defined_output_fields(bulk_data);
    write_defined_global_any_fields(m_region, m_globalAnyFields);
    end_output_step();

    return m_currentOutputStep;
}

int impl::OutputFile::write_defined_output_fields(const stk::mesh::BulkData& bulk_data, const stk::mesh::FieldState *state)
{
    Ioss::Region *region = m_region.get();
    ThrowErrorMsgIf (region==nullptr, "INTERNAL ERROR: Mesh Output Region pointer is NULL in write_defined_output_fields.");

    OutputParams params(*region, bulk_data);
    setup_output_params(params);

    const stk::mesh::MetaData& meta_data = bulk_data.mesh_meta_data();
    // Special processing for nodeblock (all nodes in model)...
    put_field_data(params, meta_data.universal_part(), stk::topology::NODE_RANK,
                   region->get_node_blocks()[0], m_namedFields, state);

    // Now handle all non-nodeblock parts...
    const stk::mesh::PartVector &all_parts = meta_data.get_parts();
    for( auto part : all_parts ) {
        // Check whether this part should be output to database.
        bool isIoPart = stk::io::is_part_io_part(*part);

        if (isIoPart) {
            stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
            if( rank == stk::topology::INVALID_RANK ) continue;
            // Get Ioss::GroupingEntity corresponding to this part...
            std::string partName = getPartName(*part);
            auto entity = region->get_entity(partName);

            // If the sideset has only a single sideblock and it
            // shares the same name as the parent sideset, then the
            // entity that we want to output on at this point is the
            // sideblock and not the sideset.  If there are multiple
            // sideblocks in the sideset, then they will be output separately...
            if ( entity != nullptr && entity->type() == Ioss::SIDESET ) {
                auto sset = dynamic_cast<Ioss::SideSet*>(entity);
                size_t block_count = sset->block_count();
                if ( block_count == 1 ) {
                    auto ssblock = sset->get_side_block(partName);
                    if ( ssblock ) {
                        entity = ssblock; // NOTE: 'entity' is reset at this point.
                    }
                }
            }

            if ( entity != nullptr && entity->type() != Ioss::SIDESET ) {
                put_field_data(params, *part, rank, entity, m_namedFields, state);
            }

            // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
            // (should probably do edges and faces also...)
            // Get Ioss::GroupingEntity corresponding to the nodes on this part...

            if(rank != stk::topology::NODE_RANK && entity != nullptr) {
                Ioss::GroupingEntity *node_entity = nullptr;

                if (rank == meta_data.side_rank() || rank == stk::topology::ELEM_RANK) {
                    bool use_nodeset = (rank == meta_data.side_rank()) ? m_useNodesetForSidesetNodesFields :
                            m_useNodesetForBlockNodesFields;

                    if (use_nodeset) {
                        std::string nodes_name = partName + s_entity_nodes_suffix;
                        node_entity = region->get_entity(nodes_name);
                    }
                }

                if(node_entity != nullptr) {
                    put_field_data(params, *part, stk::topology::NODE_RANK, node_entity, m_namedFields, state);
                }
            }
        }
    }
    return m_currentOutputStep;
}

int impl::OutputFile::write_defined_output_fields_for_selected_subset(const stk::mesh::BulkData& bulk_data,
                                                                      std::vector<stk::mesh::Part*>& selectOutputElementParts,
                                                                      const stk::mesh::FieldState *state)
{
    Ioss::Region *region = m_region.get();
    ThrowErrorMsgIf (region==nullptr, "INTERNAL ERROR: Mesh Output Region pointer is NULL in write_defined_output_fields.");

    OutputParams params(*region, bulk_data);
    setup_output_params(params);

    const stk::mesh::MetaData& meta_data = bulk_data.mesh_meta_data();
    // Special processing for nodeblock (all nodes in model)...
    put_field_data(params, meta_data.universal_part(), stk::topology::NODE_RANK,
                   region->get_node_blocks()[0], m_namedFields, state);

    for( auto part : selectOutputElementParts )
    {
        stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
        auto entity = region->get_entity(part->name());

        if ( entity != nullptr && entity->type() != Ioss::SIDESET ) {
            put_field_data(params, *part, rank, entity, m_namedFields, state);
        }
    }
    return m_currentOutputStep;
}

void impl::OutputFile::end_output_step()
{
    m_region->end_state(m_currentOutputStep);
}

void impl::OutputFile::set_subset_selector(Teuchos::RCP<stk::mesh::Selector> my_selector)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: On region named " << m_region->name() <<
                    " the subset_selector cannot be changed after the mesh has already been written.");
    m_subsetSelector = my_selector;
}

void impl::OutputFile::set_shared_selector(Teuchos::RCP<stk::mesh::Selector> my_selector)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: On region named " << m_region->name() <<
                    " the shared_selector cannot be changed after the mesh has already been written.");
    m_sharedSelector = my_selector;
}

void impl::OutputFile::set_output_selector(stk::topology::rank_t rank, Teuchos::RCP<stk::mesh::Selector> my_selector)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: On region named " << m_region->name() <<
                    " the output_selector cannot be changed after the mesh has already been written.");

    ThrowErrorMsgIf(!(rank >= stk::topology::NODE_RANK && rank <= stk::topology::ELEM_RANK),
                    "ERROR: On region named " << m_region->name() <<
                    " the output_selector must be NODE, EDGE, FACE or ELEM.");

    m_outputSelector[rank] = my_selector;
}

bool impl::OutputFile::use_nodeset_for_block_nodes_fields() const
{
    return m_useNodesetForBlockNodesFields;
}

void impl::OutputFile::use_nodeset_for_block_nodes_fields(bool true_false)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: The use_nodeset_for_block_nodes_fields setting cannot be changed after "
                    "the mesh has already been written.");
    m_useNodesetForBlockNodesFields = true_false;
}

bool impl::OutputFile::use_nodeset_for_sideset_nodes_fields() const
{
    return m_useNodesetForSidesetNodesFields;
}

void impl::OutputFile::use_nodeset_for_sideset_nodes_fields(bool true_false)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: The use_nodeset_for_sideset_nodes_fields setting cannot be changed after "
                    "the mesh has already been written.");
    m_useNodesetForSidesetNodesFields = true_false;
}

bool impl::OutputFile::check_field_existence_when_creating_nodesets() const
{
    return m_checkFieldExistenceWhenCreatingNodesets;
}

void impl::OutputFile::check_field_existence_when_creating_nodesets(bool true_false)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: The check_field_existence_when_creating_nodesets setting cannot be changed after "
                    "the mesh has already been written.");
    m_checkFieldExistenceWhenCreatingNodesets = true_false;
}

bool impl::OutputFile::use_part_id_for_output() const
{
    return m_usePartIdForOutput;
}

void impl::OutputFile::use_part_id_for_output(bool true_false)
{
    ThrowErrorMsgIf(m_meshDefined,
                    "ERROR: The use_part_id_for_output setting cannot be changed after "
                    "the mesh has already been written.");
    m_usePartIdForOutput = true_false;
}

bool impl::OutputFile::has_ghosting() const
{
    return m_hasGhosting;
}

void impl::OutputFile::has_ghosting(bool hasGhosting)
{
    m_hasGhosting = hasGhosting;
}

bool impl::OutputFile::has_adaptivity() const
{
    return m_hasAdaptivity;
}

void impl::OutputFile::has_adaptivity(bool hasAdaptivity)
{
    m_hasAdaptivity = hasAdaptivity;
}

bool impl::OutputFile::is_skin_mesh() const
{
    return m_isSkinMesh;
}

void impl::OutputFile::is_skin_mesh(bool skinMesh)
{
    m_isSkinMesh = skinMesh;
}

} // namespace io
} // namespace stk
