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
#include <stk_io/InputFile.hpp>                      // for InputFile
#include <stk_io/OutputFile.hpp>                     // for InputFile
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

template <typename T>
bool is_index_valid(const std::vector<T> &file_vector, size_t input_file_index)
{
    bool invalid = file_vector.empty() ||
            input_file_index >= file_vector.size() ||
            Teuchos::is_null(file_vector[input_file_index]);
    return !invalid;
}


bool is_skipped_attribute_field(const std::string &name, size_t numAttrFields)
{
    return name == "attribute" && numAttrFields > 1;
}

std::vector<std::string> get_ordered_attribute_field_names(const Ioss::GroupingEntity &iossElementBlock)
{
    std::vector<std::string> names;
    iossElementBlock.field_describe(Ioss::Field::ATTRIBUTE, &names);
    std::vector<std::pair<int, std::string>> indexAndName;
    for(size_t i = 0; i < names.size(); i++) {
        int attributeIndex = iossElementBlock.get_field(names[i]).get_index();
        if(!is_skipped_attribute_field(names[i], names.size()))
            indexAndName.emplace_back(attributeIndex, names[i]);
    }
    std::sort(indexAndName.begin(), indexAndName.end());

    names.resize(indexAndName.size());
    for(size_t i=0;i<names.size();++i) {
        names[i] = indexAndName[i].second;
    }
    return names;
}

stk::mesh::FieldVector get_fields_by_name(const stk::mesh::MetaData &meta, const std::vector<std::string> &names)
{
    stk::mesh::FieldVector attrFields(names.size());
    for(size_t i=0;i<attrFields.size();++i)
    {
        attrFields[i] = meta.get_field(stk::topology::ELEM_RANK, names[i]);
        ThrowRequireMsg(attrFields[i] != nullptr, "Can't find field named " << names[i]);
    }
    return attrFields;
}

StkMeshIoBroker::StkMeshIoBroker()
: m_communicator(MPI_COMM_NULL),
  m_connectivityMap(nullptr),
  m_activeMeshIndex(0),
  m_sidesetFaceCreationBehavior(STK_IO_SIDE_CREATION_USING_GRAPH_TEST),
  m_autoLoadAttributes(true)
{
    Ioss::Init::Initializer::initialize_ioss();
}

StkMeshIoBroker::StkMeshIoBroker(stk::ParallelMachine comm, const stk::mesh::ConnectivityMap * connectivity_map)
: m_communicator(comm),
  m_connectivityMap(connectivity_map),
  m_activeMeshIndex(0),
  m_sidesetFaceCreationBehavior(STK_IO_SIDE_CREATION_USING_GRAPH_TEST),
  m_autoLoadAttributes(true)
{
    Ioss::Init::Initializer::initialize_ioss();
}

StkMeshIoBroker::~StkMeshIoBroker()
{
}

void StkMeshIoBroker::property_add(const Ioss::Property &property)
{
    m_propertyManager.add(property);
    //In case there are already input/output files, put the property on them too.
    if (get_input_io_region().get() != nullptr) {
        get_input_io_region()->property_add(property);
    }
    for (size_t i=0; i<m_outputFiles.size(); ++i) {
        get_output_io_region(i)->property_add(property);
    }
}

Ioss::Property StkMeshIoBroker::property_get(const std::string &property_name) const
{
    return m_propertyManager.get(property_name);
}

bool StkMeshIoBroker::property_exists(const std::string &property_name) const
{
    return m_propertyManager.exists(property_name);
}

void StkMeshIoBroker::copy_property(const StkMeshIoBroker& src_broker, const std::string &property_name)
{
    if(src_broker.property_exists(property_name)) {
        Ioss::Property property = src_broker.property_get(property_name);
        property_add(property);
    }
}

void StkMeshIoBroker::remove_property_if_exists(const std::string &property_name)
{
    m_propertyManager.erase(property_name);
}

stk::mesh::FieldBase const& StkMeshIoBroker::get_coordinate_field()
{
    stk::mesh::FieldBase const* coord_field = meta_data().coordinate_field();
    STKIORequire( coord_field != nullptr);
    return * coord_field;
}

size_t StkMeshIoBroker::add_mesh_database(Teuchos::RCP<Ioss::Region> ioss_input_region)
{
    auto input_file = Teuchos::rcp(new InputFile(ioss_input_region));
    m_inputFiles.push_back(input_file);

    size_t index_of_input_file = m_inputFiles.size()-1;
    return index_of_input_file;
}

void StkMeshIoBroker::create_sideset_observer()
{
    ThrowRequireMsg( !Teuchos::is_null(m_bulkData), "Bulk data not initialized");
    if (!bulk_data().has_observer_type<SidesetUpdater>()) {
        stk::mesh::Selector activeSelector = get_active_selector();
        if (activeSelector == stk::mesh::Selector()) {
            activeSelector = !activeSelector;
        }
        bulk_data().register_observer(std::make_shared<SidesetUpdater>(bulk_data(), activeSelector));
    }
}

void StkMeshIoBroker::set_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data )
{
    ThrowErrorMsgIf( !Teuchos::is_null(m_bulkData),
                     "Bulk data already initialized" );
    m_bulkData = arg_bulk_data;

    if (Teuchos::is_null(m_metaData)) {
        m_metaData = Teuchos::rcpFromRef(bulk_data().mesh_meta_data());
    }

#ifdef STK_BUILT_IN_SIERRA
    m_communicator = m_bulkData->parallel();
#endif
    create_sideset_observer();
}

void StkMeshIoBroker::replace_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data )
{
    ThrowErrorMsgIf( Teuchos::is_null(m_bulkData),
                     "There is  no bulk data to replace." );
    ThrowErrorMsgIf( Teuchos::is_null(m_metaData),
                     "Meta data must be non-null when calling StkMeshIoBroker::replace_bulk_data." );

    stk::mesh::MetaData &new_meta_data = arg_bulk_data->mesh_meta_data();
    ThrowErrorMsgIf( &(*m_metaData) != &new_meta_data,
                     "Meta data for both new and old bulk data must be the same." );

    m_bulkData = arg_bulk_data;
    create_sideset_observer();
}

size_t StkMeshIoBroker::add_mesh_database(std::string filename, DatabasePurpose purpose)
{
    std::string type = "exodus";

    // See if filename contains a ":" at the beginning of the filename
    // and if the text preceding that filename specifies a valid
    // database type.  If so, set that as the file type and extract
    // the portion following the colon as the filename.
    // If no colon in name, use default type.

    size_t colon = filename.find(':');
    if (colon != std::string::npos && colon > 0) {
        type = filename.substr(0, colon);
        filename = filename.substr(colon+1);
    }
    return add_mesh_database(filename, type, purpose);
}


size_t StkMeshIoBroker::add_mesh_database(const std::string &filename,
                                          const std::string &type,
                                          DatabasePurpose purpose)
{
    auto input_file = Teuchos::rcp(new InputFile(filename, m_communicator, type, purpose, m_propertyManager));
    m_inputFiles.push_back(input_file);

    size_t index_of_input_file = m_inputFiles.size()-1;
    return index_of_input_file;
}

void StkMeshIoBroker::copy_property_manager(const Ioss::PropertyManager &properties)
{
    Ioss::NameList props;
    int num_prop = properties.describe(&props);
    for(int i = 0; i < num_prop; i++) {
        m_propertyManager.add(properties.get(props[i]));
    }
}

size_t StkMeshIoBroker::add_mesh_database(const std::string &filename,
                                          const std::string &type,
                                          DatabasePurpose purpose,
                                          Ioss::PropertyManager &properties)
{
    copy_property_manager(properties);
    return add_mesh_database(filename, type,purpose);
}

Teuchos::RCP<Ioss::Region> StkMeshIoBroker::get_input_io_region()
{
    if (is_index_valid(m_inputFiles, m_activeMeshIndex)) {
        return m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    } else {
        return Teuchos::RCP<Ioss::Region>();
    }
}

InputFile &StkMeshIoBroker::get_mesh_database(size_t input_file_index)
{
    validate_input_file_index(input_file_index);
    return *m_inputFiles[input_file_index];
}

void StkMeshIoBroker::remove_mesh_database(size_t input_file_index)
{
    validate_input_file_index(input_file_index);
    // It would be nice to be able to just delete the database, but
    // its io_region may be being used by one or more output files
    // in the Ioss::Region::synchronize... function.  Therefore, we
    // need to keep it around (Could also use some shared pointers,
    // but not allowed to use C++11 yet).  But, we want this database to be
    // inaccessible and close the file associated with the database.  Therefore,
    // we add an empty InputFile to the end of 'm_input_files' and then swap it with
    // this one--That way the 'input_file_index' points to an invalid InputFile class as
    // it should.
    m_inputFiles[input_file_index]->get_input_io_region()->get_database()->closeDatabase();
    m_inputFiles.push_back(m_inputFiles[input_file_index]);
    m_inputFiles[input_file_index] = Teuchos::RCP<InputFile>();
    assert(Teuchos::is_null(m_inputFiles[input_file_index]));
}

size_t StkMeshIoBroker::set_active_mesh(size_t input_file_index)
{
    validate_input_file_index(input_file_index);
    size_t old = m_activeMeshIndex;
    m_activeMeshIndex = input_file_index;
    m_inputFiles[m_activeMeshIndex]->create_ioss_region();
    return old;
}

void StkMeshIoBroker::create_ioss_region()
{
    validate_input_file_index(m_activeMeshIndex);
    m_inputFiles[m_activeMeshIndex]->create_ioss_region();
}

void StkMeshIoBroker::set_rank_name_vector(const std::vector<std::string> &rank_names)
{
    ThrowErrorMsgIf(!Teuchos::is_null(m_metaData),
                    "There meta data associated with this StkMeshIoBroker has already been created. "
                    "It is not permissible to set the rank_name_vector() at this time.");

    m_rankNames.clear();
    m_rankNames.insert( m_rankNames.end(),rank_names.begin(), rank_names.end());
}

bool StkMeshIoBroker::is_output_index_valid(size_t outputIndex)
{
    return is_index_valid(m_outputFiles, outputIndex);
}

bool StkMeshIoBroker::is_input_index_valid(size_t inputIndex)
{
    return is_index_valid(m_inputFiles, inputIndex);
}

std::string StkMeshIoBroker::get_output_filename(size_t outputIndex)
{
    if (!is_output_index_valid(outputIndex)) {
        return "";
    }
    Teuchos::RCP<Ioss::Region> outputRegion = m_outputFiles[outputIndex]->get_output_io_region();
    Ioss::DatabaseIO* outputDatabase = outputRegion->get_database();
    return outputDatabase->get_filename();
}

void StkMeshIoBroker::store_attribute_field_ordering()
{
    const stk::mesh::PartVector& parts = meta_data().get_parts();
    attributeFieldOrderingByPartOrdinal.clear();
    attributeFieldOrderingByPartOrdinal.resize(parts.size());
    for(const stk::mesh::Part* part : parts)  {
        const Ioss::GroupingEntity* iossGroupingEntity = part->attribute<Ioss::GroupingEntity>();
        if(iossGroupingEntity != nullptr) {
            std::vector<std::string> names = get_ordered_attribute_field_names(*iossGroupingEntity);
            const stk::mesh::FieldVector attrFields = get_fields_by_name(*m_metaData, names);
            int partOrd = part->mesh_meta_data_ordinal();
            attributeFieldOrderingByPartOrdinal[partOrd].resize(attrFields.size());
            for(size_t i = 0; i < attrFields.size(); i++) {
                attributeFieldOrderingByPartOrdinal[partOrd][i] = attrFields[i]->mesh_meta_data_ordinal();
            }
        }
    }
}

void StkMeshIoBroker::create_surface_to_block_mapping()
{
    IossBlockMembership blockMemberships = get_block_memberships(*this);
    for(IossBlockMembership::iterator iter = blockMemberships.begin(); iter != blockMemberships.end(); iter++) {
        stk::mesh::Part* sidesetPart = meta_data().get_part(iter->first);
        if(sidesetPart != nullptr && sidesetPart->primary_entity_rank() == meta_data().side_rank()) {
            std::vector<const stk::mesh::Part*> blocks;
            fill_block_parts_given_names(iter->second, meta_data(), blocks);
            meta_data().set_surface_to_block_mapping(sidesetPart, blocks);
        }
    }
}

void StkMeshIoBroker::create_input_mesh()
{
    validate_input_file_index(m_activeMeshIndex);
    if (Teuchos::is_null(m_inputFiles[m_activeMeshIndex]->get_input_io_region())) {
        m_inputFiles[m_activeMeshIndex]->create_ioss_region();
    }

    Ioss::Region *region = m_inputFiles[m_activeMeshIndex]->get_input_io_region().get();
    ThrowErrorMsgIf (region==nullptr,
                     "INTERNAL ERROR: Mesh Input Region pointer is NULL in create_input_mesh.");

    // See if meta data is null, if so, create a new one...
    if (Teuchos::is_null(m_metaData)) {
        m_metaData = Teuchos::rcp(new stk::mesh::MetaData());
    }

    size_t spatial_dimension = region->get_property("spatial_dimension").get_int();
    if (m_rankNames.empty()) {
        initialize_spatial_dimension(meta_data(), spatial_dimension, stk::mesh::entity_rank_names());
    } else {
        initialize_spatial_dimension(meta_data(), spatial_dimension, m_rankNames);
    }

    process_nodeblocks(*region,    meta_data());
    process_elementblocks(*region, meta_data());
    process_sidesets(*region,      meta_data());
    process_nodesets(*region,      meta_data());

    create_surface_to_block_mapping();
    store_attribute_field_ordering();
}


size_t StkMeshIoBroker::create_output_mesh(const std::string &filename, DatabasePurpose db_type,
                                           char const* type, bool openFileImmediately)
{
    return create_output_mesh(filename, db_type, m_propertyManager, type, openFileImmediately);
}

size_t StkMeshIoBroker::create_output_mesh(const std::string &filename, DatabasePurpose db_type,
                                           double time,
                                           char const* type, bool openFileImmediately)
{
    return create_output_mesh(filename, db_type, m_propertyManager, time, type, openFileImmediately);
}

size_t StkMeshIoBroker::create_output_mesh(const std::string &filename, DatabasePurpose db_type,
                                           Ioss::PropertyManager &properties,
                                           double time,
                                           char const* type, bool openFileImmediately)
{
    if(db_type == stk::io::APPEND_RESULTS) {
        properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND));
        properties.add(Ioss::Property("APPEND_OUTPUT_AFTER_TIME", time));
    }

    return create_output_mesh(filename, db_type, properties, type, openFileImmediately);
}

size_t StkMeshIoBroker::create_output_mesh(const std::string &filename, DatabasePurpose db_type,
                                           Ioss::PropertyManager &properties,
                                           char const* type, bool openFileImmediately)
{
    if(db_type == stk::io::APPEND_RESULTS) {
        properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND));
        db_type = stk::io::WRITE_RESULTS;
    }

    std::string out_filename = filename;
    stk::util::filename_substitution(out_filename);
    Ioss::Region *input_region = nullptr;
    if (is_index_valid(m_inputFiles, m_activeMeshIndex)) {
        input_region = get_input_io_region().get();
    }

    // Determine whether 64-bit integers are required for the output mesh...
    if (!properties.exists("INTEGER_SIZE_DB")){
        bool requires_64bit = check_integer_size_requirements() == 8;
        if (requires_64bit) {
            properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
        }
    }

    auto output_file = Teuchos::rcp(new impl::OutputFile(out_filename, m_communicator, db_type,
                                                         properties, input_region, type, openFileImmediately));
    m_outputFiles.push_back(output_file);

    size_t index_of_output_file = m_outputFiles.size()-1;
    return index_of_output_file;
}


void StkMeshIoBroker::write_output_mesh(size_t output_file_index)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->write_output_mesh(*m_bulkData, attributeFieldOrderingByPartOrdinal);
}

void StkMeshIoBroker::flush_output() const
{
    for (const auto& out_file : m_outputFiles) {
        out_file->flush_output();
    }
    for (const auto& hb : m_heartbeat) {
        hb->flush_output();
    }
}

int StkMeshIoBroker::write_defined_output_fields(size_t output_file_index, const stk::mesh::FieldState *state)
{
    validate_output_file_index(output_file_index);
    int current_output_step = m_outputFiles[output_file_index]->write_defined_output_fields(*m_bulkData, state);
    return current_output_step;
}

int StkMeshIoBroker::write_defined_output_fields_for_selected_subset(size_t output_file_index,
                                                                     std::vector<stk::mesh::Part*>& selectOutputElementParts,
                                                                     const stk::mesh::FieldState *state)
{
    validate_output_file_index(output_file_index);
    int current_output_step = m_outputFiles[output_file_index]->write_defined_output_fields_for_selected_subset(*m_bulkData,
                                                                                                                 selectOutputElementParts,
                                                                                                                 state);
    return current_output_step;
}

int StkMeshIoBroker::process_output_request(size_t output_file_index, double time)
{
    validate_output_file_index(output_file_index);
    int current_output_step = m_outputFiles[output_file_index]->process_output_request(time, *m_bulkData, attributeFieldOrderingByPartOrdinal);
    return current_output_step;
}

void StkMeshIoBroker::begin_output_step(size_t output_file_index, double time)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->begin_output_step(time, *m_bulkData, attributeFieldOrderingByPartOrdinal);
}

void StkMeshIoBroker::end_output_step(size_t output_file_index)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->end_output_step();
}

template<typename INT>
void populate_elements_and_nodes(Ioss::Region &region,
                                 stk::mesh::BulkData& bulkData,
                                 stk::ParallelMachine communicator,
                                 const bool processAllInputNodes)
{
    if(processAllInputNodes) {
#ifdef STK_BUILT_IN_SIERRA
        process_nodeblocks<INT>(region,    bulkData);
#else
        process_nodeblocks<INT>(region,    bulkData, communicator);
#endif
    }

    process_elementblocks<INT>(region, bulkData);

    if(!processAllInputNodes) {
#ifdef STK_BUILT_IN_SIERRA
        process_node_sharing<INT>(region,    bulkData);
#else
        process_node_sharing<INT>(region,    bulkData, communicator);
#endif
    }

    process_nodesets<INT>(region, bulkData);
    process_hidden_nodesets<INT>(region, bulkData);
}

bool StkMeshIoBroker::populate_mesh_elements_and_nodes(bool delay_field_data_allocation)
{
    validate_input_file_index(m_activeMeshIndex);

    create_bulk_data();

    if (delay_field_data_allocation) {
        bulk_data().deactivate_field_updating();
    }

    toggle_sideset_updaters(bulk_data(), false);

    bool i_started_modification_cycle = bulk_data().modification_begin("Mesh Read");

    Ioss::Region *region = m_inputFiles[m_activeMeshIndex]->get_input_io_region().get();
    bool ints64bit = db_api_int_size(region) == 8;
    bool processAllInputNodes = true;
    if(region->property_exists(stk::io::s_process_all_input_nodes)) {
        processAllInputNodes = region->get_property(stk::io::s_process_all_input_nodes).get_int();
    }

    if (ints64bit) {
        populate_elements_and_nodes<int64_t>(*region, bulk_data(), m_communicator, processAllInputNodes);
    } else {
        populate_elements_and_nodes<int>(*region, bulk_data(), m_communicator, processAllInputNodes);
    }

    stk_mesh_resolve_node_sharing();

    toggle_sideset_updaters(bulk_data(), true);

    return i_started_modification_cycle;
}

void StkMeshIoBroker::populate_mesh_sidesets(bool i_started_modification_cycle)
{
    validate_input_file_index(m_activeMeshIndex);

    toggle_sideset_updaters(bulk_data(), false);

    Ioss::Region *region = m_inputFiles[m_activeMeshIndex]->get_input_io_region().get();
    stk::mesh::EntityIdProcMap elemIdMovedToProc;

    if(m_sidesetFaceCreationBehavior!=STK_IO_SIDE_CREATION_USING_GRAPH_TEST) {
        process_sidesets(*region, bulk_data(), elemIdMovedToProc, m_sidesetFaceCreationBehavior);
        bool saveOption = bulk_data().use_entity_ids_for_resolving_sharing();
        bulk_data().set_use_entity_ids_for_resolving_sharing(true);
        stk_mesh_modification_end_after_node_sharing_resolution();
        bulk_data().set_use_entity_ids_for_resolving_sharing(saveOption);
    } else {
        bulk_data().initialize_face_adjacent_element_graph();
        process_sidesets(*region, bulk_data(), elemIdMovedToProc, m_sidesetFaceCreationBehavior);
        stk_mesh_modification_end_after_node_sharing_resolution();
    }

    // Not sure if this is needed anymore. Don't think it'll be called with a nested modification cycle
    if(!i_started_modification_cycle)
        bulk_data().modification_begin();

    toggle_sideset_updaters(bulk_data(), true);
}

void StkMeshIoBroker::populate_mesh(bool delay_field_data_allocation)
{
    bool i_started_modification_cycle = populate_mesh_elements_and_nodes(delay_field_data_allocation);
    populate_mesh_sidesets(i_started_modification_cycle);
}

template<typename INT>
void process_field_data(Ioss::Region &region, stk::mesh::BulkData& bulkData, const bool autoLoadAttributes)
{
    process_node_coords_and_attributes<INT>(region, bulkData);
    process_elem_attributes_and_implicit_ids<INT>(region, bulkData, autoLoadAttributes);
    process_nodesets_df<INT>(region, bulkData);
    process_sidesets_df(region, bulkData);
}

void StkMeshIoBroker::populate_field_data()
{
    validate_input_file_index(m_activeMeshIndex);

    //if field-data has already been allocated, then the allocate_field_data() method
    //is a harmless no-op.
    bulk_data().allocate_field_data();

    Ioss::Region *region = m_inputFiles[m_activeMeshIndex]->get_input_io_region().get();
    ThrowErrorMsgIf (region==nullptr,
                     "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_field_data.");

    bool ints64bit = db_api_int_size(region) == 8;
    if (ints64bit) {
        process_field_data<int64_t>(*region, bulk_data(), m_autoLoadAttributes);
    } else {
        process_field_data<int>(*region, bulk_data(), m_autoLoadAttributes);
    }
}

void StkMeshIoBroker::create_bulk_data()
{
    if (!meta_data().is_commit()) {
        meta_data().commit();
    }

    validate_input_file_index(m_activeMeshIndex);
    ThrowErrorMsgIf (Teuchos::is_null(m_inputFiles[m_activeMeshIndex]->get_input_io_region()),
                     "There is no Input mesh region associated with this Mesh Data.");

    Ioss::Region *region = m_inputFiles[m_activeMeshIndex]->get_input_io_region().get();
    ThrowErrorMsgIf (region==nullptr,
                     "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_mesh.");

    // Check if bulk_data is null; if so, create a new one...
    if (Teuchos::is_null(m_bulkData)) {
        set_bulk_data(Teuchos::rcp( new stk::mesh::BulkData(   meta_data()
                                                               , region->get_database()->util().communicator()
                                                               , stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
                                                               , false
#endif
                                                               , m_connectivityMap
        )));
    }
}

// ========================================================================
void StkMeshIoBroker::populate_bulk_data()
{
    validate_input_file_index(m_activeMeshIndex);

    create_bulk_data();

    // to preserve behavior for callers of this method, don't do the
    // delay-field-data-allocation optimization.
    // If want the optimization, call the population_mesh/populate_field_data methods separately.

    bool delay_field_data_allocation = false;
    populate_mesh(delay_field_data_allocation);

    populate_field_data();

    if(m_bulkData->is_automatic_aura_on()) {
        std::vector< const stk::mesh::FieldBase *> fields(m_metaData->get_fields().begin(), m_metaData->get_fields().end());
        stk::mesh::communicate_field_data(m_bulkData->aura_ghosting(), fields);
    }
}

void StkMeshIoBroker::add_input_field(const stk::io::MeshField &mesh_field)
{
    add_input_field(m_activeMeshIndex, mesh_field);
}

void StkMeshIoBroker::add_input_field(size_t mesh_index, const stk::io::MeshField &mesh_field)
{
    validate_input_file_index(mesh_index);
    m_inputFiles[mesh_index]->add_input_field(mesh_field);
}

void StkMeshIoBroker::validate_output_file_index(size_t output_file_index) const
{
    ThrowErrorMsgIf(!is_index_valid(m_outputFiles, output_file_index),
                    "StkMeshIoBroker::validate_output_file_index: invalid output file index of "
                    << output_file_index << ".");

    ThrowErrorMsgIf (Teuchos::is_null(m_outputFiles[output_file_index]->get_output_io_region()),
                     "StkMeshIoBroker::validate_output_file_index: There is no Output mesh region associated with this output file index: " << output_file_index << ".");
}

void StkMeshIoBroker::validate_heartbeat_file_index(size_t heartbeat_file_index) const
{
    ThrowErrorMsgIf(!is_index_valid(m_heartbeat, heartbeat_file_index),
                    "StkMeshIoBroker::validate_heartbeat_file_index: invalid heartbeat file index of "
                    << heartbeat_file_index << ".");

    ThrowErrorMsgIf (Teuchos::is_null(m_heartbeat[heartbeat_file_index]->get_heartbeat_io_region()),
                     "StkMeshIoBroker::validate_heartbeat_file_index: There is no heartbeat mesh region associated with this heartbeat file index: " << heartbeat_file_index << ".");
}

void StkMeshIoBroker::validate_input_file_index(size_t input_file_index) const
{
    ThrowErrorMsgIf(!is_index_valid(m_inputFiles, input_file_index),
                    "StkMeshIoBroker::validate_input_file_index: invalid input file index of "
                    << input_file_index << ".");
}

void StkMeshIoBroker::add_field(size_t output_file_index, stk::mesh::FieldBase &field)
{
    add_field(output_file_index, field, field.name());
}

void StkMeshIoBroker::add_field(size_t output_file_index, stk::mesh::FieldBase &field, const std::string &alternate_name)
{
    add_field(output_file_index, field, field.entity_rank(), alternate_name);
}

void StkMeshIoBroker::add_field(size_t output_file_index, stk::mesh::FieldBase &field, stk::mesh::EntityRank var_type, const std::string &alternate_name)
{
    validate_output_file_index(output_file_index);
    OutputVariableParams var(alternate_name);
    m_outputFiles[output_file_index]->add_field(field, var, var_type);
}

void StkMeshIoBroker::add_field(size_t output_file_index, stk::mesh::FieldBase &field, stk::mesh::EntityRank var_type, const OutputVariableParams &var)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_field(field, var, var_type);
}

void StkMeshIoBroker::add_attribute_field(size_t output_file_index, stk::mesh::FieldBase &field, const OutputVariableParams &var)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_attribute_field(field, var);
}

void StkMeshIoBroker::add_user_data(size_t output_file_index, const std::vector<std::string> &parts, const std::string &alternate_name, stk::io::DataLocation loc)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_user_data(parts, alternate_name, loc);
}

bool StkMeshIoBroker::has_input_global(const std::string &globalVarName) const
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    return internal_has_global(region, globalVarName);
}

void StkMeshIoBroker::get_global_variable_names(std::vector<std::string> &names)
{
    validate_input_file_index(m_activeMeshIndex);
    m_inputFiles[m_activeMeshIndex]->get_global_variable_names(names);
}

bool StkMeshIoBroker::get_global(const std::string &globalVarName,
                                 boost::any &value, stk::util::ParameterType::Type type,
                                 bool abort_if_not_found)
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    return internal_read_parameter(region, globalVarName, value, type, abort_if_not_found);
}

size_t StkMeshIoBroker::get_global_variable_length(const std::string& globalVarName)
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();

    size_t length = 0;
    if (region->field_exists(globalVarName)) {
        Ioss::Field field = region->get_field(globalVarName);
        length = field.raw_count() * field.raw_storage()->component_count();
    }
    return length;
}

bool StkMeshIoBroker::get_global(const std::string &globalVarName, std::vector<double> &globalVar,
                                 bool abort_if_not_found)
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    return internal_read_global(region, globalVarName, globalVar, Ioss::Field::REAL,
                                abort_if_not_found);
}

bool StkMeshIoBroker::get_global(const std::string &globalVarName, std::vector<int> &globalVar,
                                 bool abort_if_not_found)
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    return internal_read_global(region, globalVarName, globalVar, Ioss::Field::INTEGER,
                                abort_if_not_found);
}

bool StkMeshIoBroker::get_global(const std::string &globalVarName, int &globalVar,
                                 bool abort_if_not_found)
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    return internal_read_global(region, globalVarName, globalVar, Ioss::Field::INTEGER,
                                abort_if_not_found);
}

bool StkMeshIoBroker::get_global(const std::string &globalVarName, double &globalVar,
                                 bool abort_if_not_found)
{
    validate_input_file_index(m_activeMeshIndex);
    auto region = m_inputFiles[m_activeMeshIndex]->get_input_io_region();
    return internal_read_global(region, globalVarName, globalVar, Ioss::Field::REAL,
                                abort_if_not_found);
}

bool StkMeshIoBroker::has_global(size_t output_file_index, const std::string &globalVarName) const
{
    validate_output_file_index(output_file_index);
    return m_outputFiles[output_file_index]->has_global(globalVarName);
}

void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &name,
                                 const boost::any &value, stk::util::ParameterType::Type type)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_global(name, value, type);
}

void StkMeshIoBroker::add_global_ref(size_t output_file_index, const std::string &name,
                                     const boost::any *value, stk::util::ParameterType::Type type)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_global_ref(name, value, type);
}

void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &globalVarName, Ioss::Field::BasicType dataType)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_global(globalVarName, dataType);
}

void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &globalVarName, int component_count, Ioss::Field::BasicType dataType)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_global(globalVarName, component_count, dataType);
}

void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &globalVarName, const std::string &storage, Ioss::Field::BasicType dataType)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->add_global(globalVarName, storage, dataType);
}

void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName,
                                   const boost::any &value, stk::util::ParameterType::Type type)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->write_global(globalVarName, value, type);
}

void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, double globalVarData)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->write_global(globalVarName, globalVarData);
}

void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, int globalVarData)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->write_global(globalVarName, globalVarData);
}

void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, std::vector<double>& globalVarData)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->write_global(globalVarName, globalVarData);
}

void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, std::vector<int>& globalVarData)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->write_global(globalVarName, globalVarData);
}

FieldNameToPartVector StkMeshIoBroker::get_nodal_var_names()
{
    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->get_var_names(Ioss::NODEBLOCK, meta_data());
}

FieldNameToPartVector StkMeshIoBroker::get_elem_var_names()
{
    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->get_var_names(Ioss::ELEMENTBLOCK, meta_data());
}

FieldNameToPartVector StkMeshIoBroker::get_nodeset_var_names()
{
    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->get_var_names(Ioss::NODESET, meta_data());
}

FieldNameToPartVector StkMeshIoBroker::get_sideset_var_names()
{
    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->get_var_names(Ioss::SIDESET, meta_data());
}

void StkMeshIoBroker::add_all_mesh_fields_as_input_fields(MeshField::TimeMatchOption tmo)
{
    validate_input_file_index(m_activeMeshIndex);
    m_inputFiles[m_activeMeshIndex]->add_all_mesh_fields_as_input_fields(meta_data(), tmo);
}

bool StkMeshIoBroker::read_input_field(stk::io::MeshField &mf)
{
    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->read_input_field(mf, bulk_data());
}

double StkMeshIoBroker::read_defined_input_fields(double time,
                                                  std::vector<stk::io::MeshField> *missingFields)
{
    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->read_defined_input_fields(time, missingFields, bulk_data());
}

double StkMeshIoBroker::read_defined_input_fields(int step,
                                                  std::vector<stk::io::MeshField> *missing)
{
    if (step <= 0) {
        return 0.0;
    }

    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->read_defined_input_fields(step, missing, bulk_data());
}

double StkMeshIoBroker::read_defined_input_fields_at_step(int step,
                                                          std::vector<stk::io::MeshField> *missing)
{
    if (step <= 0) {
        return 0.0;
    }

    validate_input_file_index(m_activeMeshIndex);
    return m_inputFiles[m_activeMeshIndex]->read_defined_input_fields_at_step(step, missing, bulk_data());
}

bool StkMeshIoBroker::use_nodeset_for_block_nodes_fields(size_t output_file_index) const
{
    validate_output_file_index(output_file_index);
    return m_outputFiles[output_file_index]->use_nodeset_for_block_nodes_fields();
}

void StkMeshIoBroker::use_nodeset_for_block_nodes_fields(size_t output_file_index, bool true_false)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->use_nodeset_for_block_nodes_fields(true_false);
}

bool StkMeshIoBroker::use_nodeset_for_sideset_nodes_fields(size_t output_file_index) const
{
    validate_output_file_index(output_file_index);
    return m_outputFiles[output_file_index]->use_nodeset_for_sideset_nodes_fields();
}

void StkMeshIoBroker::use_nodeset_for_sideset_nodes_fields(size_t output_file_index, bool true_false)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->use_nodeset_for_sideset_nodes_fields(true_false);
}

bool StkMeshIoBroker::use_nodeset_for_part_nodes_fields(size_t output_file_index) const
{
    validate_output_file_index(output_file_index);
    bool useNodesetForBlocks = m_outputFiles[output_file_index]->use_nodeset_for_block_nodes_fields();
    bool useNodesetForSidesets = m_outputFiles[output_file_index]->use_nodeset_for_sideset_nodes_fields();

    return (useNodesetForBlocks || useNodesetForSidesets);
}

void StkMeshIoBroker::use_nodeset_for_part_nodes_fields(size_t output_file_index, bool true_false)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->use_nodeset_for_block_nodes_fields(true_false);
    m_outputFiles[output_file_index]->use_nodeset_for_sideset_nodes_fields(true_false);
}

bool StkMeshIoBroker::check_field_existence_when_creating_nodesets(size_t output_file_index) const
{
    validate_output_file_index(output_file_index);
    return m_outputFiles[output_file_index]->check_field_existence_when_creating_nodesets();
}

void StkMeshIoBroker::check_field_existence_when_creating_nodesets(size_t output_file_index, bool true_false)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->check_field_existence_when_creating_nodesets(true_false);
}

bool StkMeshIoBroker::use_part_id_for_output(size_t output_file_index) const
{
    validate_output_file_index(output_file_index);
    return m_outputFiles[output_file_index]->use_part_id_for_output();
}

void StkMeshIoBroker::use_part_id_for_output(size_t output_file_index, bool true_false)
{
    validate_output_file_index(output_file_index);
    m_outputFiles[output_file_index]->use_part_id_for_output(true_false);
}

void StkMeshIoBroker::set_option_to_not_collapse_sequenced_fields()
{
    property_add(Ioss::Property("ENABLE_FIELD_RECOGNITION", "NO"));
}

int StkMeshIoBroker::get_num_time_steps()
{
    int numTimeSteps = 0;
    Ioss::Region *ioRegion = get_input_io_region().get();
    if(ioRegion != nullptr) {
        Ioss::Property stateCount = ioRegion->get_implicit_property("state_count");
        numTimeSteps = stateCount.get_int();
    }
    return numTimeSteps;
}

std::vector<double> StkMeshIoBroker::get_time_steps()
{
    int numTimeSteps = get_num_time_steps();
    std::vector<double> timeSteps;

    Ioss::Region *ioRegion = get_input_io_region().get();
    if(ioRegion != nullptr)  {
        for(int istep = 0; istep < numTimeSteps; istep++) {
            double state_time = ioRegion->get_state_time(istep + 1);
            timeSteps.push_back(state_time);
        }
    }
    return timeSteps;
}

double StkMeshIoBroker::get_max_time()
{
    return get_input_io_region()->get_max_time().second;
}

void StkMeshIoBroker::set_max_num_steps_before_overwrite(size_t outputFileIndex, int maxNumStepsInFile)
{
    get_output_io_region(outputFileIndex)->get_database()->set_cycle_count(maxNumStepsInFile);
}

size_t StkMeshIoBroker::add_heartbeat_output(const std::string &filename, HeartbeatType hb_type,
                                             const Ioss::PropertyManager &properties, bool openFileImmediately)
{
    std::string out_filename = filename;
    stk::util::filename_substitution(out_filename);
    auto heartbeat = Teuchos::rcp(new impl::Heartbeat(out_filename, hb_type,
                                                      properties, m_communicator, openFileImmediately));
    m_heartbeat.push_back(heartbeat);
    return m_heartbeat.size()-1;
}

int StkMeshIoBroker::check_integer_size_requirements()
{
    // 1. If the INTEGER_SIZE_DB or _API property exists, then use its value no matter what...
    if (m_propertyManager.exists("INTEGER_SIZE_DB")) {
        return m_propertyManager.get("INTEGER_SIZE_DB").get_int();
    }

    if (m_propertyManager.exists("INTEGER_SIZE_API")) {
        return m_propertyManager.get("INTEGER_SIZE_API").get_int();
    }

    // 2. If input_region exists, then if it is using 64-bit integers, the output should
    //    use those also.
    Ioss::Region *input_region = nullptr;
    if (is_index_valid(m_inputFiles, m_activeMeshIndex)) {
        input_region = get_input_io_region().get();
    }
    if (input_region != nullptr) {
        // Get the integer size setting for the database associated with the region.
        int int_size = db_api_int_size(input_region);
        if (int_size == 8) {
            return int_size;
        }
    }

    // 3. If any entity count exceeds INT_MAX, then use 64-bit integers.
    if ( !Teuchos::is_null(m_bulkData) ) {
        std::vector<size_t> entityCounts;
        stk::mesh::comm_mesh_counts(*m_bulkData, entityCounts);
        for (size_t i=0; i < entityCounts.size(); i++) {
            if (entityCounts[i] > (size_t)std::numeric_limits<int>::max()) {
                return 8;
            }
        }
    }

    // 4. Should also check if the maximum node or element id exceeds INT_MAX.

    // 5. Default to 4-byte integers...
    return 4;
}

void StkMeshIoBroker::set_name_and_version_for_qa_record(size_t outputFileIndex, const std::string &codeName, const std::string &codeVersion)
{
    Ioss::Region *region = get_output_io_region(outputFileIndex).get();
    region->property_add(Ioss::Property(std::string("code_name"), codeName));
    region->property_add(Ioss::Property(std::string("code_version"), codeVersion));
}

void StkMeshIoBroker::add_qa_records(size_t outputFileIndex, const std::vector<QaRecord> &qaRecords)
{
    Ioss::Region *region = get_output_io_region(outputFileIndex).get();
    for(const QaRecord &qaRec : qaRecords)
        region->add_qa_record(qaRec.name, qaRec.version, qaRec.date, qaRec.time);
}

void StkMeshIoBroker::add_info_records(size_t outputFileIndex, const std::vector<std::string> &infoRecords)
{
    Ioss::Region *region = get_output_io_region(outputFileIndex).get();
    region->add_information_records(infoRecords);
}

std::vector<QaRecord> StkMeshIoBroker::get_qa_records()
{
    std::vector<QaRecord> qaRecords;
    Ioss::Region *region = get_input_io_region().get();
    const std::vector<std::string> &qa = region->get_qa_records();
    for (size_t i = 0; i < qa.size(); i += 4)
        qaRecords.push_back({qa[i + 0], qa[i + 1], qa[i + 2], qa[i + 3]});
    return qaRecords;
}

std::vector<std::string> StkMeshIoBroker::get_info_records()
{
    Ioss::Region *region = get_input_io_region().get();
    return region->get_information_records();
}

stk::mesh::FieldVector StkMeshIoBroker::get_ordered_attribute_fields(const stk::mesh::Part *blockPart) const
{
    stk::mesh::FieldVector attrFields;
    if(blockPart->mesh_meta_data_ordinal() < attributeFieldOrderingByPartOrdinal.size()) {
        const std::vector<int> &fieldOrds = attributeFieldOrderingByPartOrdinal[blockPart->mesh_meta_data_ordinal()];
        attrFields.resize(fieldOrds.size());
        const stk::mesh::FieldVector &allFields = m_metaData->get_fields();
        for(size_t i=0; i<fieldOrds.size(); i++) {
            attrFields[i] = allFields[fieldOrds[i]];
        }
    }
    return attrFields;
}

const std::vector<std::vector<int>> & StkMeshIoBroker::get_attribute_field_ordering_stored_by_part_ordinal() const
{
    return attributeFieldOrderingByPartOrdinal;
}

void StkMeshIoBroker::set_attribute_field_ordering_stored_by_part_ordinal(const std::vector<std::vector<int>> &ordering)
{
    attributeFieldOrderingByPartOrdinal = ordering;
}

void StkMeshIoBroker::fill_coordinate_frames(std::vector<int>& ids, std::vector<double>& coords, std::vector<char>& tags)
{
    Ioss::Region *ioregion = get_input_io_region().get();
    const Ioss::CoordinateFrameContainer& coordFrames = ioregion->get_coordinate_frames();

    size_t nFrames = coordFrames.size();

    ids.resize(nFrames);
    const int coordSize = 9;
    coords.resize(coordSize*nFrames);
    tags.resize(nFrames);

    for (size_t i=0;i<nFrames;i++) {
        tags[i] = coordFrames[i].tag();
        ids[i] = coordFrames[i].id();
        memcpy(&coords[coordSize*i],coordFrames[i].coordinates(),sizeof(double)*coordSize);
    }
}

Ioss::DatabaseIO *StkMeshIoBroker::get_input_database(size_t input_index)
{
    if(is_input_index_valid(input_index)) {
        return m_inputFiles[input_index]->get_input_database().get();
    }

    return nullptr;
}

Ioss::DatabaseIO *StkMeshIoBroker::get_output_database(size_t output_index)
{
    if(is_output_index_valid(output_index)) {
        return m_outputFiles[output_index]->get_output_database();
    }

    return nullptr;
}

bool StkMeshIoBroker::set_input_multistate_suffixes(size_t input_index, std::vector<std::string>& multiStateSuffixes)
{
    if(is_input_index_valid(input_index)) {
        return m_inputFiles[input_index]->set_multistate_suffixes(multiStateSuffixes);
    }

    return false;
}

bool StkMeshIoBroker::set_output_multistate_suffixes(size_t output_index, std::vector<std::string>& multiStateSuffixes)
{
    if(is_output_index_valid(output_index)) {
        return m_outputFiles[output_index]->set_multistate_suffixes(multiStateSuffixes);
    }

    return false;
}

void StkMeshIoBroker::set_reference_input_region(size_t outputIndex, StkMeshIoBroker& inputBroker)
{
    validate_output_file_index(outputIndex);

    const Ioss::Region *input_region = inputBroker.get_input_io_region().get();
    m_outputFiles[outputIndex]->set_input_region(input_region);
}

bool StkMeshIoBroker::create_named_suffix_field_type(const std::string& type_name, std::vector<std::string>& suffices) const
{
    return Ioss::VariableType::create_named_suffix_field_type(type_name, suffices);
}

bool StkMeshIoBroker::add_field_type_mapping(const std::string& field, const std::string& type) const
{
    return Ioss::VariableType::add_field_type_mapping(field, type);
}

bool StkMeshIoBroker::has_heartbeat_global(size_t heartbeat_file_index, const std::string &globalVarName) const
{
    validate_heartbeat_file_index(heartbeat_file_index);
    return m_heartbeat[heartbeat_file_index]->has_global(globalVarName);
}

size_t StkMeshIoBroker::get_heartbeat_global_component_count(size_t heartbeat_file_index, const std::string &globalVarName) const
{
    validate_heartbeat_file_index(heartbeat_file_index);

    size_t comp_count = 0;
    if(has_heartbeat_global(heartbeat_file_index, globalVarName)) {
        Teuchos::RCP<Ioss::Region> hbRegion = m_heartbeat[heartbeat_file_index]->get_heartbeat_io_region();

        if(!hbRegion.is_null()) {
            comp_count = hbRegion->get_fieldref(globalVarName).raw_storage()->component_count();
        }
    }

    return comp_count;
}

std::vector<stk::mesh::Entity> StkMeshIoBroker::get_output_entities(size_t output_index,
                                                                    const stk::mesh::BulkData& bulk_data,
                                                                    const std::string &name)
{
    std::vector<stk::mesh::Entity> entities;

    if(is_output_index_valid(output_index)) {
        entities = m_outputFiles[output_index]->get_output_entities(bulk_data, name);
    }

    return entities;
}




} // namespace io
} // namespace stk
