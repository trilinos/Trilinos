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

#ifndef STK_IO_OUTPUTFILE_HPP
#define STK_IO_OUTPUTFILE_HPP
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_Field.h>                     // for Field::BasicType, Field
#include <exception>                 // for exception
#include <cstddef>                          // for size_t
#include <stk_io/DatabasePurpose.hpp>       // for DatabasePurpose
#include <stk_io/IOHelpers.hpp>             // for impl::add_global
#include <stk_io/IossBridge.hpp>            // for GlobalAnyVariable
#include <stk_mesh/base/Selector.hpp>       // for Selector
#include <stk_util/util/ParameterList.hpp>  // for Type
#include <string>                           // for string
#include <utility>                          // for pair, swap
#include <vector>                           // for vector
#include <any>
#include "mpi.h"                            // for MPI_Comm, ompi_communicat...
#include "stk_io/FieldAndName.hpp"          // for FieldAndName, UserDataAnd...
#include "stk_io/OutputVariableParams.hpp"  // for OutputVariableParams
#include "stk_io/StkIoUtils.hpp"            // for get_io_parameter_size_and...
#include "stk_mesh/base/FieldState.hpp"     // for FieldState
#include "stk_mesh/base/Types.hpp"          // for EntityRank
#include "stk_topology/topology.hpp"        // for topology, topology::ELEM_...
namespace Ioss { class DatabaseIO; }
namespace Ioss { class PropertyManager; }
namespace Ioss { class Region; }
namespace stk { namespace io { struct OutputParams; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct Entity; } }
namespace Ioss { class Property; }
namespace stk { namespace io { class InputFile; } }
namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk {
namespace io {
namespace impl
{
class OutputFile
{
public:
    OutputFile(const std::string &filename, MPI_Comm communicator, DatabasePurpose db_type,
               Ioss::PropertyManager& property_manager, const Ioss::Region *input_region,
               char const* type = "exodus", bool openFileImmediately = true)
    : m_currentOutputStep(-1),
      m_useNodesetForBlockNodesFields(false),
      m_useNodesetForSidesetNodesFields(false),
      m_checkFieldExistenceWhenCreatingNodesets(true),
      m_usePartIdForOutput(true),
      m_meshDefined(false),
      m_fieldsDefined(false),
      m_anyGlobalVariablesDefined(false),
      m_appendingToMesh(false),
      m_hasGhosting(false),
      m_hasAdaptivity(false),
      m_isSkinMesh(false),
      m_enableEdgeIO(false),
      m_dbPurpose(db_type),
      m_inputRegion(input_region),
      m_outputParams(nullptr),
      m_subsetSelector(nullptr),
      m_sharedSelector(nullptr),
      m_skinMeshSelector(nullptr),
      m_multiStateSuffixes(nullptr),
      m_filterEmptyEntityBlocks(false),
      m_filterEmptyAssemblyEntityBlocks(false)
    {
        initialize_output_selectors();
        setup_output_file(filename, communicator, property_manager, type, openFileImmediately);

        initialize_output_selectors();
    }

    OutputFile(std::shared_ptr<Ioss::Region> ioss_output_region, MPI_Comm communicator,
               DatabasePurpose db_type, const Ioss::Region *input_region)
    : m_currentOutputStep(-1),
      m_useNodesetForBlockNodesFields(false),
      m_useNodesetForSidesetNodesFields(false),
      m_checkFieldExistenceWhenCreatingNodesets(true),
      m_usePartIdForOutput(true),
      m_meshDefined(false),
      m_fieldsDefined(false),
      m_anyGlobalVariablesDefined(false),
      m_appendingToMesh(false),
      m_hasGhosting(false),
      m_hasAdaptivity(false),
      m_isSkinMesh(false),
      m_enableEdgeIO(false),
      m_dbPurpose(db_type),
      m_inputRegion(input_region),
      m_outputParams(nullptr),
      m_subsetSelector(nullptr),
      m_sharedSelector(nullptr),
      m_skinMeshSelector(nullptr),
      m_multiStateSuffixes(nullptr),
      m_filterEmptyEntityBlocks(false),
      m_filterEmptyAssemblyEntityBlocks(false)
    {
        m_region = ioss_output_region;
        m_meshDefined = true;
        initialize_output_selectors();
    }

    std::shared_ptr<Ioss::Region> get_output_ioss_region() {
        return m_region;
    }
    ~OutputFile();

    void set_input_region(const Ioss::Region *input_region);

    void setup_output_params(OutputParams &params) const;

    bool set_multistate_suffixes(const std::vector<std::string>& multiStateSuffixes);

    void write_output_mesh(const stk::mesh::BulkData& bulk_data,
                           const std::vector<std::vector<int>> &attributeOrdering);
    void flush_output() const;
    void add_field(stk::mesh::FieldBase &field, const std::string &alternate_name, stk::mesh::EntityRank var_type);
    void add_field(stk::mesh::FieldBase &field, const OutputVariableParams &var, stk::mesh::EntityRank var_type);
    void add_attribute_field(stk::mesh::FieldBase &field, const OutputVariableParams &var);
    void add_user_data(const std::vector<std::string>& userData, const std::string &alternate_name, stk::io::DataLocation loc);
    bool has_global(const std::string &globalVarName) const;
    void add_global(const std::string &variableName, const stk::util::Parameter &value);

    template<typename T>
    void add_global(const std::string &variableName, const T& value, stk::util::ParameterType::Type type)
    {
      std::any anyValue(value);
      std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, anyValue);
      m_anyGlobalVariablesDefined = true;
      impl::add_global(m_region, variableName, parameter_type.first, parameter_type.second);
    }

    void add_global_ref(const std::string &variableName, const stk::util::Parameter &value);
    void add_global(const std::string &variableName, Ioss::Field::BasicType dataType);
    void add_global(const std::string &variableName, const std::string &type, Ioss::Field::BasicType dataType);
    void add_global(const std::string &variableName, int component_count,     Ioss::Field::BasicType dataType);

    void write_global(const std::string &variableName,
                      const stk::util::Parameter &param);
    void write_global(const std::string &variableName, double globalVarData);
    void write_global(const std::string &variableName, int globalVarData);
    void write_global(const std::string &variableName, std::vector<double>& globalVarData);
    void write_global(const std::string &variableName, std::vector<int>& globalVarData);

    void begin_output_step(double time, const stk::mesh::BulkData& bulk_data, const std::vector<std::vector<int>> &attributeOrdering);
    void end_output_step();

    int write_defined_output_fields(const stk::mesh::BulkData &bulk_data, const stk::mesh::FieldState *state = nullptr);

    int process_output_request(double time, const stk::mesh::BulkData& bulk_data, const std::vector<std::vector<int>> &attributeOrdering);

    void set_subset_selector(std::shared_ptr<stk::mesh::Selector> my_selector);
    void set_shared_selector(std::shared_ptr<stk::mesh::Selector> my_selector);
    void set_skin_mesh_selector(std::shared_ptr<stk::mesh::Selector> my_selector);

    void set_output_selector(stk::topology::rank_t rank, std::shared_ptr<stk::mesh::Selector> my_selector);

    bool use_nodeset_for_block_nodes_fields() const;
    void use_nodeset_for_block_nodes_fields(bool flag);

    bool use_nodeset_for_sideset_nodes_fields() const;
    void use_nodeset_for_sideset_nodes_fields(bool flag);

    bool check_field_existence_when_creating_nodesets() const;
    void check_field_existence_when_creating_nodesets(bool flag);

    bool use_part_id_for_output() const;
    void use_part_id_for_output(bool flag);

    bool has_ghosting() const;
    void has_ghosting(bool hasGhosting);

    bool has_adaptivity() const;
    void has_adaptivity(bool hasAdaptivity);

    bool is_skin_mesh() const;

    void set_enable_edge_io(bool enableEdgeIO);

    void set_filter_empty_entity_blocks(const bool filterEmptyEntityBlocks);
    void set_filter_empty_assembly_entity_blocks(const bool filterEmptyAssemblyEntityBlocks);

    Ioss::DatabaseIO *get_output_database();

    std::vector<stk::mesh::Entity> get_output_entities(const stk::mesh::BulkData& bulk_data, const std::string &name);

private:
    void define_output_fields(const stk::mesh::BulkData& bulk_data, const std::vector<std::vector<int>> &attributeOrdering);
    void setup_output_file(const std::string &filename, MPI_Comm communicator,
                           Ioss::PropertyManager &property_manager,
                           char const* type = "exodus", bool openFileImmediately = true);
    void initialize_output_selectors()
    {
        m_outputSelector[stk::topology::NODE_RANK].reset();
        m_outputSelector[stk::topology::EDGE_RANK].reset();
        m_outputSelector[stk::topology::FACE_RANK].reset();
        m_outputSelector[stk::topology::ELEM_RANK].reset();
    }

    int m_currentOutputStep;
    bool m_useNodesetForBlockNodesFields;
    bool m_useNodesetForSidesetNodesFields;
    bool m_checkFieldExistenceWhenCreatingNodesets;
    bool m_usePartIdForOutput;
    bool m_meshDefined;
    bool m_fieldsDefined;
    bool m_anyGlobalVariablesDefined;
    bool m_appendingToMesh;
    bool m_hasGhosting;
    bool m_hasAdaptivity;
    bool m_isSkinMesh;
    bool m_enableEdgeIO;
    DatabasePurpose m_dbPurpose;
    const Ioss::Region* m_inputRegion;
    std::shared_ptr<stk::io::OutputParams> m_outputParams;
    std::shared_ptr<stk::mesh::Selector> m_subsetSelector;
    std::shared_ptr<stk::mesh::Selector> m_sharedSelector;
    std::shared_ptr<stk::mesh::Selector> m_outputSelector[stk::topology::ELEM_RANK+1];
    std::shared_ptr<stk::mesh::Selector> m_skinMeshSelector;
    std::shared_ptr<Ioss::Region> m_region;
    std::vector<stk::io::FieldAndName> m_namedFields;
    std::vector<stk::io::FieldAndName> m_additionalAttributeFields;
    std::vector<stk::io::UserDataAndName> m_userData;

    // Global fields that can be output automatically without app calling write_global.
    std::vector<GlobalAnyVariable> m_globalAnyFields;

    std::vector<std::string>* m_multiStateSuffixes = nullptr;

    bool m_filterEmptyEntityBlocks;
    bool m_filterEmptyAssemblyEntityBlocks;

    OutputFile(const OutputFile &);
    const OutputFile & operator=(const OutputFile &);
};
}
}
}
#endif
