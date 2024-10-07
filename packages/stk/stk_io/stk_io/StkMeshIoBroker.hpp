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

#ifndef STK_IO_STKMESHIOBROKER_HPP
#define STK_IO_STKMESHIOBROKER_HPP
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_Field.h>                     // for Field, Field::REDUCTION
#include <Ioss_PropertyManager.h>           // for PropertyManager
      // file exists and is readable and will throw an exception if not.
#include <cstddef>                          // for size_t
#include <stk_io/DatabasePurpose.hpp>       // for DatabasePurpose
#include <stk_io/Heartbeat.hpp>             // for Heartbeat, HeartbeatType
#include <stk_io/IossBridge.hpp>            // for STKIORequire, FieldNameTo...
#include <stk_io/MeshField.hpp>             // for MeshField, MeshField::CLO...
#include <stk_io/OutputFile.hpp>            // for OutputFile
#include <stk_mesh/base/Selector.hpp>       // for Selector
#include <stk_util/parallel/Parallel.hpp>   // for ParallelMachine
#include <stk_util/util/ParameterList.hpp>  // for Parameter, Type
#include <memory>
#include <string>                           // for string
#include <utility>                          // for move, swap
#include <vector>                           // for vector
#include "mpi.h"                            // for ompi_communicator_t
#include "stk_io/FieldAndName.hpp"          // for DataLocation
#include "stk_io/OutputVariableParams.hpp"  // for OutputVariableParams
#include "stk_mesh/base/FieldState.hpp"     // for FieldState
#include "stk_mesh/base/Types.hpp"          // for EntityRank, FieldVector
#include "stk_topology/topology.hpp"        // for topology, topology::rank_t
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, ThrowRequire...
namespace Ioss { class DatabaseIO; }
namespace Ioss { class Property; }
namespace Ioss { class Region; }
namespace stk { namespace io { class InputFile; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class SidesetUpdater; } }

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace stk { namespace mesh { class BulkData; } }


namespace stk {
  namespace io {
    
    struct QaRecord
    {
        std::string name;
        std::string version;
        std::string date;
        std::string time;
    };

    //-BEGIN
    class StkMeshIoBroker {
    public:

      // \param[in] comm MPI Communicator to be used for all
      // parallel communication needed to generate the mesh.
      StkMeshIoBroker(stk::ParallelMachine comm);
      StkMeshIoBroker();

      virtual ~StkMeshIoBroker();

      // Add the specified 'property' to the default property
      // manager for this StkMeshIoBroker object.  The property
      // manager is (currently) used for the input mesh region and
      // the output mesh regions if no output property manager is
      // specified in the create_output_mesh() call.
      //
      // If the property already exists, it is overwritten by the
      // current value.
      void property_add(const Ioss::Property &property);

      // If a property with the specified name exists, remove it.
      // (Function provided for backward compatibility with existing
      // use in Percept following change to make m_property_manger
      // private)
      void remove_property_if_exists(const std::string &property_name);

      bool property_exists(const std::string &property_name) const;
      void copy_property(const StkMeshIoBroker& src_broker, const std::string &property_name);

      std::shared_ptr<Ioss::Region> get_input_ioss_region() const;
      std::shared_ptr<Ioss::Region> get_output_ioss_region(size_t output_file_index) const;
      std::shared_ptr<Ioss::Region> get_heartbeat_ioss_region(size_t heartbeat_file_index) const;

      void begin_define_transient_for_heartbeat(size_t heartbeat_file_index);
      void end_define_transient_for_heartbeat(size_t heartbeat_file_index);

      // The default entity rank names are "NODE", "EDGE", "FACE",
      // "ELEMENT" If an application wants to change these names
      // and/or add additional names, it should call this function
      // with a vector containing the desired names.  To use the
      // default names, but add additional names, call
      // stk::mesh::entity_rank_names() which will return a reference
      // to the default entity name vector; copy these and add the new
      // names and then call this function.
      //
      // Must be called prior to create_input_mesh() or it will have
      // no effect
      //
      // @param rank_names -- std::vector<std::string> containing the
      // names of your entities.
      void set_rank_name_vector(const std::vector<std::string> &rank_names);

      // The selector is an optional selector that can be provided
      // which will be used to filter the mesh entities (node,
      // element, face, edge) associated with an
      // Ioss Entity (element block, node block, node set, ...). The
      // optional selector will be 'anded' with the normal selector
      // (typically locally owned part) used to associate entities
      // when generating the output database.
      void set_subset_selector(size_t output_file_index, std::shared_ptr<stk::mesh::Selector> my_selector);
      void set_subset_selector(size_t output_file_index, const stk::mesh::Selector &my_selector);

      void set_skin_mesh_selector(size_t output_file_index, std::shared_ptr<stk::mesh::Selector> my_selector);
      void set_skin_mesh_selector(size_t output_file_index, stk::mesh::Selector &my_selector);

      void set_shared_selector(size_t output_file_index, std::shared_ptr<stk::mesh::Selector> my_selector);
      void set_shared_selector(size_t output_file_index, stk::mesh::Selector &my_selector);

      void set_output_selector(size_t output_file_index, stk::topology::rank_t rank, std::shared_ptr<stk::mesh::Selector> my_selector);
      void set_output_selector(size_t output_file_index, stk::topology::rank_t rank, stk::mesh::Selector &my_selector);

      void set_ghosting_filter(size_t output_file_index, bool hasGhosting);
      void set_adaptivity_filter(size_t output_file_index, bool hasAdaptivity);

      void set_filter_empty_output_entity_blocks(size_t output_file_index, const bool filterEmptyEntityBlocks);
      void set_filter_empty_output_assembly_entity_blocks(size_t output_file_index, const bool filterEmptyAssemblyEntityBlocks);

      stk::mesh::Selector get_active_selector() const;
      void set_active_selector(stk::mesh::Selector my_selector);

      // Override the default MeshBuilder if you need to have a custom STK Mesh
      // auto-generated internally.
      void set_mesh_builder(std::shared_ptr<stk::mesh::MeshBuilder> meshBuilder);

      // Set bulk data directly with your own bulk data. If this is
      // not called prior to the populate_bulk_data() call, it will be
      // created automatically using the communicator of the m_input_region.
      // If meta data is not already set, then set the meta data from the
      // bulk data's metadata
      void set_bulk_data(std::shared_ptr<stk::mesh::BulkData> arg_bulk_data);
      void set_bulk_data(stk::mesh::BulkData &arg_bulk_data);

      // Replace the current bulk data directly with your own bulk data.
      // There must be a current bulk data and the current meta data
      // must match the meta data associated with the new bulk data.
      // This is a potentially dangerous call depending on what
      // point of the usage it is made.  Typical use would be to call
      // this, only if needed, after you are completely done accessing
      // the input mesh and before any access to the output mesh and only
      // if the output mesh needs a different bulk data than the input mesh.
      void replace_bulk_data(std::shared_ptr<stk::mesh::BulkData> arg_bulk_data);
      void replace_bulk_data(stk::mesh::BulkData &arg_bulk_data);

      enum SideSetFaceCreationBehavior {
          STK_IO_SIDESET_FACE_CREATION_CLASSIC = 42,
          STK_IO_SIDESET_FACE_CREATION_CURRENT = 73,
          STK_IO_SIDE_CREATION_USING_GRAPH_TEST = 99
      };

      void set_sideset_face_creation_behavior(SideSetFaceCreationBehavior behavior)
      {
          STK_ThrowRequireWithSierraHelpMsg(behavior!=STK_IO_SIDE_CREATION_USING_GRAPH_TEST);
          m_sidesetFaceCreationBehavior = behavior;
      }

      void set_auto_load_attributes(bool shouldAutoLoadAttributes)
      {
          m_autoLoadAttributes = shouldAutoLoadAttributes;
      }

      void set_auto_load_distribution_factor_per_nodeset(bool shouldAutoLoad)
      {
          m_autoLoadDistributionFactorPerNodeSet = shouldAutoLoad;
      }

      void cache_entity_list_for_transient_steps(bool cacheEntityList)
      {
          m_cacheEntityListForTransientSteps = cacheEntityList;
      }

      bool get_filter_empty_input_entity_blocks() const;
      bool get_filter_empty_input_entity_blocks(size_t input_file_index) const;

      // Create the Ioss::DatabaseIO associated with the specified filename
      // and type (exodus by default). The routine checks that the
      // file exists and is readable and will throw an exception if not.
      // No meta_data or bulk_data is created at this time, but the
      // get_input_io_database() will be valid.
      //
      // The Ioss::DatabaseIO can be accessed via the
      // get_input_io_database() method if you need to set some
      // options on the database prior to it being read.
      //
      // \param[in] filename If the mesh type is file based ("exodus"),
      // then this contains the full pathname to the file containing the
      // mesh information.  If the mesh type is a generated type, then
      // this parameter contains data used by the generation routines.
      // See the GeneratedMesh documentation.
      //
      // \param[in] type The format of the mesh that will be "read".
      // Valid types are "exodus", "generated", "pamgen".
      //
      // \param[in] type The format of the mesh that will be "read".
      // Valid types are "exodus", "generated", "pamgen".
      //
      // \param[in] purpose The purpose of the mesh database.
      // An input of type READ_RESTART will read all-but-one of the 
      //    states of a multi-state field,
      // An input of type READ_MESH will only read the newest state 
      //    of a multi-state field.
      // Other behavioral differences may be added in the future 
      //    (e.g., dealing with adaptivity...)
      size_t add_mesh_database(const std::string &filename,
			       const std::string &type,
			       DatabasePurpose purpose);

      
      size_t add_mesh_database(const std::string &filename,
                               const std::string &type,
                               DatabasePurpose purpose,
                               const Ioss::PropertyManager &properties);

      // Create Ioss::DatabaseIO associated with the specified
      // filename using the default filetype (typically "exodus"). If
      // the filename is prepended with a type followed by a colon
      // (e.g., "generated:10x10x10"), then use that type.  Valid
      // input types are: exodus, dof, pamgen, and possibly others.
      //
      // \param[in] filename If the mesh type is file based ("exodus"),
      // then this contains the full pathname to the file containing the
      // mesh information.  If the mesh type is a generated type, then
      // this parameter contains data used by the generation routines.
      // Optionally prepended by a filetype and a colon.
      size_t add_mesh_database(std::string filename,
			       DatabasePurpose purpose);

      // Set the input Ioss::Region directly instead of letting it be
      // created by StkMeshIoBroker during the create_input_mesh(type,
      // filename) call. After setting the input io region, you would
      // then either set the metadata manually using the
      // set_meta_data() call, or call the no-argument
      // create_input_mesh() function which will then create a meta
      // data corresponding to the data in the Ioss::Region.
      //
      // [2013-11-13: GDS: Currently
      // only used in Salinas/tools/superelem/MkSuperStkMesh.C:
      // The use-case is adding new parts to a mesh]
      size_t add_mesh_database(std::shared_ptr<Ioss::Region> ioss_input_region);

      // Get a reference to an existing mesh database so it can be modified
      // Typical modifications deal with
      // times: tart/stop/offset/scale/cycle/periodlength.
      InputFile &get_mesh_database(size_t input_file_index);

      // Remove the specified mesh database from the list of mesh databases.
      // All files associated with the mesh database will be closed and destructors
      // run on all classes associated with the database. An exeception will be thrown
      // if the specified input_file_index does not refer to a valid mesh database.
      void remove_mesh_database(size_t input_file_index);

      size_t set_active_mesh(size_t input_file_index);
      size_t get_active_mesh() const {return m_activeMeshIndex;}
      
      // Read/Generate the metadata for mesh of the specified type. By
      // default, all entities in the mesh (nodeblocks, element blocks,
      // nodesets, sidesets) will have an associated stk mesh part
      // created for it.
      //
      // If the input region was previously defined by calling
      // 'set_input_io_region()', then that Ioss::Region will be used
      // instead of opening the file and creating a new Ioss::Region.
      //
      // Following this call, the 'populate_bulk_data()' function should
      // be called to read the bulk data from the mesh and generate the
      // corresponding stk mesh entities (nodes, elements, faces, ...)
      //
      // The meta_data will not be committed by this function.
      void create_input_mesh();

      // Read/Generate the bulk data for the mesh.  This function will
      // create all stk mesh entities (nodes, elements) with the
      // correct nodeal coordinates, element connectivity, element
      // attribute data, and nodeset and sideset membership.
      //
      // NOTE: this function will commit the meta data if it hasn't
      // already been committed.
      //
      // NOTE: this function internally calls the two methods
      // 'populate_mesh()' and 'populate_field_data()', declared
      // below, and does NOT do the delayed field-data allocation
      // optimization.
      void populate_bulk_data();

      // Read/generate the bulk data for the mesh, delaying (by default)
      // the allocation of any field-data such as coordinates, attributes
      // and distribution factors.
      // This field-data is read/generated by separately calling the
      // 'populate_field_data()' method declared below.
      // Note that the above-declared 'populate_bulk_data()' method
      // calls both of these methods.
      virtual void populate_mesh(bool delay_field_data_allocation = true);
      bool populate_mesh_elements_and_nodes(bool delay_field_data_allocation);
      void populate_mesh_entitysets(bool i_started_modification_cycle);

      // Read/generate the field-data for the mesh, including
      // coordinates, attributes and distribution factors.
      // Note that this call should be preceded by a call to
      // the above-declared 'populate_mesh()' method.
      void populate_field_data();

      FieldNameToPartVector get_nodal_var_names() const;
      FieldNameToPartVector get_elem_var_names() const;
      FieldNameToPartVector get_nodeset_var_names() const;
      FieldNameToPartVector get_sideset_var_names() const;

      // For all transient fields on the mesh database:
      // - declare a stk_field of the same name,
      // - add the field as in input field for the mesh 
      //        (call add_input_field()).
      //
      // Then, when 'read_defined_input_fields()' is called, the data on
      // the mesh will be read and will populate the stk field.
      //
      // NOTE: This should only be called if you want *ALL* transient
      // fields on the input mesh to be defined as stk fields and read.
      void add_all_mesh_fields_as_input_fields(MeshField::TimeMatchOption tmo = MeshField::CLOSEST);

      // For all transient input fields defined, read the data at the
      // specified database step 'step' (1-based) and populate the stk
      // data structures with those values.
      // If 'missing' is non-NULL, then any fields that are not found
      // on the input database will be put on the vector.  If 'missing'
      // is NULL, then an exception will be thrown if any fields are
      // not found.
      double read_defined_input_fields(int step,
				       std::vector<stk::io::MeshField> *missing=nullptr);
      double read_defined_input_fields_at_step(int step,
                                       std::vector<stk::io::MeshField> *missing=nullptr);
      // For all transient input fields defined, read the data at the
      // specified database time 'time' and populate the stk data
      // structures with those values.  
      //
      // If the MeshField specifies "CLOSEST" option, then the
      // database time closest to the specified time will be used; if
      // the MeshField specifies "LINEAR_INTERPOLATION" option, then
      // the field values will be interpolated based on the two
      // surrounding times on the database; if the time is less than
      // the minimum time on the database or greater than the maximum
      // time, then the field values at those extremes will be
      // returned (no extrapolation).
      // 
      // If 'missing' is non-NULL, then any fields that are not found
      // on the input database will be put on the vector.  If
      // 'missing' is NULL, then an exception will be thrown if any
      // fields are not found.
      double read_defined_input_fields(double time,
				       std::vector<stk::io::MeshField> *missing=nullptr);

      bool read_input_field(stk::io::MeshField &mf);
      bool read_input_field(stk::io::MeshField &mf, stk::io::FieldReadStatus &readStatus);
      
      void get_global_variable_names(std::vector<std::string> &names) const;
      size_t get_global_variable_length(const std::string& name) const;

      bool get_global(const std::string &variableName,
              stk::util::Parameter &param,
              bool abort_if_not_found=true) const;
      bool get_global(const std::string &variableName,
              int &globalVar,
              bool abort_if_not_found=true) const;
      bool get_global(const std::string &variableName,
              double &globalVar,
              bool abort_if_not_found=true) const;
      bool get_global(const std::string &variableName,
              std::vector<double> &globalVar,
              bool abort_if_not_found=true) const;
      bool get_global(const std::string &variableName,
              std::vector<int> &globalVar,
              bool abort_if_not_found=true) const;
      bool has_input_global(const std::string &globalVarName) const;

      void add_input_field(const stk::io::MeshField &mesh_field);
      void add_input_field(size_t mesh_index, const stk::io::MeshField &mesh_field);

      // Create an exodus mesh database with the specified
      // filename. See STK IO documentation tests for demonstrations of
      // the proper sequence of calls needed to write an exodus database
      // with transient field-data, etc.
      //
      // \param[in] filename The full pathname to the file which will be
      // created and the mesh data written to. If the file already
      // exists, it will be overwritten.
      //
      // \param[in] purpose The purpose of the output mesh.
      // An output of type WRITE_RESTART will write all-but-one 
      //    of the states of a multi-state field,
      // An output of type WRITE_RESULTS will only write 
      //    the newest state of a multi-state field.
      // Other behavioral differences may be added in the future 
      //    (e.g., dealing with adaptivity...)
      // \param[in] type The format of the mesh that will be output.
      // Valid types are "exodus", "catalyst".
      size_t create_output_mesh(const std::string &filename,
				DatabasePurpose purpose,
                                char const* type = "exodus", bool openFileImmediately = true);
      size_t create_output_mesh(const std::string &filename,
                                DatabasePurpose purpose,
                                double time,
                                char const* type = "exodus", bool openFileImmediately = true);
      size_t create_output_mesh(const std::string &filename,
				DatabasePurpose purpose,
				Ioss::PropertyManager &properties,
                                char const* type = "exodus", bool openFileImmediately = true);
      size_t create_output_mesh(const std::string &filename,
				DatabasePurpose purpose,
				Ioss::PropertyManager &properties,
                                double time,
                                char const* type = "exodus", bool openFileImmediately = true);
 
      // Free up memory by removing resouces associated with output files that will no longer be used by the run
      void close_output_mesh(size_t output_file_index);

      // write_output_mesh writes the non-transient portion
      // of the mesh, including the number and type of element blocks,
      // nodesets, and sidesets, and then outputs the mesh bulk data such as the
      // node coordinates, id maps, element connectivity.
      void write_output_mesh(size_t output_file_index);

      void add_field(size_t output_file_index, stk::mesh::FieldBase &field);
      void add_field(size_t output_file_index, stk::mesh::FieldBase &field, const std::string &db_name);
      void add_field(size_t output_file_index, stk::mesh::FieldBase &field, stk::mesh::EntityRank var_type, const std::string &db_name);
      void add_field(size_t output_file_index, stk::mesh::FieldBase &field, stk::mesh::EntityRank var_type, const OutputVariableParams &var);

      void add_attribute_field(size_t output_file_index, stk::mesh::FieldBase &field, const OutputVariableParams &var);

      void add_user_data(size_t output_file_index,
                     const std::vector<std::string> &parts,
                     const std::string &db_name,
                     stk::io::DataLocation loc);
      bool has_global(size_t output_file_index,
                      const std::string &globalVarName) const;
      void add_global_ref(size_t output_file_index,
			  const std::string &variableName,
			  const stk::util::Parameter &param);
      void add_global(size_t output_file_index,
		      const std::string &variableName,
		      const stk::util::Parameter &param);
      template<typename T>
      void add_global(size_t output_file_index,
              const std::string& variableName,
              const T& value,
              stk::util::ParameterType::Type type)
      {
        validate_output_file_index(output_file_index);
        m_outputFiles[output_file_index]->add_global(variableName, value, type);
      }

      void add_global(size_t output_file_index,
		      const std::string &variableName,
		      Ioss::Field::BasicType dataType);
      void add_global(size_t output_file_index,
		      const std::string &variableName,
		      const std::string &type,
		      Ioss::Field::BasicType dataType);
      void add_global(size_t output_file_index,
		      const std::string &variableName,
		      int component_count,
		      Ioss::Field::BasicType dataType);

      // Add a transient step to the database at time 'time'.
      void begin_output_step(size_t output_file_index, double time);
      void end_output_step(size_t output_file_index);

      // Output the data for all defined fields to the database for
      // the step added by "begin_output_step".  End step with a call
      // to "end_output_step"
      int write_defined_output_fields(size_t output_file_index, const stk::mesh::FieldState *state = nullptr) const;

      // Force all output databases to "flush" their data to disk (if possible)
      // Typically called by the application during a planned or unplanned
      // termination.
      void flush_output() const; 

      // Add a transient step to the mesh database at time 'time' and
      // output the data for all defined fields to the database.
      // Performs the same functions as:
      //
      // 	begin_output_step(output_file_index, time);
      // 	write_defined_output_fields(output_file_index);
      // 	end_output_step(output_file_index);
      //
      // NOTE: if there are any global variables defined, they will
      // *not* be output with this function and will be zero-filled on
      // the database. Call the individual functions listed above if
      // you have global variables to output.
      int process_output_request(size_t output_file_index, double time);

      template<typename T>
      void write_global(size_t output_file_index, const std::string &globalVarName,
                        const T &value, stk::util::ParameterType::Type type) const
      {
          validate_output_file_index(output_file_index);
          stk::util::Parameter p;
          p.value = value;
          p.type = type;
          m_outputFiles[output_file_index]->write_global(globalVarName, p);
      }

      void write_global(size_t output_file_index,
                        const std::string& variableName,
                        const stk::util::Parameter& param) const;

      void write_global(size_t output_file_index,
			const std::string &variableName,
      double data) const;
      void write_global(size_t output_file_index,
			const std::string &variableName,
      int data) const;
      void write_global(size_t output_file_index,
			const std::string &variableName,
      std::vector<double>& data) const;
      void write_global(size_t output_file_index,
			const std::string &variableName,
      std::vector<int>& data) const;

      // Add a history or heartbeat output...
      size_t add_heartbeat_output(const std::string &filename,
                                  HeartbeatType db_type,
                                  const Ioss::PropertyManager &properties = Ioss::PropertyManager(),
                                  bool openFileImmediately = true);
  
      void define_heartbeat_global(size_t index,
                                   const std::string &name,
                                   const stk::util::Parameter &param,
                                   int copies = 1,
                                   Ioss::Field::RoleType role = Ioss::Field::REDUCTION);

      void define_heartbeat_global(size_t index,
                                   const std::string &globalVarName,
                                   const stk::util::Parameter &param,
                                   const std::string &storage,
                                   Ioss::Field::BasicType dataType,
                                   int copies = 1,
                                   Ioss::Field::RoleType role = Ioss::Field::REDUCTION);

      void add_heartbeat_global(size_t index,
                                const std::string &name,
                                const stk::util::Parameter &param,
                                int copies = 1,
                                Ioss::Field::RoleType role = Ioss::Field::REDUCTION);
  
      void add_heartbeat_global(size_t index,
                                const std::string &globalVarName,
                                const stk::util::Parameter &param,
                                const std::string &storage,
                                Ioss::Field::BasicType dataType,
                                int copies = 1,
                                Ioss::Field::RoleType role = Ioss::Field::REDUCTION);

      bool has_heartbeat_global(size_t output_file_index,
                                const std::string &globalVarName) const;
      size_t get_heartbeat_global_component_count(size_t output_file_index,
                                                  const std::string &globalVarName) const;
      void process_heartbeat_output(size_t index, int step, double time);

      void process_heartbeat_output_pre_write(size_t index, int step, double time);
      void process_heartbeat_output_write(size_t index, int step, double time);
      void process_heartbeat_output_post_write(size_t index, int step, double time);

      void use_simple_fields();

      bool is_meta_data_null() const;
      bool is_bulk_data_null() const;
      stk::mesh::MetaData &meta_data();
      stk::mesh::BulkData &bulk_data();
      const stk::mesh::MetaData &meta_data() const;
      const stk::mesh::BulkData &bulk_data() const;

      std::shared_ptr<stk::mesh::MetaData> meta_data_ptr() { return m_metaData; }
      std::shared_ptr<stk::mesh::BulkData> bulk_data_ptr() { return m_bulkData; }

      std::shared_ptr<const stk::mesh::MetaData> meta_data_ptr() const { return m_metaData; }
      std::shared_ptr<const stk::mesh::BulkData> bulk_data_ptr() const { return m_bulkData; }

      // Return the coordinate field for this mesh.
      stk::mesh::FieldBase const& get_coordinate_field() const;

      // If there are nodal fields defined on parts with
      // higher-rank. For example, a nodal field on all the nodes of
      // the elements in an element block, then there are two ways to
      // output that field:
      // 1. Create a nodeset of those nodes and output the field on it
      // 2. Output the fields as nodal fields (on all nodes in the model; 
      //        zero-fill if not defined on a node)
      //
      // The first method is the most space-efficient; however, many
      // visualization packages do not handle nodeset fields as well
      // as they handle nodal field.  This function specifies how the
      // user application wants the fields output.  
      bool use_nodeset_for_block_nodes_fields(size_t outputIndex) const;
      void use_nodeset_for_block_nodes_fields(size_t outputIndex,
					     bool flag);
      bool use_nodeset_for_sideset_nodes_fields(size_t outputIndex) const;
      void use_nodeset_for_sideset_nodes_fields(size_t outputIndex,
                                             bool flag);
      bool use_nodeset_for_part_nodes_fields(size_t outputIndex) const;
      void use_nodeset_for_part_nodes_fields(size_t outputIndex,
                                             bool flag);
      bool check_field_existence_when_creating_nodesets(size_t outputIndex) const;
      void check_field_existence_when_creating_nodesets(size_t outputIndex,
                                                        bool flag);
      void use_part_id_for_output(size_t output_file_index, bool flag);
      bool use_part_id_for_output(size_t output_file_index) const;

      void set_throw_on_missing_input_fields(bool flag);
      bool get_throw_on_missing_input_fields() const;

      void set_option_to_not_collapse_sequenced_fields();
      int get_num_time_steps() const;
      double get_max_time() const;
      std::vector<double> get_time_steps() const;
      void set_max_num_steps_before_overwrite(size_t outputFileIndex, int maxNumStepsInFile);

      void set_name_and_version_for_qa_record(size_t outputFileIndex, const std::string &codeName, const std::string &codeVersion);
      void add_qa_records(size_t outputFileIndex, const std::vector<QaRecord> &qaRecords);
      void add_info_records(size_t outputFileIndex, const std::vector<std::string> &infoRecords);
      std::vector<QaRecord> get_qa_records() const;
      std::vector<std::string> get_info_records() const;
      stk::mesh::FieldVector get_ordered_attribute_fields(const stk::mesh::Part *blockPart) const;
      const std::vector<std::vector<int>> & get_attribute_field_ordering_stored_by_part_ordinal() const;
      void set_attribute_field_ordering_stored_by_part_ordinal(const std::vector<std::vector<int>> &ordering);
      void fill_coordinate_frames(std::vector<int>& ids, std::vector<double>& coords, std::vector<char>& tags) const;
      std::string get_output_filename(size_t outputIndex) const;

      bool is_output_index_valid(size_t outputIndex) const;
      bool is_input_index_valid(size_t outputIndex) const;

      Ioss::DatabaseIO *get_input_database(size_t inputIndex) const;
      Ioss::DatabaseIO *get_output_database(size_t outputIndex) const;

      bool set_input_multistate_suffixes(size_t inputIndex, const std::vector<std::string>& multiStateSuffixes);
      bool set_output_multistate_suffixes(size_t outputIndex, const std::vector<std::string>& multiStateSuffixes);

      void set_reference_input_region(size_t outputIndex, const StkMeshIoBroker& inputBroker);

      bool create_named_suffix_field_type(const std::string& type_name, const std::vector<std::string>& suffices) const;

      bool add_field_type_mapping(const std::string& field, const std::string& type) const;

      // Returns 4 or 8 based on several hueristics to determine
      // the integer size required for an output database.

      int check_integer_size_requirements_serial() const;
      int check_integer_size_requirements_parallel() const;
      int check_integer_size_requirements() const;
      void create_surface_to_block_mapping();

      std::vector<stk::mesh::Entity> get_output_entities(size_t outputIndex, const stk::mesh::BulkData& bulk_data, const std::string &name) const;
      //-END
    protected:
      void set_sideset_face_creation_behavior_for_testing(SideSetFaceCreationBehavior behavior)
      {
          m_sidesetFaceCreationBehavior = behavior;
      }

    protected:
      void create_bulk_data();
      void validate_input_file_index(size_t input_file_index) const;
      void stk_mesh_resolve_node_sharing() { bulk_data().resolve_node_sharing(); }
      void stk_mesh_modification_end_after_node_sharing_resolution() { bulk_data().modification_end_after_node_sharing_resolution(); }
    private:
      void create_sideset_observer();
      void create_ioss_region();
      void validate_output_file_index(size_t output_file_index) const;
      void validate_heartbeat_file_index(size_t heartbeat_file_index) const;
      
      void check_for_missing_input_fields(std::vector<stk::io::MeshField> *missingFields);

      void copy_property_manager(const Ioss::PropertyManager &properties);

      Ioss::Property property_get(const std::string &property_name) const;

      // The `m_property_manager` member data contains properties that
      // can be used to set database-specific options in the
      // Ioss::DatabaseIO class.  Examples include compression, name
      // lengths, integer sizes, floating point sizes. By convention,
      // the property name is all uppercase. Some existing properties
      // recognized by the Exodus Ioex::DatabaseIO class are:
      //
      // | Property              | Value
      // |-----------------------|-------------------
      // | COMPRESSION_LEVEL     | In the range [0..9]. A value of 0 indicates no compression
      // | COMPRESSION_SHUFFLE   | (true/false) to enable/disable hdf5's shuffle compression algorithm.
      // | FILE_TYPE             | netcdf4
      // | MAXIMUM_NAME_LENGTH   | Maximum length of names that will be returned/passed via api call.
      // | INTEGER_SIZE_DB       | 4 or 8 indicating byte size of integers stored on the database.
      // | INTEGER_SIZE_API      | 4 or 8 indicating byte size of integers used in api functions.
      // | LOGGING               | (true/false) to enable/disable logging of field input/output
      Ioss::PropertyManager m_propertyManager;

      stk::ParallelMachine m_communicator;
      std::vector<std::string>       m_rankNames; // Optional rank name vector.

      std::shared_ptr<stk::mesh::MeshBuilder> m_meshBuilder;
      std::shared_ptr<stk::mesh::MetaData>  m_metaData;
      std::shared_ptr<stk::mesh::BulkData>  m_bulkData;


      stk::mesh::Selector m_activeSelector;
      stk::mesh::Selector m_subsetSelector;
      std::shared_ptr<stk::mesh::Selector> m_deprecatedSelector;

    protected:
      std::vector<std::vector<int>> attributeFieldOrderingByPartOrdinal;
      std::vector<std::shared_ptr<impl::OutputFile> > m_outputFiles;
    private:
      std::vector<std::shared_ptr<impl::Heartbeat> > m_heartbeat;
    protected:
      std::vector<std::shared_ptr<InputFile> > m_inputFiles;
    private:
      StkMeshIoBroker(const StkMeshIoBroker&); // Do not implement
      StkMeshIoBroker& operator=(const StkMeshIoBroker&); // Do not implement
    void store_attribute_field_ordering();
    void update_sidesets();

    protected:
      size_t m_activeMeshIndex;
      SideSetFaceCreationBehavior m_sidesetFaceCreationBehavior;
      bool m_autoLoadAttributes;
      bool m_autoLoadDistributionFactorPerNodeSet;
      bool m_enableEdgeIO;
      bool m_cacheEntityListForTransientSteps;
      bool m_throwOnMissingInputFields{false};
    };

    inline std::shared_ptr<Ioss::Region> StkMeshIoBroker::get_output_ioss_region(size_t output_file_index) const {
      validate_output_file_index(output_file_index);
      return m_outputFiles[output_file_index]->get_output_ioss_region();
    }

    inline std::shared_ptr<Ioss::Region> StkMeshIoBroker::get_heartbeat_ioss_region(size_t heartbeat_file_index) const {
      validate_heartbeat_file_index(heartbeat_file_index);
      return m_heartbeat[heartbeat_file_index]->get_heartbeat_ioss_region();
    }

    inline void StkMeshIoBroker::begin_define_transient_for_heartbeat(size_t heartbeat_file_index) {
      validate_heartbeat_file_index(heartbeat_file_index);
      return m_heartbeat[heartbeat_file_index]->begin_define_transient();
    }

    inline void StkMeshIoBroker::end_define_transient_for_heartbeat(size_t heartbeat_file_index) {
      validate_heartbeat_file_index(heartbeat_file_index);
      return m_heartbeat[heartbeat_file_index]->end_define_transient();
    }

    inline void StkMeshIoBroker::set_subset_selector(size_t output_file_index,
						     std::shared_ptr<stk::mesh::Selector> my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_subset_selector(my_selector);
    }

    inline void StkMeshIoBroker::set_subset_selector(size_t output_file_index,
						     const stk::mesh::Selector &my_selector) {
      m_subsetSelector = my_selector;
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_subset_selector(std::make_shared<stk::mesh::Selector>(m_subsetSelector));
    }

    inline void StkMeshIoBroker::set_skin_mesh_selector(size_t output_file_index,
						     std::shared_ptr<stk::mesh::Selector> my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_skin_mesh_selector(my_selector);
    }

    inline void StkMeshIoBroker::set_skin_mesh_selector(size_t output_file_index,
						     stk::mesh::Selector &my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_skin_mesh_selector(std::make_shared<stk::mesh::Selector>(my_selector));
    }

    inline void StkMeshIoBroker::set_shared_selector(size_t output_file_index,
                                                     std::shared_ptr<stk::mesh::Selector> my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_shared_selector(my_selector);
    }

    inline void StkMeshIoBroker::set_shared_selector(size_t output_file_index,
                                                     stk::mesh::Selector &my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_shared_selector(std::make_shared<stk::mesh::Selector>(my_selector));
    }

    inline void StkMeshIoBroker::set_output_selector(size_t output_file_index,
                                                     stk::topology::rank_t rank,
                                                     std::shared_ptr<stk::mesh::Selector> my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_output_selector(rank, my_selector);
    }

    inline void StkMeshIoBroker::set_output_selector(size_t output_file_index,
                                                     stk::topology::rank_t rank,
                                                     stk::mesh::Selector &my_selector) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_output_selector(rank, std::make_shared<stk::mesh::Selector>(my_selector));
    }

    inline void StkMeshIoBroker::set_ghosting_filter(size_t output_file_index, bool hasGhosting) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->has_ghosting(hasGhosting);
    }

    inline void StkMeshIoBroker::set_adaptivity_filter(size_t output_file_index, bool hasAdaptivity) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->has_adaptivity(hasAdaptivity);
    }

    inline void StkMeshIoBroker::set_filter_empty_output_entity_blocks(size_t output_file_index, const bool filterEmptyEntityBlocks) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_filter_empty_entity_blocks(filterEmptyEntityBlocks);
    }

    inline void StkMeshIoBroker::set_filter_empty_output_assembly_entity_blocks(size_t output_file_index, const bool filterEmptyAssemblyEntityBlocks) {
      validate_output_file_index(output_file_index);
      m_outputFiles[output_file_index]->set_filter_empty_assembly_entity_blocks(filterEmptyAssemblyEntityBlocks);
    }

    inline stk::mesh::Selector StkMeshIoBroker::get_active_selector() const {
      return m_activeSelector;
    }

    inline void StkMeshIoBroker::set_active_selector(stk::mesh::Selector my_selector) {
      m_activeSelector = my_selector;
    }

    inline void StkMeshIoBroker::set_bulk_data(stk::mesh::BulkData &arg_bulk_data)
    {
      set_bulk_data(std::shared_ptr<stk::mesh::BulkData>(&arg_bulk_data, [](auto pointerWeWontDelete){}));
    }

    inline void StkMeshIoBroker::replace_bulk_data(stk::mesh::BulkData &arg_bulk_data)
    {
      replace_bulk_data(std::shared_ptr<stk::mesh::BulkData>(&arg_bulk_data, [](auto pointerWeWontDelete){}));
    }

    inline void StkMeshIoBroker::define_heartbeat_global(size_t index,
                                                         const std::string &name,
                                                         const stk::util::Parameter &param,
                                                         int copies,
                                                         Ioss::Field::RoleType role)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->define_global_ref(name, param, copies, role);
    }

    inline void StkMeshIoBroker::define_heartbeat_global(size_t index,
                                                         const std::string &globalVarName,
                                                         const stk::util::Parameter &param,
                                                         const std::string &storage,
                                                         Ioss::Field::BasicType dataType,
                                                         int copies,
                                                         Ioss::Field::RoleType role)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->define_global_ref(globalVarName, param, storage, dataType, copies, role);
    }

    inline void StkMeshIoBroker::add_heartbeat_global(size_t index,
						      const std::string &name,
						      const stk::util::Parameter &param,
						      int copies,
						      Ioss::Field::RoleType role)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->add_global_ref(name, param, copies, role);
    }
  
    inline void StkMeshIoBroker::add_heartbeat_global(size_t index,
                                                      const std::string &globalVarName,
                                                      const stk::util::Parameter &param,
                                                      const std::string &storage,
                                                      Ioss::Field::BasicType dataType,
                                                      int copies,
                                                      Ioss::Field::RoleType role)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->add_global_ref(globalVarName, param, storage, dataType, copies, role);
    }

    inline void StkMeshIoBroker::process_heartbeat_output(size_t index, int step, double time)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->process_output(step, time);
    }
  
    inline void StkMeshIoBroker::process_heartbeat_output_pre_write(size_t index, int step, double time)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->process_output_pre_write(step, time);
    }

    inline void StkMeshIoBroker::process_heartbeat_output_write(size_t index, int step, double time)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->process_output_write(step, time);
    }

    inline void StkMeshIoBroker::process_heartbeat_output_post_write(size_t index, int step, double time)
    {
      STKIORequire(index < m_heartbeat.size());
      m_heartbeat[index]->process_output_post_write(step, time);
    }

    inline bool StkMeshIoBroker::is_meta_data_null() const
    {
      return m_metaData == nullptr;
    }

    inline bool StkMeshIoBroker::is_bulk_data_null() const
    {
      return m_bulkData == nullptr;
    }

    inline stk::mesh::MetaData &StkMeshIoBroker::meta_data()
    {
      STK_ThrowAssert(!is_meta_data_null());
      return *m_metaData;
    }

    inline stk::mesh::BulkData &StkMeshIoBroker::bulk_data()
    {
      STK_ThrowAssert(!is_bulk_data_null());
      return *m_bulkData;
    }

    inline const stk::mesh::MetaData &StkMeshIoBroker::meta_data() const
    {
      STK_ThrowAssert(!is_meta_data_null());
      return *m_metaData;
    }

    inline const stk::mesh::BulkData &StkMeshIoBroker::bulk_data() const
    {
      STK_ThrowAssert(!is_bulk_data_null());
      return *m_bulkData;
    }
  }
}
#endif
