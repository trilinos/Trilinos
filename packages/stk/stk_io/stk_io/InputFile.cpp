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

#include <stk_io/InputFile.hpp>
#include <math.h>                       // for fmod
#include <stddef.h>                     // for NULL, size_t
#include <algorithm>                    // for sort, swap
#include <limits>                       // for numeric_limits
#include <ostream>                      // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_io/DbStepTimeInterval.hpp>  // for DBStepTimeInterval
#include <stk_io/IossBridge.hpp>        // for get_field_role, etc
#include <stk_io/MeshField.hpp>         // for MeshField, etc
#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase, etc
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowErrorMsgIf, etc
#include <utility>                      // for pair
#include "Ioss_DBUsage.h"               // for DatabaseUsage, etc
#include "Ioss_DatabaseIO.h"            // for DatabaseIO
#include "Ioss_EntityType.h"            // for EntityType::SIDESET
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_GroupingEntity.h"        // for GroupingEntity
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include "Teuchos_Ptr.hpp"              // for Ptr::get
#include "Teuchos_PtrDecl.hpp"          // for Ptr
#include "Teuchos_RCP.hpp"              // for is_null, RCP::operator->, etc
#include "stk_io/DatabasePurpose.hpp"
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc

#include "SidesetTranslator.hpp"
#include "StkIoUtils.hpp"


namespace {
  bool meshFieldSort(const stk::io::MeshField& f1, const stk::io::MeshField &f2) {
    return f1.field()->mesh_meta_data_ordinal() < f2.field()->mesh_meta_data_ordinal();
  }
}

namespace stk {
  namespace io {

    InputFile::InputFile(std::string mesh_filename,
			 MPI_Comm communicator, 
			 const std::string &mesh_type,
			 DatabasePurpose purpose,
			 Ioss::PropertyManager &properties)
      : m_db_purpose(purpose), m_region(NULL),
	m_startupTime(0.0),
	m_periodLength(0.0),
	m_scaleTime(1.0),
	m_offsetTime(0.0),
	m_startTime(-std::numeric_limits<double>::max()),
	m_stopTime(std::numeric_limits<double>::max()),
	m_periodType(CYCLIC),
	m_fieldsInitialized(false)
    {
      Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL;
      if (m_db_purpose == stk::io::READ_RESTART)
        db_usage = Ioss::READ_RESTART;
        
      stk::util::filename_substitution(mesh_filename);
      m_database = Teuchos::rcp(Ioss::IOFactory::create(mesh_type, mesh_filename,
							db_usage, communicator,
							properties));
      ThrowErrorMsgIf(Teuchos::is_null(m_database) || !m_database->ok(true), 
		      "ERROR: Could not open database '" << mesh_filename
		      << "' of type '" << mesh_type << "'");
    }


    InputFile::InputFile(Teuchos::RCP<Ioss::Region> ioss_input_region)
      : m_database(ioss_input_region->get_database()), m_region(ioss_input_region),
	m_startupTime(0.0),
	m_periodLength(0.0),
	m_scaleTime(1.0),
	m_offsetTime(0.0),
	m_startTime(-std::numeric_limits<double>::max()),
	m_stopTime(std::numeric_limits<double>::max()),
	m_periodType(CYCLIC),
	m_fieldsInitialized(false)
    {
      ThrowErrorMsgIf(Teuchos::is_null(m_database) || !m_database->ok(true), 
		      "ERROR: Invalid Ioss region detected in add_mesh_database");

      Ioss::DatabaseUsage db_usage = m_database->usage();
      if (db_usage == Ioss::READ_RESTART) {
	m_db_purpose = stk::io::READ_RESTART;
      }
      else if (db_usage == Ioss::READ_MODEL) {
	m_db_purpose = stk::io::READ_MESH;
      }
      else {
        std::ostringstream msg;
        msg << "ERROR: Unrecognized database usage for Ioss region named "
	    << ioss_input_region->name()
	    << ". Must be READ_RESTART or READ_MODEL";
        throw std::runtime_error( msg.str() );
      }

      ThrowErrorMsgIf(m_region->mesh_type() != Ioss::MeshType::UNSTRUCTURED,
		      "Mesh type is '" << m_region->mesh_type_string() << "' which is not supported. "
		      "Only 'Unstructured' mesh is currently supported.");

      m_database.release(); // The m_region will delete the m_database pointer.
    }

    void InputFile::create_ioss_region()
    {
      // If the m_region is null, try to create it from
      // the m_input_database. If that is null, throw an error.
      if (Teuchos::is_null(m_region)) {
        ThrowErrorMsgIf(Teuchos::is_null(m_database),
                        "There is no input mesh database associated with this StkMeshIoBroker. Please call open_mesh_database() first.");
        // The Ioss::Region takes control of the m_input_database pointer, so we need to make sure the
        // RCP doesn't retain ownership...
        m_region = Teuchos::rcp(new Ioss::Region(m_database.release().get(), "input_model"));

	ThrowErrorMsgIf(m_region->mesh_type() != Ioss::MeshType::UNSTRUCTURED,
			"Mesh type is '" << m_region->mesh_type_string() << "' which is not supported. "
			"Only 'Unstructured' mesh is currently supported.");
      }
    }

    void InputFile::add_input_field(const stk::io::MeshField &mesh_field)
    {
      bool fieldAlreadyExists=false;
      for (size_t i=0; i <m_fields.size(); i++) {
        if (mesh_field == m_fields[i]) {
          fieldAlreadyExists = true;
          break;
        }
      }

      if (!fieldAlreadyExists) {
        m_fields.push_back(mesh_field);
        stk::io::set_field_role(*mesh_field.field(), Ioss::Field::TRANSIENT);
	m_fieldsInitialized = false;
      }
    }

    InputFile& InputFile::set_offset_time(double offset_time)
    {
      m_offsetTime = offset_time;
      return *this;
    }

    InputFile& InputFile::set_scale_time(double scale_time)
    {
      m_scaleTime = scale_time;
      return *this;
    }

    InputFile& InputFile::set_periodic_time(double period_length, double startup_time,
					    PeriodType ptype)
    {
      m_periodLength = period_length;
      m_startupTime = startup_time;
      m_periodType = ptype;
      return *this;
    }

    InputFile& InputFile::set_start_time(double start_time)
    {
      m_startTime = start_time;
      return *this;
    }

    InputFile& InputFile::set_stop_time(double stop_time)
    {
      m_stopTime = stop_time;
      return *this;
    }

    void InputFile::get_global_variable_names(std::vector<std::string> &names)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_region),
		       "Attempt to read global variables before restart initialized.");
      m_region->field_describe(Ioss::Field::TRANSIENT, &names);
    }

    FieldNameToPartVector InputFile::get_var_names(Ioss::EntityType type, stk::mesh::MetaData& meta)
    {
        return stk::io::get_var_names(*m_region.get(), type, meta);
    }

    void InputFile::add_all_mesh_fields_as_input_fields(stk::mesh::MetaData &meta, MeshField::TimeMatchOption tmo)
    {
      stk::io::define_input_fields(*m_region.get(),  meta);

      // Iterate all fields and set them as input fields...
      const stk::mesh::FieldVector &fields = meta.get_fields();
      for (size_t i=0; i < fields.size(); i++) {
        const Ioss::Field::RoleType* role = stk::io::get_field_role(*fields[i]);
        if ( role && *role == Ioss::Field::TRANSIENT ) {
	  add_input_field(MeshField(fields[i], fields[i]->name(), tmo));
	}
      }
    }

    bool InputFile::read_input_field(stk::io::MeshField &mf, stk::mesh::BulkData &bulk)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_region),
		       "ERROR: There is no Input mesh/restart region associated with this Mesh Data.");
      Ioss::Region *region = m_region.get();
      double time = mf.m_timeToRead;
      if (time < m_startTime || time > m_stopTime)
	return false;
      
      // Map analysis time to database time using offset, periodic, ...
      // See details in header file.
      double db_time = map_analysis_to_db_time(time);

      // Get struct containing interval of database time(s) containing 'time'
      DBStepTimeInterval sti(region, db_time);

      ThrowErrorMsgIf(!sti.exists_before && !sti.exists_after,
		      "ERROR: Input database '" << region->get_database()->get_filename()
		      << "' has no transient data.");

      const stk::mesh::FieldBase *f = mf.field();
      std::vector<const stk::mesh::Part*>::iterator P = mf.m_subsetParts.begin();
      while (P != mf.m_subsetParts.end()) {
	// Find the Ioss::GroupingEntity corresponding to this part...
	mf.set_inactive();
	const stk::mesh::Part *part = *P; ++P;
	stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
	if (f->entity_rank() == rank) {
	  Ioss::GroupingEntity *io_entity = region->get_entity(part->name());
	  ThrowErrorMsgIf( io_entity == NULL,
			   "ERROR: For field '" << mf.field()->name()
			   << "' Could not find database entity corresponding to the part named '"
			   << part->name() << "'.");
	  build_field_part_associations(mf, *part, rank, io_entity);
	} 

	// If rank is != NODE_RANK, then see if field is defined on the nodes of this part
	if (rank != stk::topology::NODE_RANK && f->entity_rank() == stk::topology::NODE_RANK) {
	  Ioss::GroupingEntity *node_entity = NULL;
	  std::string nodes_name = part->name() + "_nodes";
	  node_entity = region->get_entity(nodes_name);
	  if (node_entity == NULL) {
	    node_entity = region->get_entity("nodeblock_1");
	  }
	  if (node_entity != NULL) {
	    build_field_part_associations(mf, *part, stk::topology::NODE_RANK, node_entity);
	  }
	}

	if (mf.is_active()) {
	  mf.restore_field_data(bulk, sti);
	}
      }
      return true;
    }

    double InputFile::read_defined_input_fields(int step,
						std::vector<stk::io::MeshField> *missingFields,
						stk::mesh::BulkData &bulk)
    {
      ThrowErrorMsgIf(step <= 0, 
		      "ERROR: Invalid step (" << step << ") requested. Value must be greater than zero.");

      ThrowErrorMsgIf (Teuchos::is_null(m_region),
                       "There is no Input mesh/restart region associated with this Mesh Data.");

      Ioss::Region *region = m_region.get();

      int step_count = region->get_property("state_count").get_int();

      ThrowErrorMsgIf(step_count == 0, 
		      "ERROR: Input database '" << region->get_database()->get_filename()
		      << "' has no transient data.");

      ThrowErrorMsgIf(step > step_count,
		      "ERROR: Input database '" << region->get_database()->get_filename()
		      << "'. Step " << step << " was specified, but database only has "
		      << step_count << " steps.");

      double state_time = region->get_state_time(step);
      return read_defined_input_fields(state_time, missingFields, bulk);
    }

    bool InputFile::build_field_part_associations(stk::io::MeshField &mesh_field,
						  const stk::mesh::Part &part,
						  const stk::mesh::EntityRank rank,
						  Ioss::GroupingEntity *io_entity,
						  std::map<stk::mesh::FieldBase *, const stk::io::MeshField *> *missing_fields_collector)
    {
        bool field_is_missing = false;
        stk::mesh::FieldBase *f = mesh_field.field();
        // Only add TRANSIENT Fields -- check role; if not present assume transient...
        const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
        if (role == NULL || *role == Ioss::Field::TRANSIENT) {
            if (stk::io::is_field_on_part(f, rank, part)) {
                const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, rank, part);
                std::pair<std::string, Ioss::Field::BasicType> field_type;
                stk::io::get_io_field_type(f, res, &field_type);
                if (field_type.second != Ioss::Field::INVALID) {

                    const std::string &db_name = mesh_field.db_name();
                    unsigned num_states = f->number_of_states();
                    std::vector<stk::mesh::FieldState> missing_states;
                    if (num_states > 1 && !all_field_states_exist_on_io_entity(db_name, f, io_entity, missing_states)) {
                        field_is_missing = true;
                        if (missing_fields_collector) {
                            for (stk::mesh::FieldState missing_state : missing_states)
                                (*missing_fields_collector)[f->field_state(missing_state)] = &mesh_field;
                        }
                    }

                    bool field_exists = io_entity->field_exists(db_name);
                    if (!field_exists) {
                        field_is_missing = true;
                        if (missing_fields_collector) {
                            (*missing_fields_collector)[f] = &mesh_field;
                        }
                    }

                    // See if field with that name exists on io_entity...
                    if (field_exists) {
                        mesh_field.add_part(rank, part, io_entity);
                        mesh_field.set_single_state((m_db_purpose == stk::io::READ_RESTART) ? false : true);
                        mesh_field.set_active();
                    }
                }
            }
        }
        return field_is_missing;
    }

    void InputFile::build_field_part_associations(stk::mesh::BulkData &bulk, std::vector<stk::io::MeshField> *missingFields)
    {
        std::map<stk::mesh::FieldBase *, const stk::io::MeshField *> missing_fields_collector;
        std::map<stk::mesh::FieldBase *, const stk::io::MeshField *> *missing_fields_collector_ptr =
                (missingFields ? &missing_fields_collector : 0);

        // Each input field will have a list of the Parts that the field exists on...
        // Create this list.
        Ioss::Region *region = m_region.get();
        size_t num_missing_fields = 0;
        // First handle any fields that are subsetted (restricted to a specified list of parts)
        {
            std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
            while (I != m_fields.end()) {
                const stk::mesh::FieldBase *f = (*I).field();
                std::vector<const stk::mesh::Part*>::iterator P = (*I).m_subsetParts.begin();
                while (P != (*I).m_subsetParts.end()) {
                    // Find the Ioss::GroupingEntity corresponding to this part...
                    const stk::mesh::Part *part = *P; ++P;
                    stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
                    bool field_is_missing = false;
                    if (f->entity_rank() == rank) {
                        Ioss::GroupingEntity *io_entity = region->get_entity(part->name());
                        ThrowErrorMsgIf( io_entity == NULL,
                                "ERROR: For field '" << (*I).field()->name()
                                << "' Could not find database entity corresponding to the part named '"
                                << part->name() << "'.");
                        field_is_missing = build_field_part_associations(*I, *part, rank, io_entity, missing_fields_collector_ptr);
                    }

                    // If rank is != NODE_RANK, then see if field is defined on the nodes of this part
                    if (rank != stk::topology::NODE_RANK && f->entity_rank() == stk::topology::NODE_RANK) {
                        Ioss::GroupingEntity *node_entity = NULL;
                        std::string nodes_name = part->name() + "_nodes";
                        node_entity = region->get_entity(nodes_name);
                        if (node_entity == NULL) {
                            node_entity = region->get_entity("nodeblock_1");
                        }
                        if (node_entity != NULL) {
                            field_is_missing = build_field_part_associations(*I, *part, stk::topology::NODE_RANK, node_entity,
                                                                             missing_fields_collector_ptr);
                        }
                    }

                    if (field_is_missing) {
                        ++num_missing_fields;
                    }
                }
                ++I;
            }
        }

        // Now handle the non-subsetted fields...

        // Check universal_part() NODE_RANK first...
        const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);
        {
            std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
            while (I != m_fields.end()) {
                if ((*I).m_subsetParts.empty()) {
                    const stk::mesh::FieldBase *f = (*I).field();
                    if (f->entity_rank() == stk::topology::NODE_RANK) {
                        bool field_is_missing = build_field_part_associations(*I, meta.universal_part(), stk::topology::NODE_RANK,
                                                                              region->get_node_blocks()[0], missing_fields_collector_ptr);
                        if (field_is_missing) {
                            ++num_missing_fields;
                        }
                    }
                }
                ++I;
            }
        }

        // Now handle all non-nodeblock parts...
        const stk::mesh::PartVector &all_parts = meta.get_parts();
        for ( stk::mesh::PartVector::const_iterator
                ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

            stk::mesh::Part * const part = *ip;

            // Check whether this part is an input part...
            if (stk::io::is_part_io_part(*part)) {
                stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
                // Get Ioss::GroupingEntity corresponding to this part...
                Ioss::GroupingEntity *entity = region->get_entity(part->name());
                if (entity != NULL && !m_fields.empty() && entity->type() != Ioss::SIDESET) {
                    std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
                    while (I != m_fields.end()) {
                        if ((*I).m_subsetParts.empty()) {
                            const stk::mesh::FieldBase *f = (*I).field();
                            bool field_is_missing = false;
                            if (f->entity_rank() == rank) {
                                field_is_missing = build_field_part_associations(*I, *part, rank, entity, missing_fields_collector_ptr);
                            }

                            // If rank is != NODE_RANK, then see if field is defined on the nodes of this part
                            if (rank != stk::topology::NODE_RANK && f->entity_rank() == stk::topology::NODE_RANK) {
                                Ioss::GroupingEntity *node_entity = NULL;
                                std::string nodes_name = part->name() + "_nodes";
                                node_entity = region->get_entity(nodes_name);
                                if (node_entity == NULL) {
                                    node_entity = region->get_entity("nodeblock_1");
                                }
                                if (node_entity != NULL) {
                                    field_is_missing = build_field_part_associations(*I, *part, stk::topology::NODE_RANK, node_entity,
                                                                                     missing_fields_collector_ptr);
                                }
                            }

                            if (field_is_missing) {
                                ++num_missing_fields;
                            }
                        }
                        ++I;
                    }
                }
            }
        }

        if (num_missing_fields > 0 && missingFields==NULL) {
            std::ostringstream msg;
            msg << "ERROR: Input field processing could not find " << num_missing_fields << " fields.\n";
            throw std::runtime_error( msg.str() );
        }

        if (missingFields)
        {
            std::vector<stk::io::MeshField> discoveredMissingFields;
            for (auto missingStatedFieldIter : missing_fields_collector)
            {
                discoveredMissingFields.push_back(stk::io::MeshField(missingStatedFieldIter.first,
                                                                     missingStatedFieldIter.second->db_name()));
            }
            std::sort(discoveredMissingFields.begin(), discoveredMissingFields.end(),
                      [](const stk::io::MeshField &a, const stk::io::MeshField &b) {
                            return (a.db_name() < b.db_name())
                                    || ((a.db_name() == b.db_name()) && (a.field()->name() < b.field()->name())); });
            missingFields->insert(missingFields->end(), discoveredMissingFields.begin(), discoveredMissingFields.end());
        }
    }

    double InputFile::map_analysis_to_db_time(double time) const
    {
      double db_time = time;
      if (time > m_startupTime && m_periodLength > 0.0) {
	if (m_periodType == CYCLIC) {
	  db_time = m_startupTime + fmod(time-m_startupTime, m_periodLength);
	} else if (m_periodType == REVERSING) {
	  double pmod = fmod((time-m_startupTime), 2.0*m_periodLength);
	  if (pmod <= m_periodLength) {
	    db_time = m_startupTime + pmod;
	  }
	  else {
	    db_time = m_startupTime + 2.0 * m_periodLength - pmod;
	  }
	}
      }
      db_time *= m_scaleTime;
      db_time += m_offsetTime;
      return db_time;
    }

    double InputFile::read_defined_input_fields(double time,
						std::vector<stk::io::MeshField> *missingFields,
						stk::mesh::BulkData &bulk)
    {
      // Sort fields to ensure they are iterated in the same order on all processors.
      std::sort(m_fields.begin(), m_fields.end(), meshFieldSort);

      bool ignore_missing_fields = (missingFields != NULL);

      if (!m_fieldsInitialized) {
	std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
	while (I != m_fields.end()) {
	  (*I).set_inactive(); ++I;
	}

	build_field_part_associations(bulk, missingFields);
	  
	m_fieldsInitialized = true;
      }

      if (time < m_startTime || time > m_stopTime)
	return 0.0;

      // Map analysis time to database time using offset, periodic, ...
      // See details in header file.
      double db_time = map_analysis_to_db_time(time);
	
      ThrowErrorMsgIf (Teuchos::is_null(m_region),
		       "ERROR: There is no Input mesh/restart region associated with this Mesh Data.");

      Ioss::Region *region = m_region.get();

      // Get struct containing interval of database time(s) containing 'time'
      DBStepTimeInterval sti(region, db_time);

      ThrowErrorMsgIf(!sti.exists_before && !sti.exists_after,
		      "ERROR: Input database '" << region->get_database()->get_filename()
		      << "' has no transient data.");

      std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
      double time_read = -1.0;
      while (I != m_fields.end()) {
	// NOTE: If the fields being restored have different settings, the time
	// value can be different for each field and this will return the value
	// of the last field.  For example, if one field is CLOSEST, one is SPECFIED,
	// and one is TIME_INTERPOLATION, then the time value to return is
	// ambiguous.  Also an issue if some of the fields are inactive.
	double time_t = (*I).restore_field_data(bulk, sti, ignore_missing_fields);
	if ((*I).is_active()) {
	  time_read = time_t > time_read ? time_t : time_read;
	}
	++I;
      }

      size_t step = sti.get_closest_step();
      int current_step = region->get_current_state();
      if (current_step != -1 && current_step != static_cast<int>(step))
	region->end_state(current_step);
      if (current_step != static_cast<int>(step))
	region->begin_state(step);

      return time_read;
    }


    double InputFile::read_defined_input_fields_at_step(int step,
                                                std::vector<stk::io::MeshField> *missingFields,
                                                stk::mesh::BulkData &bulk)
    {
        ThrowErrorMsgIf(step <= 0,
                        "ERROR: Invalid step (" << step << ") requested. Value must be greater than zero.");

        ThrowErrorMsgIf (Teuchos::is_null(m_region),
                         "There is no Input mesh/restart region associated with this Mesh Data.");

        Ioss::Region *region = m_region.get();

        int step_count = region->get_property("state_count").get_int();

        ThrowErrorMsgIf(step_count == 0,
                        "ERROR: Input database '" << region->get_database()->get_filename()
                        << "' has no transient data.");

        ThrowErrorMsgIf(step > step_count,
                        "ERROR: Input database '" << region->get_database()->get_filename()
                        << "'. Step " << step << " was specified, but database only has "
                        << step_count << " steps.");

        // Sort fields to ensure they are iterated in the same order on all processors.
        std::sort(m_fields.begin(), m_fields.end(), meshFieldSort);

        bool ignore_missing_fields = (missingFields != NULL);

        if (!m_fieldsInitialized) {
            std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
            while (I != m_fields.end()) {
                (*I).set_inactive(); ++I;
            }

            build_field_part_associations(bulk, missingFields);

            m_fieldsInitialized = true;
        }

        double time  = region->get_state_time(step);
        if (time < m_startTime || time > m_stopTime)
            return 0.0;

        std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
        double time_read = -1.0;
        while (I != m_fields.end()) {
            // NOTE: If the fields being restored have different settings, the time
            // value can be different for each field and this will return the value
            // of the last field.  For example, if one field is CLOSEST, one is SPECFIED,
            // and one is TIME_INTERPOLATION, then the time value to return is
            // ambiguous.  Also an issue if some of the fields are inactive.
            double time_t = (*I).restore_field_data_at_step(region, bulk, step, ignore_missing_fields);
            if ((*I).is_active()) {
                time_read = time_t > time_read ? time_t : time_read;
            }
            ++I;
        }

        int current_step = region->get_current_state();
        if (current_step != -1 && current_step != static_cast<int>(step))
            region->end_state(current_step);
        if (current_step != static_cast<int>(step))
            region->begin_state(step);

        return time_read;
    }

  }
}
