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
#include <stk_io/InputFile.hpp>
#include <stk_io/InputQuery.hpp>
#include <exception>                    // for exception
#include <algorithm>                           // for copy, sort, max, find
#include <cmath>                               // for fmod
#include <cstddef>                             // for size_t
#include <iostream>                            // for operator<<, basic_ostream
#include <limits>                              // for numeric_limits
#include <stdexcept>                           // for runtime_error
#include <stk_io/DbStepTimeInterval.hpp>       // for DBStepTimeInterval
#include <stk_io/IossBridge.hpp>               // for is_part_io_part, all_f...
#include <stk_io/MeshField.hpp>                // for MeshField, MeshField::...
#include <stk_mesh/base/FieldBase.hpp>         // for FieldBase, FieldBase::...
#include <stk_mesh/base/FindRestriction.hpp>   // for find_restriction
#include <stk_mesh/base/MetaData.hpp>          // for MetaData
#include <stk_util/environment/FileUtils.hpp>  // for filename_substitution
#include <stk_util/util/ReportHandler.hpp>     // for ThrowErrorMsgIf, Throw...
#include <utility>                             // for move, pair
#include "Ioss_DBUsage.h"                      // for DatabaseUsage, READ_MODEL
#include "Ioss_DatabaseIO.h"                   // for DatabaseIO
#include "Ioss_EntityType.h"                   // for SIDESET, EntityType
#include "Ioss_Field.h"                        // for Field, Field::TRANSIENT
#include "Ioss_GroupingEntity.h"               // for GroupingEntity
#include "Ioss_IOFactory.h"                    // for IOFactory
#include "Ioss_MeshType.h"                     // for MeshType, MeshType::UN...
#include "Ioss_NodeBlock.h"                    // for NodeBlock
#include "Ioss_NodeSet.h"                      // for NodeSet
#include "Ioss_Property.h"                     // for Property
#include "Ioss_Region.h"                       // for Region, NodeBlockConta...
#include "Ioss_SideBlock.h"                    // for SideBlock
#include "Ioss_SideSet.h"                      // for SideSet
#include "StkIoUtils.hpp"                      // for part_primary_entity_rank
#include "stk_io/DatabasePurpose.hpp"          // for READ_RESTART, Database...
#include "stk_mesh/base/BulkData.hpp"          // for BulkData
#include "stk_mesh/base/FieldState.hpp"        // for FieldState
#include "stk_mesh/base/Part.hpp"              // for Part
#include "stk_mesh/base/Types.hpp"             // for PartVector, EntityRank
#include "stk_topology/topology.hpp"           // for topology, topology::NO...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################



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
    : m_db_purpose(purpose), m_region(nullptr),
	    m_startupTime(0.0),
	    m_periodLength(0.0),
	    m_scaleTime(1.0),
	    m_offsetTime(0.0),
	    m_startTime(-std::numeric_limits<double>::max()),
	    m_stopTime(std::numeric_limits<double>::max()),
      m_maxTime(-std::numeric_limits<double>::max()),
	    m_periodType(CYCLIC),
	    m_fieldsInitialized(false),
	    m_haveCachedEntityList(false),
	    m_multiStateSuffixes(nullptr)
  {
    Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL;
    if (m_db_purpose == stk::io::READ_RESTART)
      db_usage = Ioss::READ_RESTART;

    stk::util::filename_substitution(mesh_filename);
    m_database = std::shared_ptr<Ioss::DatabaseIO>(Ioss::IOFactory::create(mesh_type, mesh_filename,
 	    				                                                             db_usage, communicator,
      	    				                                                       properties), [](auto /*pointerWeWontDelete*/){});

    if (m_database.get() == nullptr || !m_database->ok(true)) {
      delete m_database.get();
      STK_ThrowErrorMsg("ERROR: Could not open database '" << mesh_filename
                        << "' of type '" << mesh_type << "'");
    }
  }

  InputFile::InputFile(std::shared_ptr<Ioss::Region> ioss_input_region)
    : m_database(ioss_input_region->get_database(), [](auto /*pointerWeWontDelete*/){}), m_region(ioss_input_region),
	    m_startupTime(0.0),
	    m_periodLength(0.0),
	    m_scaleTime(1.0),
	    m_offsetTime(0.0),
	    m_startTime(-std::numeric_limits<double>::max()),
	    m_stopTime(std::numeric_limits<double>::max()),
      m_maxTime(-std::numeric_limits<double>::max()),
	    m_periodType(CYCLIC),
	    m_fieldsInitialized(false),
	    m_haveCachedEntityList(false),
      m_multiStateSuffixes(nullptr)
  {
    STK_ThrowErrorMsgIf(m_database == nullptr || !m_database->ok(true), 
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

    STK_ThrowErrorMsgIf(m_region->mesh_type() != Ioss::MeshType::UNSTRUCTURED,
                        "Mesh type is '" << m_region->mesh_type_string() << "' which is not supported. "
                        "Only 'Unstructured' mesh is currently supported.");

  }

    void InputFile::create_ioss_region()
    {
      // If the m_region is null, try to create it from
      // the m_input_database. If that is null, throw an error.
      if (m_region == nullptr) {
        STK_ThrowErrorMsgIf(m_database.get() == nullptr,
                        "There is no input mesh database associated with this StkMeshIoBroker. Please call open_mesh_database() first.");
        // The Ioss::Region takes control of the m_input_database pointer, so we need to make sure the
        // RCP doesn't retain ownership...
        Ioss::Region *region = nullptr;
        try {
          region = new Ioss::Region(m_database.get(), "input_model");
        } catch (...) {
          m_database.reset();
          throw;
        }
        m_region = std::shared_ptr<Ioss::Region>(region);

        STK_ThrowErrorMsgIf(m_region->mesh_type() != Ioss::MeshType::UNSTRUCTURED,
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
      STK_ThrowErrorMsgIf (m_region.get() == nullptr,
		       "Attempt to read global variables before restart initialized.");
      m_region->field_describe(Ioss::Field::REDUCTION, &names);
    }

    FieldNameToPartVector InputFile::get_var_names(Ioss::EntityType type, const stk::mesh::MetaData& meta)
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
      STK_ThrowErrorMsgIf (m_region == nullptr,
                       "ERROR: There is no Input mesh/restart region associated with this Mesh Data.");

      mf.m_fieldRestored = false;
      mf.m_timeRestored = std::numeric_limits<double>::max();

      Ioss::Region *region = m_region.get();
      double time = mf.m_timeToRead;
      if (time < m_startTime || time > m_stopTime)
        return false;

      // Map analysis time to database time using offset, periodic, ...
      // See details in header file.
      double db_time = map_analysis_to_db_time(time);

      // Get struct containing interval of database time(s) containing 'time'
      DBStepTimeInterval sti(region, db_time);

      STK_ThrowErrorMsgIf(!sti.exists_before && !sti.exists_after,
                      "ERROR: Input database '" << region->get_database()->get_filename()
                      << "' has no transient data.");

      InputQuery iq(*region, bulk.mesh_meta_data(), m_db_purpose, m_multiStateSuffixes);

      const stk::mesh::FieldBase *f = mf.field();

      for (const stk::mesh::Part* part : mf.m_subsetParts) {
        // Find the Ioss::GroupingEntity corresponding to this part...
        mf.set_inactive();

        iq.build_field_part_associations_for_part(mf, part);

        if (mf.is_active()) {
          mf.restore_field_data(bulk, sti, false, m_multiStateSuffixes);
        }
      }

      if(mf.m_subsetParts.empty()) {
        mf.set_inactive();
        // Now handle the non-subsetted fields...

        // Check universal_part() NODE_RANK first...
        const stk::mesh::MetaData &meta = bulk.mesh_meta_data();
        {
          if (f->entity_rank() == stk::topology::NODE_RANK) {
            Ioss::NodeBlock* nb = region->get_node_blocks()[0];
            iq.build_field_part_associations(mf, meta.universal_part(), stk::topology::NODE_RANK, nb);
            iq.process_fields_for_grouping_entity(mf, meta.universal_part(), nb);
          }
        }

        // Now handle all non-nodeblock parts...
        for ( const stk::mesh::Part * part : meta.get_parts() ) {
          // Check whether this part is an input part...
          if (stk::io::is_part_io_part(*part)) {
            iq.build_field_part_associations_for_part(mf, part);
          }
        }

        if (mf.is_active()) {
          mf.restore_field_data(bulk, sti, false, m_multiStateSuffixes);
        }
      }

      return mf.is_active();
    }

    double InputFile::read_defined_input_fields(int step,
						std::vector<stk::io::MeshField> *missingFields,
						stk::mesh::BulkData &bulk)
    {
      STK_ThrowErrorMsgIf(step <= 0, 
		      "ERROR: Invalid step (" << step << ") requested. Value must be greater than zero.");

      STK_ThrowErrorMsgIf (m_region.get() == nullptr,
                       "There is no Input mesh/restart region associated with this Mesh Data.");

      Ioss::Region *region = m_region.get();

      int step_count = region->get_property("state_count").get_int();

      STK_ThrowErrorMsgIf(step_count == 0, 
		      "ERROR: Input database '" << region->get_database()->get_filename()
		      << "' has no transient data.");

      STK_ThrowErrorMsgIf(step > step_count,
		      "ERROR: Input database '" << region->get_database()->get_filename()
		      << "'. Step " << step << " was specified, but database only has "
		      << step_count << " steps.");

      double state_time = region->get_state_time(step);
      return read_defined_input_fields(state_time, missingFields, bulk);
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

      bool ignore_missing_fields = (missingFields != nullptr);
      Ioss::Region *region = m_region.get();

      if (!m_fieldsInitialized) {
        InputQuery iq(*region, bulk.mesh_meta_data(), m_db_purpose, m_multiStateSuffixes);

        for (stk::io::MeshField& mf : m_fields) {
          mf.set_inactive();
          iq.build_field_part_associations(mf, missingFields);
          iq.build_field_part_associations_from_grouping_entity(mf, missingFields);
        }

        m_fieldsInitialized = true;
      }

      if (time < m_startTime || time > m_stopTime)
        return 0.0;

      // Map analysis time to database time using offset, periodic, ...
      // See details in header file.
      double db_time = map_analysis_to_db_time(time);
	
      STK_ThrowErrorMsgIf (m_region.get() == nullptr,
		       "ERROR: There is no Input mesh/restart region associated with this Mesh Data.");


      // Get struct containing interval of database time(s) containing 'time'
      DBStepTimeInterval sti(region, db_time);

      STK_ThrowErrorMsgIf(!sti.exists_before && !sti.exists_after,
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
        double time_t = (*I).restore_field_data(bulk, sti, ignore_missing_fields, m_multiStateSuffixes);
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
                                                        stk::mesh::BulkData &bulk, bool useEntityListCache)
    {
      STK_ThrowErrorMsgIf(step <= 0,
                      "ERROR: Invalid step (" << step << ") requested. Value must be greater than zero.");

        STK_ThrowErrorMsgIf (m_region.get() == nullptr,
                         "There is no Input mesh/restart region associated with this Mesh Data.");

      Ioss::Region *region = m_region.get();

      int step_count = region->get_property("state_count").get_int();

      STK_ThrowErrorMsgIf(step_count == 0,
                      "ERROR: Input database '" << region->get_database()->get_filename()
                      << "' has no transient data.");

      STK_ThrowErrorMsgIf(step > step_count,
                      "ERROR: Input database '" << region->get_database()->get_filename()
                      << "'. Step " << step << " was specified, but database only has "
                      << step_count << " steps.");

      // Sort fields to ensure they are iterated in the same order on all processors.
      std::sort(m_fields.begin(), m_fields.end(), meshFieldSort);

      bool ignore_missing_fields = (missingFields != nullptr);


      if (!m_fieldsInitialized) {
        InputQuery iq(*region, bulk.mesh_meta_data(), m_db_purpose, m_multiStateSuffixes);

        for (auto & meshField : m_fields) {
          meshField.set_inactive();
          iq.build_field_part_associations(meshField, missingFields);
        }

        m_fieldsInitialized = true;
      }

      double time  = region->get_state_time(step);
      if (time < m_startTime || time > m_stopTime)
        return 0.0;

      if (useEntityListCache && not m_haveCachedEntityList) {
        for (auto & meshField : m_fields) {
          meshField.fill_entity_list_cache(bulk);
        }
        m_haveCachedEntityList = true;
      }

      double time_read = -1.0;
      for (auto & meshField : m_fields) {
        // NOTE: If the fields being restored have different settings, the time
        // value can be different for each field and this will return the value
        // of the last field.  For example, if one field is CLOSEST, one is SPECFIED,
        // and one is TIME_INTERPOLATION, then the time value to return is
        // ambiguous.  Also an issue if some of the fields are inactive.
        double time_t = meshField.restore_field_data_at_step(region, bulk, step, ignore_missing_fields,
                                                             m_multiStateSuffixes, useEntityListCache);
        if (meshField.is_active()) {
          time_read = time_t > time_read ? time_t : time_read;
        }
      }

      int current_step = region->get_current_state();
      if (current_step != -1 && current_step != static_cast<int>(step))
        region->end_state(current_step);
      if (current_step != static_cast<int>(step))
        region->begin_state(step);

      return time_read;
    }

    void InputFile::initialize_input_fields()
    {
      for (auto & meshField : m_fields) {
        meshField.set_inactive();
        meshField.clear_field_parts();
      }

      m_fieldsInitialized = false;
    }

    double InputFile::get_max_time(bool useCache)
    {
      if(m_maxTime == -std::numeric_limits<double>::max()) {
        m_maxTime = get_input_ioss_region()->get_max_time().second;
        return m_maxTime;
      }

      auto mode = get_input_ioss_region()->get_database()->open_create_behavior();
      bool openForReadAndWrite = (mode == Ioss::DB_APPEND || mode == Ioss::DB_APPEND_GROUP);
      if(useCache && !openForReadAndWrite) {
        return m_maxTime;
      }

      m_maxTime = get_input_ioss_region()->get_max_time().second;
      return m_maxTime;
    }
  }
}
