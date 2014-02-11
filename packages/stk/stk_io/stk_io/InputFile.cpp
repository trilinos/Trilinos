#include <stk_io/InputFile.hpp>
#include <stk_io/MeshField.hpp>
#include <stk_io/IossBridge.hpp>  
#include <stk_io/DbStepTimeInterval.hpp>

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldBase.hpp> 
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction

#include <stk_util/diag/FileUtils.hpp>  // for filename_substitution

#include "Ioss_DBUsage.h"               // for DatabaseUsage, etc
#include "Ioss_PropertyManager.h"       // for PropertyManager
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_Field.h"                 // for Field, Field::BasicType, etc
#include "Ioss_GroupingEntity.h"        // for GroupingEntity
#include "Ioss_NodeBlock.h"
#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc

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
      : m_db_purpose(purpose), m_fieldsInitialized(false)
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
      : m_database(ioss_input_region->get_database()), m_region(ioss_input_region)
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
      }
    }

    void InputFile::add_input_field(const stk::io::MeshField &mesh_field)
    {
      bool fieldAlreadyExists=false;
      for (size_t i=0; i <m_fields.size(); i++) {
        if (mesh_field.field() == m_fields[i].field()) {
          fieldAlreadyExists = true;
          break;
        }
      }

      if (!fieldAlreadyExists) {
        m_fields.push_back(mesh_field);
        stk::io::set_field_role(*mesh_field.field(), Ioss::Field::TRANSIENT);
      }
    }

    void InputFile::get_global_variable_names(std::vector<std::string> &names)
    {
      ThrowErrorMsgIf (Teuchos::is_null(m_region),
		       "Attempt to read global variables before restart initialized.");
      m_region->field_describe(Ioss::Field::TRANSIENT, &names);
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

    double InputFile::read_defined_input_fields(int step,
						std::vector<stk::io::MeshField> *missing,
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
      return read_defined_input_fields(state_time, missing, bulk);
    }

    void InputFile::build_field_part_associations(stk::io::MeshField &mesh_field,
						  const stk::mesh::Part &part,
						  const stk::mesh::EntityRank rank,
						  Ioss::GroupingEntity *io_entity)
    {
      const stk::mesh::FieldBase *f = mesh_field.field();
      // Only add TRANSIENT Fields -- check role; if not present assume transient...
      const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
      if (role == NULL || *role == Ioss::Field::TRANSIENT) {
	if (stk::io::is_field_on_part(f, rank, part)) {
	  const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, rank, part);
	  std::pair<std::string, Ioss::Field::BasicType> field_type;
	  stk::io::get_io_field_type(f, res, &field_type);
	  if (field_type.second != Ioss::Field::INVALID) {
	  
	    // See if field with that name exists on io_entity...
	    const std::string &db_name = mesh_field.db_name();
	    if (io_entity->field_exists(db_name)) {
	      mesh_field.add_part(rank, io_entity);
	      mesh_field.m_singleState = (m_db_purpose == stk::io::READ_RESTART) ? false : true;
	      mesh_field.m_wasFound = true;
	    }
	  }
	}
      }
    }

    void InputFile::build_field_part_associations(stk::mesh::BulkData &bulk)
    {
      // Each input field will have a list of the Parts that the field exists on...
      // Create this list.
      
      Ioss::Region *region = m_region.get();

      // Check universal_part() NODE_RANK first...
      const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);
      {
	std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
	while (I != m_fields.end()) {
	  build_field_part_associations(*I, meta.universal_part(), stk::topology::NODE_RANK,
					region->get_node_blocks()[0]);
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
          if (entity != NULL && entity->type() != Ioss::SIDESET) {
	    std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
	    while (I != m_fields.end()) {
	      build_field_part_associations(*I, *part, rank, entity);

	      // If rank is != NODE_RANK, then see if field is defined on the nodes of this part
	      if (rank != stk::topology::NODE_RANK) {
		Ioss::GroupingEntity *node_entity = NULL;
		std::string nodes_name = part->name() + "_nodes";
		node_entity = region->get_entity(nodes_name);
		if (node_entity == NULL) {
		  node_entity = region->get_entity("nodeblock_1");
		}
		if (node_entity != NULL) {
		  build_field_part_associations(*I, *part, stk::topology::NODE_RANK, node_entity);
		}
	      }
	      ++I;
	    }
	  }
	}
      }
    }

    void InputFile::report_missing_fields(std::vector<stk::io::MeshField> *missing) const
    {
      size_t missing_fields = 0;
      std::ostringstream msg ;
      std::vector<stk::io::MeshField>::const_iterator I = m_fields.begin();
      while (I != m_fields.end()) {
	if (!(*I).m_wasFound) {
	  ++missing_fields;
	  if (missing) {
	    missing->push_back(*I);
	  }
	  else {
	    msg << "ERROR: Could not find input field '" << (*I).db_name() << "'.\n";
	  }
	}
	++I;
      }

      ThrowAssert(missing==NULL || missing_fields == missing->size());
      if (missing_fields > 0 && missing==NULL) {
	msg << "ERROR: Input field processing could not find " << missing_fields << " fields.\n";
	throw std::runtime_error( msg.str() );
      }
    }

      double InputFile::read_defined_input_fields(double time,
						  std::vector<stk::io::MeshField> *missing,
						  stk::mesh::BulkData &bulk)
      {
	// Sort fields to ensure they are iterated in the same order on all processors.
	std::sort(m_fields.begin(), m_fields.end(), meshFieldSort);

	if (!m_fieldsInitialized) {
	  std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
	  while (I != m_fields.end()) {
	    (*I).m_wasFound = false; ++I;
	  }

	  build_field_part_associations(bulk);
	  report_missing_fields(missing);
	  
	  m_fieldsInitialized = true;
	}

	ThrowErrorMsgIf (Teuchos::is_null(m_region),
			 "ERROR: There is no Input mesh/restart region associated with this Mesh Data.");

	Ioss::Region *region = m_region.get();

	// Get struct containing interval of database time(s) containing 'time'
	DBStepTimeInterval sti(region, time);

	ThrowErrorMsgIf(!sti.exists_before && !sti.exists_after,
			"ERROR: Input database '" << region->get_database()->get_filename()
			<< "' has no transient data.");

	std::vector<stk::io::MeshField>::iterator I = m_fields.begin();
	while (I != m_fields.end()) {
	  (*I).restore_field_data(bulk, sti);
	  ++I;
	}

	size_t step = sti.get_closest_step();
	int current_step = region->get_current_state();
	if (current_step != -1 && current_step != (int)step)
	  region->end_state(current_step);
	if (current_step != (int)step) 
	  region->begin_state(step);

	return region->get_state_time(step);
      }

  }
}
