/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_io_IossBridge_hpp
#define stk_io_IossBridge_hpp

#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

// TODO: remove this and uses of USE_FEMMETADATA once the migration is done (srkenno@sandia.gov)

#define USE_FEMMETADATA
#ifdef USE_FEMMETADATA
#include <stk_mesh/fem/FEMHelpers.hpp>
#endif

#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <Ioss_DBUsage.h>
#include <Ioss_Field.h>
#include <Ioss_SideBlock.h>
#include <Ioss_ElementTopology.h>

namespace Ioss {
  class Region;
  class GroupingEntity;
  class EntityBlock;
  class SideSet;
  class Field;
  class ElementTopology;
}

class CellTopologyData;

namespace stk {

/**
 * The stk::io namespace contains functions related to the
 * transfer of data between the Ioss classes and the stk::mesh
 * classes.  These functions do not provide a total turnkey mesh
 * reading or results writing capability; rather, they provide
 * helper functions for the application to use which make it
 * easier to read and/or write the data.  The application has full
 * control over the mesh reading and results/restart writing.
 *
 * The basic flow chart for reading mesh data from the Ioss and
 * defining and populating an stk::mesh structure is shown in the
 * use_cases/io_example.cpp file.
 */
namespace io {

/** \addtogroup stk_io_module
 * \{
 */
/** Returns true if the Ioss 'entity' should be a 'part' in the
 * analysis mesh.  Returns false if the application is only using
 * a subset of the database entities and this entity is not to be
 * used. The function checks whether the "omitted" property exists
 * and if it does, whether the value of the property is "1". The
 * "omitted" property is set by the application during parsing or
 * pre-mesh reading time.
 */
bool include_entity(Ioss::GroupingEntity *entity);

/** This is the primary function used by an application to define
 *	the stk::mesh which corresponds to the Ioss mesh read from the
 *	finite element model (e.g. exodusII file). For all entities in
 *	the passed in 'entities' list, the function will determine
 *	whether the entity should be included (see stk::io::include_entity()),
 *	and it will then delcare a part corresponding to the
 *	entity. It also adds the io_part attribute (see
 *	stk::io::define_output_db()) which will cause the part to be output to a
 *	results or restart file.
 */
template <typename T>
void default_part_processing(const std::vector<T*> &entities,
			     stk::mesh::fem::FEMMetaData &fem_meta)
{
  for(size_t i=0; i < entities.size(); i++) {
    T* entity = entities[i];
    internal_part_processing(entity, fem_meta);
  }
}

//! \deprecated
template <typename T>
void default_part_processing(const std::vector<T*> &entities, stk::mesh::MetaData &meta,
                             stk::mesh::EntityRank)
{
  stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(meta);
  default_part_processing (entities, fem_meta);
}

/** Given the newly created Ioss::Region 'io_region', define the
 * model corresponding to the stk::mesh 'bulk_data'.  If the
 * optional 'input_region' is passed as an argument, then
 * synchronize all names and ids found on 'input_region' to the
 * output region 'io_region'.  The routine will query all parts
 * in 'bulk_data' and if they are io_parts (define by the existance
 * of the IOPartAttribute attribute on the part), then a
 * corresponding Ioss entity will be defined.  This routine only
 * deals with the non-transient portion of the model; no transient
 * fields are defined at this point.
 */
void define_output_db( Ioss::Region & io_region,
                       const mesh::BulkData& bulk_data,
                       const Ioss::Region *input_region = NULL);

/** Given an Ioss::Region 'io_region' which has already had its
 * metadata defined via 'define_output_db()' call; transfer all bulk
 * data (node coordinates, element connectivity, ...) to the
 * output database that corresponds to this Ioss::Region. At
 * return, all non-transient portions of the output database will
 * have been output.
 */
void write_output_db( Ioss::Region & io_region ,
                      const mesh::BulkData& bulk);


//----------------------------------------------------------------------

/** Determine whether the field is defined on the specified part
 * and should also be associated with an Ioss::GroupingEntity for
 * input or output.
 *
 * The 'universal' part is argument currently only used for an
 * EntityRank of type 'node'. In this case, the field is only
 * valid if it does exist on 'part', but does not exist on
 * 'universal' unless 'part' and 'universal' are the same part.
 * The motivation for this is that an Ioss::NodeBlock corresponds
 * to the universal part and an Ioss::NodeSet corresponds to the
 * sub part. Without the check, all nodeblock fields would also be
 * defined on the nodeset which is not what is desired.
 *
 * The 'filter_role' only selects fields with the specified role
 * (e.g., TRANSIENT, ATTRIBUTE, ..., see Ioss documentation for
 * valid roles) unless 'add_all == true' is specified.
 *
 * The field's role is defined via a call to
 * 'stk::io::set_field_role'
 */
bool is_valid_part_field(const stk::mesh::FieldBase *field,
                         stk::mesh::EntityRank part_type,
                         stk::mesh::Part &part,
                         stk::mesh::Part &universal,
                         Ioss::Field::RoleType filter_role,
                         bool add_all = false);

/** Add all stk::Fields on the specified part of the specified
 * filter_role to the specified Ioss::GroupingEntity.  Retrieves
 * all fields; calls 'is_valid_part_field'; and adds those that
 * return true.
 */
void ioss_add_fields(stk::mesh::Part &part,
                     stk::mesh::EntityRank part_type,
                     Ioss::GroupingEntity *entity,
                     const Ioss::Field::RoleType filter_role,
                     bool add_all = false);

/**
 * For the given Ioss::GroupingEntity "entity", find all fields that
 * exist on the input database of type "role" and declare them on
 * the give stk::mesh::Part "part". The "part_type" argument
 * specifies the entity type (node, element, ...) that the field
 * should be declared on for this "part"
 *
 * The "role" will typically be either "ATTRIBUTE" or "TRANSIENT"
 *
 * This is essentially the complement of the 'ioss_add_fields'
 * function.
 */
void define_io_fields(Ioss::GroupingEntity *entity,
                      Ioss::Field::RoleType role,
                      stk::mesh::Part &part,
                      stk::mesh::EntityRank part_type);

/** Given an Ioss::ElementTopolgy, return the corresponding
 *  CellTopologyData. If a corresponding topology is not found, a
 *  runtime error exception will be thrown.
 */
const CellTopologyData *map_topology_ioss_to_cell(const Ioss::ElementTopology *topology);

/** Given a CellTopologyData and a spatial_dimension, return the
 *	corresponding Ioss::ElementTopology. If a corresponding
 *	topology is not found, a runtime error exception will be
 *	thrown.
 */
std::string map_topology_cell_to_ioss( const CellTopologyData *cell_top,
                                       int spatial_dimension);

/**
 * For the given Ioss entity, create a vector of stk::mesh::Entity
 * pointers such that the entities in the 'entities' list match
 * the order of the entities in the Ioss entity. If there is not a
 * corresponding stk::mesh::Entity, the entry at that location
 * will be NULL. Upon return, the size of the 'entities' list
 * should match the number of entities in 'io_entity'. The returned
 * list is typically used to get/put field data from/to an
 * Ioss::GroupingEntity to/from an stk::mesh::Field.  See
 * stk::io::field_data_from_ioss() and stk::io::field_data_to_ioss() for examples.
 */
void get_entity_list(Ioss::GroupingEntity *io_entity,
                     stk::mesh::EntityRank part_type,
                     const stk::mesh::BulkData &bulk,
                     std::vector<stk::mesh::Entity*> &entities);

/**
 * Fill the specified 'field' with data from the Ioss field named
 * 'io_fld_name' on the Ioss entity 'io_entity'. The mapping from
 * the meshobjects in the Ioss io_entity to the
 * stk::mesh::Entities is given by the 'entities' list.
 */
void field_data_from_ioss(const stk::mesh::FieldBase *field,
                          std::vector<stk::mesh::Entity*> &entities,
                          Ioss::GroupingEntity *io_entity,
                          const std::string &io_fld_name);

/**
 * Extract data from the specified 'field' and put it to the Ioss
 * field named 'io_fld_name' on the Ioss entity 'io_entity'. The
 * mapping from the meshobjects in the Ioss io_entity to the
 * stk::mesh::Entities is given by the 'entities' list.
 */
void field_data_to_ioss(const stk::mesh::FieldBase *field,
                        std::vector<stk::mesh::Entity*> &entities,
                        Ioss::GroupingEntity *io_entity,
                        const std::string &io_fld_name,
			Ioss::Field::RoleType filter_role);


/** Returns the stk::mesh::Field which contains the distribution
 *	factors for the specified part 'p'. Returns NULL if there is
 *	no such field.
 */
const mesh::Field<double, mesh::ElementNode> *get_distribution_factor_field(const mesh::Part &p);

/** Defines the stk::mesh::Field which contains the distribution
 *	factors for the specified part 'p'.
 */
void set_distribution_factor_field(mesh::Part &p,
                                   const mesh::Field<double, mesh::ElementNode> &df_field);

/** Returns the Ioss::Field::RoleType of the mesh::Field 'f'.
 *	This must have earlier been defined using
 *	stk::io::set_field_role().  Returns NULL if the role was not
 *	defined.
 */
const Ioss::Field::RoleType* get_field_role(const mesh::FieldBase &f);

/** Defines the Ioss::Field::RoleType of the mesh::Field 'f' to be
 *	'role'.
 */
void set_field_role(mesh::FieldBase &f, const Ioss::Field::RoleType &role);

/** Returns whether the mesh::Part 'p' should be output to a
 *	results or restart database. Or, in other words, whether the
 *	part should have an Ioss::GroupingEntity of the correct type
 *	defined in an Ioss::Region.  The function will return true if
 *	the part 'p' has an IOPartAttribute defined on it. The
 *	attributed is defined via the stk::io::put_io_part_attribute()
 *	function.
 */
bool is_part_io_part(mesh::Part &part);

/** Define an attribute on the specified part 'part' indicating that
 * this part should be used for io.  \see is_part_io_part()
 */
void put_io_part_attribute( mesh::Part &part, Ioss::GroupingEntity *entity = NULL);

const Ioss::GroupingEntity *get_associated_ioss_entity(const mesh::Part &part);

void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::fem::FEMMetaData &meta);

void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::fem::FEMMetaData &meta);

//! \deprecated
void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::MetaData &meta);

//! \deprecated
void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::MetaData &meta);

// To minimize ifdefs for the deprecated code:
bool invalid_rank(mesh::EntityRank rank);
mesh::EntityRank part_primary_entity_rank(const mesh::Part &part);
mesh::EntityRank element_rank(const mesh::MetaData &meta);
mesh::EntityRank side_rank(const mesh::MetaData &meta);
mesh::EntityRank face_rank(const mesh::MetaData &meta);
mesh::EntityRank edge_rank(const mesh::MetaData &meta);
mesh::EntityRank node_rank(const mesh::MetaData &meta);
void set_cell_topology(mesh::Part &part, const CellTopologyData * const cell_topology);
const CellTopologyData *get_cell_topology(const mesh::Part &part);

void initialize_spatial_dimension(mesh::fem::FEMMetaData &fem_meta, size_t spatial_dimension, const std::vector<std::string> &entity_rank_names);

//! \deprecated
void initialize_spatial_dimension(mesh::MetaData &meta, size_t spatial_dimension, const std::vector<std::string> &entity_rank_names);

void get_io_field_type(const stk::mesh::FieldBase *field,
                       int num_comp, std::pair<std::string, Ioss::Field::BasicType> *result);
/**
 * \}
 */

}//namespace io
}//namespace stk
#endif

