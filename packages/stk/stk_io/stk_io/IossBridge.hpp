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

#ifndef stk_io_IossBridge_hpp
#define stk_io_IossBridge_hpp

#include <Ioss_DBUsage.h>               // for DatabaseUsage
#include <Ioss_Field.h>                 // for Field, Field::RoleType, etc
#include <stddef.h>                     // for size_t, NULL
#include <stk_mesh/base/Types.hpp>      // for EntityRank
#include <stk_topology/topology.hpp>    // for topology
#include <string>                       // for string, basic_string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
namespace Ioss { class ElementTopology; }
namespace Ioss { class EntityBlock; }
namespace Ioss { class GroupingEntity; }
namespace Ioss { class Region; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class FieldRestriction; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }

namespace Ioss {
class SideSet;
class Field;
}

void STKIORequire(bool cond);
void STKIORequireMsg(bool cond, const std::string &msg);

namespace stk {
  namespace mesh {
  }

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
bool include_entity(const Ioss::GroupingEntity *entity);

void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::MetaData &meta);

void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::MetaData &meta);

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
                             stk::mesh::MetaData &meta)
{
  for(size_t i=0; i < entities.size(); i++) {
    T* entity = entities[i];
    internal_part_processing(entity, meta);
  }
}

//! \deprecated
template <typename T>
void default_part_processing(const std::vector<T*> &entities, stk::mesh::MetaData &meta,
                             const stk::mesh::EntityRank)
{
  default_part_processing (entities, meta);
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
 *
 * \param[in] sort_stk_parts Force a sorted order on the stk_mesh
 * parts so all pieces of a decomposed mesh have consistent part
 * ordering.  Normally not necessary, since MetaData is created
 * from a parallel-consistent data base.  It is useful in cases
 * such as streaming refinement where each piece of a parallel-
 * decomposed mesh is read in sequentially and parts are thus
 * possibly created in different orderings.
 *
 * \param[in] use_nodeset_for_part_node_fields is used to determine
 * how to output nodal fields that may be defined on a higher-rank
 * part (e.g., element block or sideset).  If the argument is true,
 * then a nodeset containing the nodes of that part is defined and
 * the fields will be output on that nodeset. If false, the fields
 * are output on all nodes in the model and zero-filled where the
 * field is not defined.
 */
void define_output_db( Ioss::Region & io_region,
                       const mesh::BulkData& bulk_data,
                       const Ioss::Region *input_region = NULL,
                       const stk::mesh::Selector *subset_selector = NULL,
                       bool sort_stk_parts = false,
                       bool use_nodeset_for_part_node_fields = true);


/** Given an Ioss::Region 'io_region' which has already had its
 * metadata defined via 'define_output_db()' call; transfer all bulk
 * data (node coordinates, element connectivity, ...) to the
 * output database that corresponds to this Ioss::Region. At
 * return, all non-transient portions of the output database will
 * have been output.
 */
void write_output_db( Ioss::Region & io_region ,
                      const mesh::BulkData& bulk,
                      const stk::mesh::Selector *subset_selector = NULL);


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
 * valid roles)
 *
 * The field's role is defined via a call to
 * 'stk::io::set_field_role'
 */
bool is_valid_part_field(const stk::mesh::FieldBase *field,
                         const stk::mesh::EntityRank part_type,
                         const stk::mesh::Part &part,
                         const Ioss::Field::RoleType filter_role);

bool is_field_on_part(const stk::mesh::FieldBase *field,
		      const stk::mesh::EntityRank part_type,
		      const stk::mesh::Part &part);

/** Add all stk::Fields on the specified part of the specified
 * filter_role to the specified Ioss::GroupingEntity.  Retrieves
 * all fields; calls 'is_valid_part_field'; and adds those that
 * return true.
 */
struct FieldAndName
{
public:
  FieldAndName(stk::mesh::FieldBase *my_field, const std::string& my_db_name) :
    m_field(my_field), m_dbName(my_db_name),
    m_wasFound(false), m_forceNodeblockOutput(false) {}
  stk::mesh::FieldBase *field() const {return m_field;};
  std::string db_name() const {return m_dbName;}
  void set_db_name(const std::string &name) {m_dbName = name;}
private:
  stk::mesh::FieldBase *m_field;
  std::string m_dbName;
public:
  bool m_wasFound;
  // Field is not defined on UNIVERSAL part, but we still want to output it on the nodeblock.
  // This is done to output, for example, nodal fields that exist on an element block without
  // creating a nodeset for the nodes of the element block.
  mutable bool m_forceNodeblockOutput; 
};

std::string get_field_name(const stk::mesh::FieldBase &f, Ioss::DatabaseUsage dbUsage);
void getNamedFields(const stk::mesh::MetaData &meta, Ioss::GroupingEntity *io_entity, std::vector<FieldAndName> &namedFields);

void ioss_add_fields(const stk::mesh::Part &part,
                     const stk::mesh::EntityRank part_type,
                     Ioss::GroupingEntity *entity,
                     const Ioss::Field::RoleType filter_role);

void ioss_add_fields(const stk::mesh::Part &part,
                     const stk::mesh::EntityRank part_type,
                     Ioss::GroupingEntity *entity,
                     const std::vector<FieldAndName> &namedFields);

/**
 * Iterate over all Ioss entities in the input mesh Ioss Region and
 * define a stk field for each transient field found.  The stk field
 * will have the same name as the field on the database.
 *
 * Note that all transient fields found on the mesh database will
 * have a corresponding stk field defined.  If you want just a
 * selected subset of the database fields defined in the stk mesh,
 * you need to define the fields manually.
 *
 * To populate the stk field with data from the database, call
 * StkMeshIoBroker::process_input_request().
 *
 */
void define_input_fields(Ioss::Region &region,  stk::mesh::MetaData &meta);

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
 *  stk::topology. If a corresponding topology is not found, a
 *  runtime error exception will be thrown.
 */
stk::topology map_ioss_topology_to_stk(const Ioss::ElementTopology *topology);

/** Given a stk::topology, return the
 *	corresponding Ioss::ElementTopology string. If a corresponding
 *	topology is not found, a runtime error exception will be
 *	thrown.
 */
std::string map_stk_topology_to_ioss(stk::topology topo);

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
                     std::vector<stk::mesh::Entity> &entities);


/**
 * Delete the selector property (if it exists) which is used to get
 * the entity list on an output database's GroupingEntitys
 */
template <typename T> void delete_selector_property(std::vector<T> &entities);
void delete_selector_property(Ioss::Region &io_region);
void delete_selector_property(Ioss::GroupingEntity *io_entity);

std::string get_stated_field_name(const std::string &field_base_name, stk::mesh::FieldState state_identifier);

void multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                          const stk::mesh::FieldBase *field,
                          std::vector<stk::mesh::Entity> &entity_list,
                          Ioss::GroupingEntity *io_entity,
                          const std::string &name,
                          const size_t state_count);

void subsetted_multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
					       const stk::mesh::FieldBase *field,
					       std::vector<stk::mesh::Entity> &entity_list,
					       Ioss::GroupingEntity *io_entity,
					       const stk::mesh::Part *stk_part,
					       const std::string &name,
					       const size_t state_count);

/**
 * Fill the specified 'field' with data from the Ioss field named
 * 'io_fld_name' on the Ioss entity 'io_entity'. The mapping from
 * the meshobjects in the Ioss io_entity to the
 * stk::mesh::Entities is given by the 'entities' list.
 */
void field_data_from_ioss(const stk::mesh::BulkData& mesh,
                          const stk::mesh::FieldBase *field,
                          std::vector<stk::mesh::Entity> &entities,
                          Ioss::GroupingEntity *io_entity,
                          const std::string &io_fld_name);

void subsetted_field_data_from_ioss(const stk::mesh::BulkData& mesh,
				    const stk::mesh::FieldBase *field,
				    std::vector<stk::mesh::Entity> &entities,
				    Ioss::GroupingEntity *io_entity,
				    const stk::mesh::Part *stk_part,
				    const std::string &io_fld_name);

void multistate_field_data_to_ioss(const stk::mesh::BulkData& mesh,
                        const stk::mesh::FieldBase *field,
                        std::vector<stk::mesh::Entity> &entities,
                        Ioss::GroupingEntity *io_entity,
                        const std::string &io_fld_name,
                        Ioss::Field::RoleType filter_role,
                        const size_t state_count);
/**
 * Extract data from the specified 'field' and put it to the Ioss
 * field named 'io_fld_name' on the Ioss entity 'io_entity'. The
 * mapping from the meshobjects in the Ioss io_entity to the
 * stk::mesh::Entities is given by the 'entities' list.
 */
void field_data_to_ioss(const stk::mesh::BulkData& mesh,
                        const stk::mesh::FieldBase *field,
                        std::vector<stk::mesh::Entity> &entities,
                        Ioss::GroupingEntity *io_entity,
                        const std::string &io_fld_name,
                        Ioss::Field::RoleType filter_role);


/** Returns the stk::mesh::Field which contains the distribution
 *	factors for the specified part 'p'. Returns NULL if there is
 *	no such field.
 */
const mesh::FieldBase *get_distribution_factor_field(const mesh::Part &p);

/** Defines the stk::mesh::Field which contains the distribution
 *	factors for the specified part 'p'.
 */
void set_distribution_factor_field(mesh::Part &p,
                                   const mesh::FieldBase &df_field);

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
bool is_part_io_part(const mesh::Part &part);

/** Define an attribute on the specified part 'part' indicating that
 * this part should be used for io.  \see is_part_io_part()
 */
void put_io_part_attribute( mesh::Part &part, Ioss::GroupingEntity *entity = NULL);

/** Remove the existing attribute on the specified part 'part' that indicates that
 * this part should be used for io.  \see is_part_io_part()
 */
void remove_io_part_attribute(mesh::Part &part);

const Ioss::GroupingEntity *get_associated_ioss_entity(const mesh::Part &part);

size_t db_api_int_size(const Ioss::GroupingEntity *entity);

mesh::EntityRank part_primary_entity_rank(const mesh::Part &part);

void initialize_spatial_dimension(mesh::MetaData &meta, size_t spatial_dimension, const std::vector<std::string> &entity_rank_names);

void get_io_field_type(const stk::mesh::FieldBase *field,
                       const stk::mesh::FieldRestriction &res,
                       std::pair<std::string, Ioss::Field::BasicType> *result);
/**
 * \}
 */

}//namespace io
}//namespace stk
#endif

