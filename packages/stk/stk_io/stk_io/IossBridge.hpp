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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_DBUsage.h>                // for DatabaseUsage
#include <Ioss_Field.h>                  // for Field, Field::RoleType, etc
#include <stddef.h>                      // for size_t
#include <stk_mesh/base/Types.hpp>       // for EntityId, EntityRank
#include <stk_topology/topology.hpp>     // for topology
#include <string>                        // for string, operator<, etc
#include <utility>                       // for pair
#include <vector>                        // for vector
#include "Ioss_EntityType.h"             // for EntityType
#include "Ioss_GroupingEntity.h"
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
#include "stk_mesh/base/FieldBase.hpp"  // for FieldState
#include "stk_mesh/base/Part.hpp"        // for Part
#include "Ioss_GroupingEntity.h"                     // for GroupingEntity
#include <stk_mesh/base/MetaData.hpp>                // for MetaData, etc
#include "SidesetTranslator.hpp"
#include "stk_io/OutputParams.hpp"

namespace Ioss { class ElementTopology; }
namespace Ioss { class EntityBlock; }
namespace Ioss { class Region; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class FieldRestriction; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class Selector; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk { namespace mesh { class Part; } }

namespace Ioss {
class SideSet;
class SideBlock;
class NodeBlock;
class Field;
class GroupingEntity;
class Region;
class ElementTopology;
}

void STKIORequire(bool cond);

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
static const std::string s_ignore_disconnected_nodes("ignore_disconnected_nodes");
static const std::string s_sort_stk_parts("sort_stk_parts");

typedef std::pair<stk::mesh::EntityId, int> EntityIdToProcPair;
typedef std::vector<EntityIdToProcPair> EntitySharingInfo;

typedef std::pair<std::string, stk::mesh::Part*> FieldNameToPart;

struct FieldNameToPartLess
{
  inline bool operator()(const FieldNameToPart& lhs, const FieldNameToPart& rhs) const
  {
     if (lhs.first < rhs.first || rhs.first < lhs.first)
     {
        return lhs.first < rhs.first;
     }
     return lhs.second->name() < rhs.second->name();
  }
};

typedef std::vector<FieldNameToPart> FieldNameToPartVector;

stk::mesh::Part *getPart(const stk::mesh::MetaData& meta_data, const std::string& name);

bool is_valid_for_output(const stk::mesh::Part &part, const stk::mesh::Selector *output_selector = nullptr);
void get_selected_nodes(const stk::mesh::BulkData &bulk,
                        Ioss::Region &io_region,
                        const stk::mesh::Selector &selector,
                        stk::mesh::EntityVector &nodes);
size_t count_selected_nodes(const stk::mesh::BulkData &bulk,
                            Ioss::Region &io_region,
                            const stk::mesh::Selector &selector);
bool node_is_connected_to_local_element(const stk::mesh::BulkData &bulk, stk::mesh::Entity node);

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
void default_part_processing(const std::vector<T*> &entities, stk::mesh::MetaData &meta)
{
  for(size_t i=0; i < entities.size(); i++) {
    T* entity = entities[i];
    internal_part_processing(entity, meta);
  }
}

//! \deprecated
template <typename T>
void default_part_processing(const std::vector<T*> &entities, stk::mesh::MetaData &meta, const stk::mesh::EntityRank)
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
void define_output_db( stk::io::OutputParams &params,
                       const std::vector<std::vector<int>> &attributeOrdering,
                       const Ioss::Region *input_region = nullptr);

/** Given an Ioss::Region 'io_region' which has already had its
 * metadata defined via 'define_output_db()' call; transfer all bulk
 * data (node coordinates, element connectivity, ...) to the
 * output database that corresponds to this Ioss::Region. At
 * return, all non-transient portions of the output database will
 * have been output.
 */
void write_output_db( stk::io::OutputParams &params);

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
    m_field(my_field),
    m_dbName(my_db_name),
    m_variableType(my_field != nullptr ? my_field->entity_rank() : stk::topology::INVALID_RANK),
    m_useAlias(true),
    m_wasFound(false),
    m_forceNodeblockOutput(false) {}

  FieldAndName(stk::mesh::FieldBase *my_field, const std::string& my_db_name, stk::mesh::EntityRank my_var_type) :
    m_field(my_field),
    m_dbName(my_db_name),
    m_variableType(my_var_type),
    m_useAlias(true),
    m_wasFound(false),
    m_forceNodeblockOutput(false) {}

  stk::mesh::FieldBase *field() const {return m_field;};
  std::string db_name() const {return m_dbName;}
  void set_db_name(const std::string &name) {m_dbName = name;}
  stk::mesh::EntityRank type() const {return m_variableType;}
  void set_use_alias(bool useAlias) { m_useAlias = useAlias; }
  bool get_use_alias() const { return m_useAlias; }
private:
  stk::mesh::FieldBase *m_field;
  std::string m_dbName;
  stk::mesh::EntityRank m_variableType;
  bool m_useAlias;
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
                     std::vector<FieldAndName> &namedFields);

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

FieldNameToPartVector get_var_names(Ioss::Region &region, Ioss::EntityType type, stk::mesh::MetaData& meta);

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
stk::topology map_ioss_topology_to_stk(const Ioss::ElementTopology *topology, unsigned mesh_spatial_dimension);

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

std::string get_stated_field_name(const std::string &field_base_name, stk::mesh::FieldState state_identifier,
                                  std::vector<std::string>* multiStateSuffixes=nullptr);

bool field_state_exists_on_io_entity(const std::string& db_name, const stk::mesh::FieldBase* field, stk::mesh::FieldState state_identifier,
                                     Ioss::GroupingEntity *io_entity, std::vector<std::string>* multiStateSuffixes=nullptr);

bool all_field_states_exist_on_io_entity(const std::string& db_name, const stk::mesh::FieldBase* field, Ioss::GroupingEntity *io_entity,
                                         std::vector<stk::mesh::FieldState> &missing_states, std::vector<std::string>* multiStateSuffixes=nullptr);

void multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                     const stk::mesh::FieldBase *field,
                                     std::vector<stk::mesh::Entity> &entity_list,
                                     Ioss::GroupingEntity *io_entity,
                                     const std::string &name,
                                     const size_t state_count,
                                     bool ignore_missing_fields = false,
                                     std::vector<std::string>* multiStateSuffixes=nullptr);

void subsetted_multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
					       const stk::mesh::FieldBase *field,
					       std::vector<stk::mesh::Entity> &entity_list,
					       Ioss::GroupingEntity *io_entity,
					       const stk::mesh::Part *stk_part,
					       const std::string &name,
					       const size_t state_count,
					       bool ignore_missing_fields = false,
					       std::vector<std::string>* multiStateSuffixes=nullptr);

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

/** Define an alternate name to use for the part on output
 */
void set_alternate_part_name(stk::mesh::Part& part, const std::string& altPartName);
std::string get_alternate_part_name(stk::mesh::Part& part);
bool has_alternate_part_name(stk::mesh::Part& part);

/** Define original topology type to use for the part on output
 */
void set_original_topology_type(stk::mesh::Part& part);
void set_original_topology_type(stk::mesh::Part& part, const std::string& origTopo);
std::string get_original_topology_type(stk::mesh::Part& part);
bool has_original_topology_type(stk::mesh::Part& part);

/** Define an attribute on the specified part 'part' indicating that
 * this part should be used for io.  \see is_part_io_part()
 */
void put_io_part_attribute( mesh::Part &part);

/** Remove the existing attribute on the specified part 'part' that indicates that
 * this part should be used for io.  \see is_part_io_part()
 */
void remove_io_part_attribute(mesh::Part &part);

size_t db_api_int_size(const Ioss::GroupingEntity *entity);

void initialize_spatial_dimension(mesh::MetaData &meta, size_t spatial_dimension, const std::vector<std::string> &entity_rank_names);

void write_file_for_subdomain(const std::string &baseFilename,
                              int index_subdomain,
                              int num_subdomains,
                              int global_num_nodes,
                              int global_num_elems,
                              stk::mesh::BulkData& bulkData,
                              const EntitySharingInfo &nodeSharingInfo,
                              int numSteps = -1,
                              double timeStep = 0.0);

void write_output_db_node_block(stk::io::OutputParams &params);

void write_output_db_element_blocks(stk::io::OutputParams &params);

struct FieldType {
    std::string name;
    Ioss::Field::BasicType type;
    unsigned copies;

    FieldType() : name(""), type(Ioss::Field::INVALID), copies(0) {}

    FieldType(std::string inputName, Ioss::Field::BasicType inputType, unsigned inputCopies) : name(inputName), type(inputType), copies(inputCopies) {}
};

void get_io_field_type(const stk::mesh::FieldBase *field,
                       const stk::mesh::FieldRestriction &res,
                       FieldType *result);

const stk::mesh::Part* get_parent_element_block(const stk::mesh::BulkData &bulk,
                                                const Ioss::Region &ioRegion,
                                                const std::string& name);

template <typename INT>
void fill_data_for_side_block( OutputParams &params,
                               Ioss::GroupingEntity & io ,
                               mesh::Part * const part ,
                               const Ioss::ElementTopology *element_topology,
                               std::vector<INT> &elem_side_ids,
                               stk::mesh::EntityVector &sides)
{
    ThrowRequireMsg(io.type() == Ioss::SIDEBLOCK, "Input GroupingEntity must be of type Ioss::SIDEBLOCK");

    stk::topology stk_elem_topology = map_ioss_topology_to_stk(element_topology, params.bulk_data().mesh_meta_data().spatial_dimension());

    const stk::mesh::Part *parentElementBlock = get_parent_element_block(params.bulk_data(), params.io_region(), part->name());

    fill_element_and_side_ids(params, io, part, parentElementBlock, stk_elem_topology, sides, elem_side_ids);
}

}//namespace io
}//namespace stk
#endif

