/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string.h>
#include <iostream>
#include <complex>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_NullEntity.h>

#include <stk_util/util/tokenize.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_topology/topology.hpp>
#include <Shards_BasicTopologies.hpp>

namespace {

  const std::string internal_selector_name = "_stk_io_internal_selector";
  const std::string base_stk_part_name = "_base_stk_part_name";
  const std::string block_nodes_suffix = "_nodes";

stk::mesh::EntityRank get_entity_rank(const Ioss::GroupingEntity *entity,
                                      const stk::mesh::MetaData &meta)
{
  switch (entity->type()) {
  case Ioss::NODEBLOCK:
    return stk::mesh::MetaData::NODE_RANK;

  case Ioss::NODESET:
    return stk::mesh::MetaData::NODE_RANK;

  case Ioss::ELEMENTBLOCK:
    return stk::mesh::MetaData::ELEMENT_RANK;

  case Ioss::SUPERELEMENT:
    return stk::mesh::MetaData::ELEMENT_RANK;

  case Ioss::SIDESET:
    {
      const Ioss::SideSet *sset = dynamic_cast<const Ioss::SideSet*>(entity);
      assert(sset != NULL);
      int my_rank = sset->max_parametric_dimension();
      if (my_rank == 2)
	return stk::mesh::MetaData::FACE_RANK;
      if (my_rank == 1)
	return stk::mesh::MetaData::EDGE_RANK;
      if (my_rank == 0)
	return stk::mesh::MetaData::NODE_RANK;
      else
        return stk::mesh::InvalidEntityRank;
    }

  case Ioss::SIDEBLOCK:
    {
      const Ioss::SideBlock *sblk = dynamic_cast<const Ioss::SideBlock*>(entity);
      assert(sblk != NULL);
      int rank = sblk->topology()->parametric_dimension();
      if (rank == 2)
	return stk::mesh::MetaData::FACE_RANK;
      if (rank == 1)
	return stk::mesh::MetaData::EDGE_RANK;
      if (rank == 0)
	return stk::mesh::MetaData::NODE_RANK;
      else
        return stk::mesh::InvalidEntityRank;
    }
  default:
    return stk::mesh::InvalidEntityRank;
  }
}


template <typename T>
const stk::mesh::FieldBase *declare_ioss_field_internal(stk::mesh::MetaData &meta,
                                                        stk::mesh::EntityRank type,
                                                        stk::mesh::Part &part,
                                                        const Ioss::Field &io_field,
                                                        bool use_cartesian_for_scalar, T /*dummy*/)
{
  stk::mesh::FieldBase *field_ptr = NULL;
  std::string field_type = io_field.transformed_storage()->name();
  std::string name = io_field.get_name();
  size_t num_components = io_field.transformed_storage()->component_count();

  if (field_type == "scalar" || num_components == 1) {
    if (!use_cartesian_for_scalar) {
      stk::mesh::Field<double> & field = meta.declare_field<stk::mesh::Field<double> >(name);
      stk::mesh::put_field(field, type, part);
      field_ptr = &field;
    } else {
      stk::mesh::Field<double, stk::mesh::Cartesian> & field =
        meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(name);
      stk::mesh::put_field(field, type, part, 1);
      field_ptr = &field;
    }
  }
  else if (field_type == "vector_2d") {
    stk::mesh::Field<double, stk::mesh::Cartesian> & field =
      meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(name);
    stk::mesh::put_field(field, type, part, 2);
    field_ptr = &field;
  }
  else if (field_type == "vector_3d") {
    stk::mesh::Field<double, stk::mesh::Cartesian> & field =
      meta.declare_field<stk::mesh::Field<double,
      stk::mesh::Cartesian> >(name);
    stk::mesh::put_field(field, type, part, 3);
    field_ptr = &field;
  }
  else if (field_type == "sym_tensor_33") {
    stk::mesh::Field<double, stk::mesh::SymmetricTensor> & field =
      meta.declare_field<stk::mesh::Field<double,
      stk::mesh::SymmetricTensor> >(name);
    stk::mesh::put_field(field, type, part, 6);
    field_ptr = &field;
  }
  else if (field_type == "full_tensor_36") {
    stk::mesh::Field<double, stk::mesh::FullTensor> & field =
      meta.declare_field<stk::mesh::Field<double,
      stk::mesh::FullTensor> >(name);
    stk::mesh::put_field(field, type, part, 9);
    field_ptr = &field;
  }
  else {
    // Just create a field with the correct number of components...
    stk::mesh::Field<double,shards::ArrayDimension> & field =
      meta.declare_field<stk::mesh::Field<double,shards::ArrayDimension> >(name);
    stk::mesh::put_field(field, type, part, num_components);
    field_ptr = &field;
  }

  if (field_ptr != NULL) {
    stk::io::set_field_role(*field_ptr, io_field.get_role());
  }
  return field_ptr;
}

const stk::mesh::FieldBase *declare_ioss_field(stk::mesh::MetaData &meta,
                                               stk::mesh::EntityRank type,
                                               stk::mesh::Part &part,
                                               const Ioss::Field &io_field,
                                               bool use_cartesian_for_scalar)
{
  const stk::mesh::FieldBase *field_ptr = NULL;
  if (io_field.get_type() == Ioss::Field::INTEGER) {
    field_ptr = declare_ioss_field_internal(meta, type, part, io_field, use_cartesian_for_scalar, (int)1);
  } else if (io_field.get_type() == Ioss::Field::REAL) {
    field_ptr = declare_ioss_field_internal(meta, type, part, io_field, use_cartesian_for_scalar, (double)1.0);
  } else if (io_field.get_type() == Ioss::Field::COMPLEX) {
    field_ptr = declare_ioss_field_internal(meta, type, part, io_field, use_cartesian_for_scalar, std::complex<double>(0.0,0.0));
  } else {
    std::ostringstream errmsg;
    errmsg << "ERROR: Unrecognized field type for IO field '"
           << io_field.get_name() << "'.";
    throw std::runtime_error(errmsg.str());
  }
  return field_ptr;
}

template <typename T>
void internal_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                   const Ioss::Field &io_field,
                                   const stk::mesh::FieldBase *field,
                                   std::vector<stk::mesh::Entity> &entities,
                                   Ioss::GroupingEntity *io_entity,
                                   T /*dummy */)
{
  size_t field_component_count = io_field.transformed_storage()->component_count();

  std::vector<T> io_field_data;
  size_t io_entity_count = io_entity->get_field_data(io_field.get_name(), io_field_data);
  assert(io_field_data.size() == entities.size() * field_component_count);

  size_t entity_count = entities.size();

  if (io_entity_count != entity_count) {
    std::ostringstream errmsg;
    errmsg << "ERROR: Field count mismatch for IO field '"
           << io_field.get_name() << "'. The IO system has " << io_entity_count
           << " entries, but the stk:mesh system has " << entity_count
           << " entries. The two counts must match.";
    throw std::runtime_error(errmsg.str());
  }

  for (size_t i=0; i < entity_count; ++i) {
    if (mesh.is_valid(entities[i])) {
      T *fld_data = (T*)mesh.field_data(*field, entities[i]);
      assert(fld_data != NULL);
      for(size_t j=0; j<field_component_count; ++j) {
        fld_data[j] = io_field_data[i*field_component_count+j];
      }
    }
  }
}

template <typename T>
void internal_field_data_to_ioss(const stk::mesh::BulkData& mesh,
                                 const Ioss::Field &io_field,
                                 const stk::mesh::FieldBase *field,
                                 std::vector<stk::mesh::Entity> &entities,
                                 Ioss::GroupingEntity *io_entity,
                                 T /*dummy */)
{
  size_t field_component_count = io_field.transformed_storage()->component_count();
  size_t entity_count = entities.size();

  std::vector<T> io_field_data(entity_count*field_component_count);

  for (size_t i=0; i < entity_count; ++i) {
    if (mesh.is_valid(entities[i])) {
      T *fld_data = (T*)mesh.field_data(*field, entities[i]);
      if (fld_data != NULL) {
        for(size_t j=0; j<field_component_count; ++j) {
          io_field_data[i*field_component_count+j] = fld_data[j];
        }
      }
    }
  }

  size_t io_entity_count = io_entity->put_field_data(io_field.get_name(), io_field_data);
  assert(io_field_data.size() == entities.size() * field_component_count);

  if (io_entity_count != entity_count) {
    std::ostringstream errmsg;
    errmsg << "ERROR: Field count mismatch for IO field '"
           << io_field.get_name() << "'. The IO system has " << io_entity_count
           << " entries, but the stk:mesh system has " << entity_count
           << " entries. The two counts must match.";
    throw std::runtime_error(errmsg.str());
  }
}

bool will_output_lower_rank_fields(const stk::mesh::Part &part, stk::mesh::EntityRank rank)
{
  // Determine whether the part 'part' needs to output any fields of rank 'rank'
  // to the output database.

  const std::vector<stk::mesh::FieldBase *> &fields = stk::mesh::MetaData::get(part).get_fields();
  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I ; ++I ;
    if (stk::io::is_valid_part_field(f, rank, part, Ioss::Field::TRANSIENT, false)) {
      return true;
    }
  }
  return false;
}
}//namespace <empty>

namespace stk {
namespace io {

size_t db_api_int_size(const Ioss::GroupingEntity *entity)
{
  return entity->get_database()->int_byte_size_api();
}

stk::mesh::EntityRank part_primary_entity_rank(const stk::mesh::Part &part)
{
  if (mesh::MetaData::get(part).universal_part() == part) {
    return stk::mesh::MetaData::NODE_RANK;
  }
  else {
    return part.primary_entity_rank();
  }
}

void initialize_spatial_dimension(stk::mesh::MetaData & meta, size_t spatial_dimension,
                                  const std::vector<std::string> &entity_rank_names)
{
  if (!meta.is_initialized() ) {
    meta.initialize(spatial_dimension, entity_rank_names);
  }
}

bool is_valid_part_field(const stk::mesh::FieldBase *field,
                         const stk::mesh::EntityRank part_type,
                         const stk::mesh::Part &part,
                         const Ioss::Field::RoleType filter_role,
                         bool add_all)
{
  const Ioss::Field::RoleType *role = stk::io::get_field_role(*field);

  if (!add_all && role == NULL) {
	return false;
  }

  if (role != NULL && *role != filter_role)
	return false;

  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(part);
  const stk::mesh::FieldBase::Restriction &res = field->restriction(part_type, part);
  if (res.dimension() > 0) {
    // The field exists on the current 'part'.  Now check (for
    // node types only) whether the 'part' is *either* the
    // universal_part() *or* the field *doesn't* exist on the
    // universal part...
    // Note that for "node" type parts, the IO database has a part
    // (nodeblock) that corresponds to the universal_part(), so
    // fields defined on all nodes are output on the nodeblock and
    // fields defined on only a subset of the parts should be
    // output on that subset which maps to a nodeset.  For other
    // part types, there is no IO entity that corresponds to a
    // universal part, so fields must be output on the part they
    // exist on.  There may be a problem if we start using element
    // sets ..., but wait until we get to that point; current code
    // works with current entity set.
    if (part_type != stk::mesh::MetaData::NODE_RANK || part == meta.universal_part()) {
      return true;
    }

    const stk::mesh::FieldBase::Restriction &res_universe = field->restriction(part_type, meta.universal_part());
    if (res_universe.dimension() <= 0) {
      // Field exists on current part, but not on the universal
      // set (and this part is not the universal part)
      return true;
    }
  }
  return false;
}

void get_io_field_type(const stk::mesh::FieldBase *field,
                       const stk::mesh::FieldRestriction &res,
                       std::pair<std::string, Ioss::Field::BasicType> *result)
{
  static const std::string invalid("invalid");
  static const std::string scalar("scalar");
  static const std::string vector_2d("vector_2d");
  static const std::string vector_3d("vector_3d");
  static const std::string quaternion_2d("quaternion_2d");
  static const std::string quaternion_3d("quaternion_3d");
  static const std::string full_tensor_36("full_tensor_36");
  static const std::string full_tensor_32("full_tensor_32");
  static const std::string full_tensor_22("full_tensor_22");
  static const std::string full_tensor_16("full_tensor_16");
  static const std::string full_tensor_12("full_tensor_12");
  static const std::string sym_tensor_33("sym_tensor_33");
  static const std::string sym_tensor_31("sym_tensor_31");
  static const std::string sym_tensor_21("sym_tensor_21");
  static const std::string sym_tensor_13("sym_tensor_13");
  static const std::string sym_tensor_11("sym_tensor_11");
  static const std::string sym_tensor_10("sym_tensor_10");
  static const std::string asym_tensor_03("asym_tensor_03");
  static const std::string asym_tensor_02("asym_tensor_02");
  static const std::string asym_tensor_01("asym_tensor_01");
  static const std::string matrix_22("matrix_22");
  static const std::string matrix_33("matrix_33");

  const unsigned rank = field->rank();
  const shards::ArrayDimTag * const * const tags = field->dimension_tags();

  result->second = Ioss::Field::INVALID;

  if ( field->type_is<double>() ) {
	result->second = Ioss::Field::REAL;
  }
  else if ( field->type_is<int>() ) {
	result->second = Ioss::Field::INTEGER;
  }

  if ( 0 == rank ) {
	result->first = scalar ;
  }
  else if ( 1 == rank ) {
	size_t num_comp = res.stride(0);
	if ( tags[0] == & stk::mesh::Cartesian::tag() && 1 == num_comp ) {
	  result->first = scalar ;
	}
	else if ( tags[0] == & stk::mesh::Cartesian::tag() && 2 == num_comp ) {
	  result->first = vector_2d ;
	}
	else if ( tags[0] == & stk::mesh::Cartesian::tag() && 3 == num_comp ) {
	  result->first = vector_3d ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 9 == num_comp ) {
	  result->first = full_tensor_36 ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 5 == num_comp ) {
	  result->first = full_tensor_32 ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 4 == num_comp ) {
	  result->first = full_tensor_22 ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 3 == num_comp ) {
	  result->first = full_tensor_12 ;
	}
	else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 6 == num_comp ) {
	  result->first = sym_tensor_33 ;
	}
	else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 4 == num_comp ) {
	  result->first = sym_tensor_31 ;
	}
	else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 3 == num_comp ) {
	  result->first = sym_tensor_21 ;
	}
  }

  if ( result->first.empty() ) {
	size_t num_comp = res.stride(rank-1);
	std::ostringstream tmp ;
	tmp << "Real[" << num_comp << "]" ;
	result->first = tmp.str();
  }
}

//----------------------------------------------------------------------
const Ioss::GroupingEntity *get_associated_ioss_entity(const mesh::Part &part)
{
  const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
  if (!entity || entity->type() == Ioss::INVALID_TYPE) {
	return NULL;
  } else {
	return entity;
  }
}

void put_io_part_attribute(mesh::Part & part, Ioss::GroupingEntity *entity)
{
  if (part.attribute<Ioss::GroupingEntity>() != NULL) {
	std::string msg = "stk::io::put_io_part_attribute( ";
	msg += part.name();
	msg += " ) FAILED:";
	msg += " io_part_attribute is already defined";
	throw std::runtime_error( msg );
  }

  mesh::MetaData & meta = mesh::MetaData::get(part);
  if (entity) {
	meta.declare_attribute_no_delete(part, entity);
  } else {
	Ioss::GroupingEntity *attr = new Ioss::NullEntity();
	meta.declare_attribute_with_delete(part, attr);
  }
}

void remove_io_part_attribute(mesh::Part & part)
{
  const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
  if (entity != NULL) {
    mesh::MetaData & meta = mesh::MetaData::get(part);
    bool success = meta.remove_attribute(part, entity);
    if (!success) {
      std::string msg = "stk::io::remove_io_part_attribute( ";
      msg += part.name();
      msg += " ) FAILED:";
      msg += " meta.remove_attribute(..) returned failure.";
      throw std::runtime_error( msg );
    }

    if (entity->type() == Ioss::INVALID_TYPE) {
      delete entity;
    }
  }
}

stk::topology map_ioss_topology_to_stk( const Ioss::ElementTopology *topology)
{
  for (stk::topology topo=stk::topology::BEGIN_TOPOLOGY; topo < stk::topology::END_TOPOLOGY; ++topo) {
    if (topology->is_alias(topo.name())) {
      return topo;
    }
  }
  std::string tmpCopy = topology->name().substr(0,5);
  sierra::make_lower(tmpCopy);
  if (tmpCopy=="super")
  {
      return stk::create_superelement_topology(topology->number_nodes());
  }

  return stk::topology::INVALID_TOPOLOGY;
}

std::string map_stk_topology_to_ioss(stk::topology topo)
{
  Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(topo.name(), true);
  return ioss_topo != NULL ? ioss_topo->name() : "invalid";
}

void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::MetaData &meta)
{
  if (include_entity(entity)) {
    mesh::EntityRank type = get_entity_rank(entity, meta);
    stk::mesh::Part & part = meta.declare_part(entity->name(), type);
    if (entity->property_exists("id")) {
        part.set_id(entity->get_property("id").get_int());
    }
    stk::io::put_io_part_attribute(part, entity);
  }
}

void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::MetaData &meta)
{
  if (include_entity(entity)) {
    mesh::EntityRank type = get_entity_rank(entity, meta);
    stk::mesh::Part * part = NULL;
    part = &meta.declare_part(entity->name(), type);
    if (entity->property_exists("id")) {
        part->set_id(entity->get_property("id").get_int());
    }
    stk::io::put_io_part_attribute(*part, entity);

    const Ioss::ElementTopology *topology = entity->topology();
    // Check spatial dimension of the element topology here so we
    // can issue a more meaningful error message.  If the
    // dimension is bad and we continue to the following calls,
    // there is an exception and we get unintelligible (to the
    // user) error messages.  Could also do a catch...

    if (entity->type() == Ioss::ELEMENTBLOCK) {
      assert(topology != NULL);
      if (topology->spatial_dimension() < (int)meta.spatial_dimension()) {
	// NOTE: The comparison is '<' and not '!=' since a 2D mesh
	// can contain a "3d" element -- a Beam is both a 2D and
	// 3D element...

	std::ostringstream msg ;
	msg << "\n\nERROR: Element Block " << entity->name()
            << " contains " << topology->name() << " elements with spatial dimension "
            << topology->spatial_dimension()
            << "\n       which does not match the spatial dimension of the model which is "
            << meta.spatial_dimension() << "\n\n";
	    throw std::runtime_error( msg.str() );
      }
    }

    stk::topology stk_topology = map_ioss_topology_to_stk(topology);
    stk::mesh::set_topology(*part, stk_topology);
    stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE, *part, type);
  }
}

//----------------------------------------------------------------------
/** Add all stk::Fields on the entities of the specified part_type
 *  on the specified part of the specified role * to the specified
 *  Ioss::GroupingEntity
 */
void ioss_add_fields(const stk::mesh::Part &part,
                     const stk::mesh::EntityRank part_type,
                     Ioss::GroupingEntity *entity,
                     const Ioss::Field::RoleType filter_role,
                     const bool add_all)
{
  const stk::mesh::MetaData & meta = mesh::MetaData::get(part);
  const std::vector<mesh::FieldBase*> &fields = meta.get_fields();

  std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
	const stk::mesh::FieldBase *f = *I; ++I;
	if (stk::io::is_valid_part_field(f, part_type, part, filter_role, add_all)) {
	  const stk::mesh::FieldBase::Restriction &res = f->restriction(part_type, part);
	  std::pair<std::string, Ioss::Field::BasicType> field_type;
	  get_io_field_type(f, res, &field_type);
	  if (field_type.second != Ioss::Field::INVALID) {
	    size_t entity_size = entity->get_property("entity_count").get_int();
	    const std::string& name = f->name();
	    entity->field_add(Ioss::Field(name, field_type.second, field_type.first,
                                      filter_role, entity_size));
	  }
	}
  }
}

/**
 * For the given Ioss::GroupingEntity "entity", find all fields that
 * exist on the input database of type "role" and declare them on
 * the give stk::mesh::Part "part". The "part_type" argument
 * specifies the entity type (node, element, ...) that the field
 * should be declared on for this "part"
 *
 * The "role" will typically be either "ATTRIBUTE" or "TRANSIENT"
 */
void define_io_fields(Ioss::GroupingEntity *entity,
                      Ioss::Field::RoleType role,
                      stk::mesh::Part &part,
                      stk::mesh::EntityRank part_type)
{
  stk::mesh::MetaData &meta = mesh::MetaData::get(part);

  bool use_cartesian_for_scalar = false;
  if (role == Ioss::Field::ATTRIBUTE)
    use_cartesian_for_scalar = true;

  Ioss::NameList names;
  entity->field_describe(role, &names);

  for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    // \todo IMPLEMENT Need a field selection mechanism and a field naming
    // (ioss_name -> stk::name)  For now, select all and give the
    // stk field the same name as the ioss field.

    // Skip the attribute field that is named "attribute"
    if (*I == "attribute" && names.size() > 1)
      continue;

    // \todo IMPLEMENT Need to determine whether these are
    // multi-state fields or constant, or interpolated, or ...
    Ioss::Field io_field = entity->get_field(*I);
    declare_ioss_field(meta, part_type, part, io_field, use_cartesian_for_scalar);
  }
}

template <typename T>
void delete_selector_property(const std::vector<T> &entities)
{
  for(size_t i=0; i < entities.size(); i++) {
    delete_selector_property(entities[i]);
  }
}

void delete_selector_property(Ioss::Region &region)
{
  // Iterate all Ioss::GroupingEntity types on the io_region and
  // if the have a property named 'selector' of type 'pointer',
  // delete the pointer and remove the property.
  delete_selector_property(region.get_node_blocks());
  delete_selector_property(region.get_element_blocks());
  delete_selector_property(region.get_nodesets());
  delete_selector_property(region.get_commsets());

  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *sset = *it;
    delete_selector_property(*it);
    delete_selector_property(sset->get_side_blocks());
  }
}

void delete_selector_property(Ioss::GroupingEntity *io_entity)
{
  // If the Ioss::GroupingEntity has a property named 'selector' of
  // type 'pointer', delete the pointer and remove the property.
  if (io_entity->property_exists(internal_selector_name)) {
    mesh::Selector *select = reinterpret_cast<mesh::Selector*>(io_entity->get_property(internal_selector_name).get_pointer());
    delete select;
    io_entity->property_erase(internal_selector_name);
  }
}

template <typename INT>
void get_entity_list(Ioss::GroupingEntity *io_entity,
                     stk::mesh::EntityRank part_type,
                     const stk::mesh::BulkData &bulk,
                     std::vector<stk::mesh::Entity> &entities, INT /*dummy*/)
{
  std::vector<INT> ids ;
  io_entity->get_field_data("ids", ids);

  size_t count = ids.size();
  entities.reserve(count);

  for(size_t i=0; i<count; ++i) {
    entities.push_back(bulk.get_entity( part_type, ids[i] ));
  }
}

void get_entity_list(Ioss::GroupingEntity *io_entity,
                     stk::mesh::EntityRank part_type,
                     const stk::mesh::BulkData &bulk,
                     std::vector<stk::mesh::Entity> &entities)
{
  if (io_entity->get_database()->is_input()) {
    if (db_api_int_size(io_entity) == 4) {
      get_entity_list(io_entity, part_type, bulk, entities, (int)0);
    } else {
      get_entity_list(io_entity, part_type, bulk, entities, (int64_t)0);
    }
  } else {
    // Output database...
    assert(io_entity->property_exists(internal_selector_name));

    mesh::Selector *select = reinterpret_cast<mesh::Selector*>(io_entity->get_property(internal_selector_name).get_pointer());
    get_selected_entities(*select, bulk.buckets(part_type), entities);
  }
}

void field_data_from_ioss(const stk::mesh::BulkData& mesh,
                          const stk::mesh::FieldBase *field,
                          std::vector<stk::mesh::Entity> &entities,
                          Ioss::GroupingEntity *io_entity,
                          const std::string &io_fld_name)
{
  /// \todo REFACTOR Need some additional compatibility checks between
  /// Ioss field and stk::mesh::Field; better error messages...

  if (field != NULL && io_entity->field_exists(io_fld_name)) {
	const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
	if (field->type_is<double>()) {
	  internal_field_data_from_ioss(mesh, io_field, field, entities, io_entity,
                                    static_cast<double>(1.0));
	} else if (field->type_is<int>()) {
	  // Make sure the IO field type matches the STK field type.
	  // By default, all IO fields are created of type 'double'
	  if (db_api_int_size(io_entity) == 4) {
	    if (io_field.get_type() != Ioss::Field::INTEGER) {
	      Ioss::Field &tmp = const_cast<Ioss::Field&>(io_field);
	      tmp.reset_type(Ioss::Field::INTEGER);
	    }
	    internal_field_data_from_ioss(mesh, io_field, field, entities, io_entity,
                                      static_cast<int>(1));
	  } else {
	    if (io_field.get_type() != Ioss::Field::INT64) {
	      Ioss::Field &tmp = const_cast<Ioss::Field&>(io_field);
	      tmp.reset_type(Ioss::Field::INT64);
	    }
	    internal_field_data_from_ioss(mesh, io_field, field, entities, io_entity,
                                      static_cast<int64_t>(1));
	  }
	}
  }
}

void field_data_to_ioss(const stk::mesh::BulkData& mesh,
                        const stk::mesh::FieldBase *field,
                        std::vector<stk::mesh::Entity> &entities,
                        Ioss::GroupingEntity *io_entity,
                        const std::string &io_fld_name,
                        Ioss::Field::RoleType filter_role)
{
  /// \todo REFACTOR Need some additional compatibility checks between
  /// Ioss field and stk::mesh::Field; better error messages...

  if (field != NULL && io_entity->field_exists(io_fld_name)) {
	const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
	if (io_field.get_role() == filter_role) {
	  if (field->type_is<double>()) {
	    internal_field_data_to_ioss(mesh, io_field, field, entities, io_entity,
                                    static_cast<double>(1.0));
	  } else if (field->type_is<int>()) {
	    if (io_field.get_type() != Ioss::Field::INTEGER) {
	      Ioss::Field &tmp = const_cast<Ioss::Field&>(io_field);
	      tmp.reset_type(Ioss::Field::INTEGER);
	    }
	    // FIX 64?
	    internal_field_data_to_ioss(mesh, io_field, field, entities, io_entity,
                                    static_cast<int>(1));
	  }
	}
  }
}

//----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Returns true if 'entity' should be a 'part' in the analysis mesh.
// Returns false if the application is only using a subset of the
// database entities and this entity is not to be used. The
// "omitted" property is set by the application during parsing or
// pre-mesh reading time.
bool include_entity(const Ioss::GroupingEntity *entity)
{
  assert(entity);

  // Check whether entity has "omitted" property...
  bool omitted = (entity->property_exists("omitted")) &&
	(entity->get_property("omitted").get_int() == 1);

  return !omitted;
}

namespace {

void define_side_block(const stk::mesh::BulkData &bulk,
                       const stk::mesh::Selector &selector,
                       stk::mesh::Part &part,
                       Ioss::SideSet *sset,
                       int spatial_dimension)
{
  stk::mesh::EntityRank type = part.primary_entity_rank();
  const stk::mesh::EntityRank siderank = stk::mesh::MetaData::get(part).side_rank();
  const stk::mesh::EntityRank edgerank = stk::mesh::MetaData::EDGE_RANK;
  ThrowRequire(type == siderank || type == edgerank);

  stk::topology side_topology = part.topology();
  std::string io_topo = map_stk_topology_to_ioss(side_topology);
  std::string element_topo_name = "unknown";

  // Get sideblock parent element topology quantities...
  // Try to decode from part name...
  std::vector<std::string> tokens;
  stk::util::tokenize(part.name(), "_", tokens);
  if (tokens.size() >= 4) {
    // Name of form: "name_eltopo_sidetopo_id" or
    //               "name_block_id_sidetopo_id"
    // "name" is typically "surface".
    const Ioss::ElementTopology *element_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-3], true);
    if (element_topo != NULL) {
      element_topo_name = element_topo->name();
    }
  }

  size_t side_count = count_selected_entities(selector, bulk.buckets(type));

  Ioss::SideBlock *side_block = new Ioss::SideBlock( sset->get_database() ,
                                                     part.name() ,
                                                     io_topo, element_topo_name, side_count);
  assert(sset->get_side_block(part.name()) == NULL);
  sset->add(side_block);

  const mesh::FieldBase *df = get_distribution_factor_field(part);
  if (df != NULL) {
    int nodes_per_side = side_topology.num_nodes();
    std::string storage_type = "Real[";
    storage_type += Ioss::Utils::to_string(nodes_per_side);
    storage_type += "]";
    side_block->field_add(Ioss::Field("distribution_factors", Ioss::Field::REAL, storage_type,
                                      Ioss::Field::MESH, side_count));
  }

  mesh::Selector *select = new mesh::Selector(selector);
  side_block->property_add(Ioss::Property(internal_selector_name, select, false));
  side_block->property_add(Ioss::Property(base_stk_part_name, part.name()));

  // Add the attribute fields.
  ioss_add_fields(part, part_primary_entity_rank(part), side_block, Ioss::Field::ATTRIBUTE);
}

void define_side_blocks(stk::mesh::Part &part,
                        const stk::mesh::BulkData &bulk_data,
                        Ioss::SideSet *sset,
                        stk::mesh::EntityRank type,
                        int spatial_dimension,
                        const stk::mesh::Selector *anded_selector)
{
  mesh::MetaData & meta = mesh::MetaData::get(part);
  ThrowRequire(type == stk::mesh::MetaData::FACE_RANK || stk::mesh::MetaData::EDGE_RANK);

  const stk::mesh::PartVector &blocks = part.subsets();
  if (blocks.size() > 0) {
    for (size_t j = 0; j < blocks.size(); j++) {
      mesh::Part & side_block_part = *blocks[j];
      mesh::Selector selector = meta.locally_owned_part() & side_block_part;
      if (anded_selector) selector &= *anded_selector;

      define_side_block(bulk_data, selector, side_block_part,
			sset, spatial_dimension);
    }
  } else {
    mesh::Selector selector = meta.locally_owned_part() & part;
    if (anded_selector) selector &= *anded_selector;
    define_side_block(bulk_data, selector, part, sset, spatial_dimension);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
#ifndef NDEBUG
// Used in an assert statement...
bool field_has_expected_size(const stk::mesh::FieldBase* field, unsigned expected_size)
{
  const stk::mesh::FieldBase::RestrictionVector& restrictions = field->restrictions();
  for(size_t i=0; i<restrictions.size(); ++i) {
    unsigned field_size = restrictions[i].stride(field->rank());
    if (field_size != expected_size) {
      return false;
    }
  }

  return true;
}
#endif

void define_node_block(stk::mesh::Part &part,
                       const stk::mesh::BulkData &bulk,
                       Ioss::Region &io_region,
                       const stk::mesh::Selector *anded_selector)
{
  //--------------------------------
  // Set the spatial dimension:
  mesh::MetaData & meta = mesh::MetaData::get(part);

  //We now get spatial-dim from meta.spatial_dimension() rather than getting
  //it from the coordinate-field's restriction onto the universal part.
  //This is because some codes (sierra framework) don't put the coordinate
  //field on the universal part. (framework puts it on active and inactive parts)
  //As a debug-mode error-check, we'll assert that 'field_has_expected_size'.

  const int spatial_dim = meta.spatial_dimension();
  io_region.property_add( Ioss::Property("spatial_dimension", spatial_dim));

  ThrowAssertMsg(field_has_expected_size(meta.coordinate_field(), spatial_dim), "IossBridge define_node_block ERROR, coordinate field doesn't have the correct size. (Should match spatial-dimension.)");

  //--------------------------------
  // Create the special universal node block:

  mesh::Selector all_selector = meta.globally_shared_part() | meta.locally_owned_part();
  if (anded_selector) all_selector &= *anded_selector;

  mesh::Selector own_selector = meta.locally_owned_part();
  if (anded_selector) own_selector &= *anded_selector;

  int64_t all_nodes = count_selected_entities(all_selector, bulk.buckets(stk::mesh::MetaData::NODE_RANK));
  int64_t own_nodes = count_selected_entities(own_selector, bulk.buckets(stk::mesh::MetaData::NODE_RANK));

  const std::string name("nodeblock_1");

  Ioss::NodeBlock * const nb = new Ioss::NodeBlock(io_region.get_database(),
                                                   name, all_nodes, spatial_dim);
  io_region.add( nb );

  mesh::Selector *node_select = new mesh::Selector(all_selector);
  nb->property_add(Ioss::Property(internal_selector_name, node_select, false));
  nb->property_add(Ioss::Property(base_stk_part_name, part.name()));

  // Add locally-owned property...
  nb->property_add(Ioss::Property("locally_owned_count", own_nodes));
  // Add the attribute fields.
  ioss_add_fields(part, part_primary_entity_rank(part), nb, Ioss::Field::ATTRIBUTE);
}

void define_node_set(stk::mesh::Part &part,
                     const std::string &name,
                     const stk::mesh::BulkData &bulk,
                     Ioss::Region &io_region,
                     const stk::mesh::Selector *anded_selector)
{
  mesh::MetaData & meta = mesh::MetaData::get(part);

  mesh::Selector all_selector = (meta.globally_shared_part() | meta.locally_owned_part()) & part;
  if (anded_selector) all_selector &= *anded_selector;

  mesh::Selector own_selector = meta.locally_owned_part() & part;
  if (anded_selector) own_selector &= *anded_selector;

  int64_t all_nodes = count_selected_entities(all_selector, bulk.buckets(stk::mesh::MetaData::NODE_RANK));
  int64_t own_nodes = count_selected_entities(own_selector, bulk.buckets(stk::mesh::MetaData::NODE_RANK));

  Ioss::NodeSet * const ns = new Ioss::NodeSet( io_region.get_database(), name, all_nodes);
  io_region.add(ns);

  ns->property_add(Ioss::Property("locally_owned_count", own_nodes));

  mesh::Selector *select = new mesh::Selector(all_selector);
  ns->property_add(Ioss::Property(internal_selector_name, select, false));
  ns->property_add(Ioss::Property(base_stk_part_name, part.name()));

  // Add the attribute fields.
  ioss_add_fields(part, stk::mesh::MetaData::NODE_RANK, ns, Ioss::Field::ATTRIBUTE);
}

void define_element_block(stk::mesh::Part &part,
                          const stk::mesh::BulkData &bulk,
                          Ioss::Region &io_region,
                          const stk::mesh::Selector *anded_selector,
                          bool use_nodeset_for_nodal_fields)
{
  mesh::MetaData & meta = mesh::MetaData::get(part);

  stk::topology topo = part.topology();
  if (topo == stk::topology::INVALID_TOPOLOGY) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part.name() << " returned INVALID from get_topology()";
    throw std::runtime_error( msg.str() );
  }

  mesh::Selector selector = meta.locally_owned_part() & part;
  if (anded_selector) selector &= *anded_selector;
  const size_t num_elems = count_selected_entities( selector, bulk.buckets(stk::mesh::MetaData::ELEMENT_RANK));

  // Defer the counting of attributes until after we define the
  // element block so we can count them as we add them as fields to
  // the element block
  Ioss::ElementBlock *eb = new Ioss::ElementBlock(io_region.get_database() ,
                                                  part.name() ,
                                                  map_stk_topology_to_ioss(part.topology()),
                                                  num_elems);
  io_region.add(eb);

  mesh::Selector *select = new mesh::Selector(selector);
  eb->property_add(Ioss::Property(internal_selector_name, select, false));
  eb->property_add(Ioss::Property(base_stk_part_name, part.name()));

  // Add the attribute fields.
  ioss_add_fields(part, part_primary_entity_rank(part), eb, Ioss::Field::ATTRIBUTE);

  // Check whether there are any transient fields defined on the nodes of this elementblock
  // that are to be output.  If so, create a nodeset named "part.name()"+block_nodes_suffix
  // and output the fields on that nodeset...
  if (use_nodeset_for_nodal_fields &&
      will_output_lower_rank_fields(part, stk::mesh::MetaData::NODE_RANK)) {
    std::string nodes_name = part.name() + block_nodes_suffix;
    define_node_set(part, nodes_name, bulk, io_region, anded_selector);
  }
}

void define_communication_maps(const stk::mesh::BulkData &bulk,
                               Ioss::Region &io_region,
                               const stk::mesh::Selector *anded_selector)
{
  if (bulk.parallel_size() > 1) {
    const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);
    const std::string cs_name("node_symm_comm_spec");

    mesh::Selector selector = meta.globally_shared_part();
    if (anded_selector) selector &= *anded_selector;

    std::vector<mesh::Entity> entities;
    get_selected_entities(selector, bulk.buckets(stk::mesh::MetaData::NODE_RANK), entities);

    size_t size = 0;
    for (size_t i=0; i < entities.size(); i++) {
      for ( stk::mesh::PairIterEntityComm ec = bulk.entity_comm(bulk.entity_key(entities[i])); ! ec.empty() ; ++ec ) {
	if (ec->ghost_id == 0) {
	  size++;
	}
      }
    }

    Ioss::DatabaseIO *dbo = io_region.get_database();
    Ioss::CommSet *io_cs = new Ioss::CommSet(dbo, cs_name, "node", size);
    io_region.add(io_cs);

    mesh::Selector *select = new mesh::Selector(selector);
    io_cs->property_add(Ioss::Property(internal_selector_name, select, false));
  }
}

void define_side_set(stk::mesh::Part &part,
                     const stk::mesh::BulkData &bulk,
                     Ioss::Region &io_region,
                     const stk::mesh::Selector *anded_selector,
                     bool use_nodeset_for_nodal_fields)
{
  const stk::mesh::EntityRank si_rank = mesh::MetaData::get(part).side_rank();

  bool create_sideset = true;
  if (part.subsets().empty()) {
    // Only define a sideset for this part if its superset part is
    // not a side-containing part..  (i.e., this part is not a subset part
    // in a surface...)
    const stk::mesh::PartVector &supersets = part.supersets();
    for (size_t i=0; i < supersets.size(); i++) {
      if (is_part_io_part(*supersets[i]) &&
          (supersets[i]->primary_entity_rank() == stk::topology::FACE_RANK ||
           supersets[i]->primary_entity_rank() == stk::topology::EDGE_RANK)) {
        create_sideset = false;
        break;
      }
    }
  }

  if (create_sideset) {
    Ioss::SideSet * const ss = new Ioss::SideSet(io_region.get_database(), part.name());

    io_region.add(ss);
    int spatial_dim = io_region.get_property("spatial_dimension").get_int();
    define_side_blocks(part, bulk, ss, si_rank, spatial_dim, anded_selector);

    if (use_nodeset_for_nodal_fields &&
        will_output_lower_rank_fields(part, stk::mesh::MetaData::NODE_RANK)) {
      std::string nodes_name = part.name() + block_nodes_suffix;
      define_node_set(part, nodes_name, bulk, io_region, anded_selector);
    }
  }
}

} // namespace <blank>

struct part_compare {
  bool operator() (stk::mesh::Part *i, stk::mesh::Part *j) { return (i->name() < j->name()); }
};

void define_output_db(Ioss::Region & io_region ,
                      const mesh::BulkData &bulk_data,
                      const Ioss::Region *input_region,
                      const stk::mesh::Selector *anded_selector,
                      const bool sort_stk_parts,
                      const bool use_nodeset_for_part_node_fields)
{
  io_region.begin_mode( Ioss::STATE_DEFINE_MODEL );

  const mesh::MetaData & meta_data = mesh::MetaData::get(bulk_data);
  define_node_block(meta_data.universal_part(), bulk_data, io_region, anded_selector);

  // All parts of the meta data:
  const mesh::PartVector *parts = NULL;
  mesh::PartVector all_parts_sorted;

  const mesh::PartVector & all_parts = meta_data.get_parts();
  // sort parts so they go out the same on all processors (srk: this was induced by streaming refine)
  if (sort_stk_parts) {
    all_parts_sorted = all_parts;
    std::sort(all_parts_sorted.begin(), all_parts_sorted.end(), part_compare());
    parts = &all_parts_sorted;
  } else {
    parts = &all_parts;
  }

  for (mesh::PartVector::const_iterator i = parts->begin(); i != parts->end(); ++i) {
    mesh::Part * const part = *i;

    if (is_part_io_part(*part)) {
      if (part->primary_entity_rank() == mesh::InvalidEntityRank)
        continue;
      else if (part->primary_entity_rank() == stk::mesh::MetaData::NODE_RANK)
        define_node_set(*part, part->name(), bulk_data, io_region, anded_selector);
      else if (part->primary_entity_rank() == stk::mesh::MetaData::ELEMENT_RANK)
        define_element_block(*part, bulk_data, io_region, anded_selector,
            use_nodeset_for_part_node_fields);
      else if (part->primary_entity_rank() == stk::mesh::MetaData::FACE_RANK)
        define_side_set(*part, bulk_data, io_region, anded_selector,
            use_nodeset_for_part_node_fields);
      else if (part->primary_entity_rank() == stk::mesh::MetaData::EDGE_RANK)
        define_side_set(*part, bulk_data, io_region, anded_selector,
            use_nodeset_for_part_node_fields);
    }
  }

  define_communication_maps(bulk_data, io_region, anded_selector);

  if (input_region != NULL)
    io_region.synchronize_id_and_name(input_region, true);

  // for streaming refinement, each "pseudo-processor" doesn't know about others, so we pick a sort order
  //   and use it for all pseudo-procs - the original_block_order property is used to set the order
  //   on all procs.
  if (sort_stk_parts) {
    int offset=0;
    for (mesh::PartVector::const_iterator i = parts->begin(); i != parts->end(); ++i) {

      mesh::Part * const part = *i ;

      if (is_part_io_part(*part)) {
        if (part->primary_entity_rank() == mesh::InvalidEntityRank)
          continue;
        else if (part->primary_entity_rank() == stk::mesh::MetaData::ELEMENT_RANK) {
          Ioss::GroupingEntity *element_block = io_region.get_entity(part->name());
          if (element_block) {
            if (element_block->property_exists("original_block_order")) {
              element_block->property_erase("original_block_order");
            }
            element_block->property_add(Ioss::Property("original_block_order", offset));
            ++offset;
          }
        }
      }
    }
  }

  io_region.end_mode( Ioss::STATE_DEFINE_MODEL );
}

//----------------------------------------------------------------------

namespace {

size_t get_entities(stk::mesh::Part &part, stk::mesh::EntityRank type,
                    const stk::mesh::BulkData &bulk,
                    std::vector<mesh::Entity> &entities,
                    bool include_shared,
                    const stk::mesh::Selector *anded_selector)
{
  mesh::MetaData & meta = mesh::MetaData::get(part);

  mesh::Selector own_share = meta.locally_owned_part();
  if (include_shared)
    own_share |= meta.globally_shared_part();

  mesh::Selector selector = part & own_share;
  if (anded_selector) selector &= *anded_selector;

  get_selected_entities(selector, bulk.buckets(type), entities);
  return entities.size();
}


template <typename INT>
void write_side_data_to_ioss( Ioss::GroupingEntity & io ,
                              mesh::Part * const part ,
                              const mesh::BulkData & bulk_data,
                              const stk::mesh::Selector *anded_selector, INT /*dummy*/ )
{
  //std::cout << "tmp write_side_data_to_ioss part= " << part->name() << std::endl;
  const mesh::MetaData & meta_data = mesh::MetaData::get(*part);

  std::vector<mesh::Entity> sides ;
  stk::mesh::EntityRank type = part_primary_entity_rank(*part);
  size_t num_sides = get_entities(*part, type, bulk_data, sides, false, anded_selector);

  std::vector<INT> elem_side_ids; elem_side_ids.reserve(num_sides*2);

  for(size_t i=0; i<num_sides; ++i) {

    const mesh::Entity side = sides[i] ;
    mesh::Entity const *side_elem = bulk_data.begin_elements(side);
    mesh::ConnectivityOrdinal const *side_ordinal = bulk_data.begin_element_ordinals(side);

    // Which element to use?
    // Any locally owned element that has the "correct" orientation

    const size_t num_side_elem = bulk_data.num_elements(side);
    size_t suitable = std::numeric_limits<size_t>::max();

    for ( size_t j = 0 ; (j < num_side_elem) && (suitable >= num_side_elem) ; ++j )
    {
      const mesh::Entity elem = side_elem[j];

      if ( bulk_data.bucket(elem).member( meta_data.locally_owned_part() ) &&
           (num_side_elem == 1 || bulk_data.element_side_polarity(elem, side, side_ordinal[j])) )
      {
        suitable = j;
      }
    }

    if (suitable >= num_side_elem)
    {
      std::ostringstream oss;
      oss << "ERROR, no suitable element found";
      throw std::runtime_error(oss.str());
    }

    elem_side_ids.push_back(bulk_data.identifier(side_elem[suitable]));
    elem_side_ids.push_back(side_ordinal[suitable] + 1) ; // Ioss is 1-based, mesh is 0-based.
  }

  const size_t num_side_written = io.put_field_data("element_side",elem_side_ids);

  if ( num_sides != num_side_written ) {
    std::ostringstream msg ;

    msg << "stk::io::write_side_data_to_ioss FAILED for " ;
    msg << io.name();
    msg << " in Ioss::GroupingEntity::put_field_data:" ;
    msg << " num_sides = " << num_sides ;
    msg << " , num_side_written = " << num_side_written ;
    throw std::runtime_error( msg.str() );
  }

  const mesh::FieldBase *df = get_distribution_factor_field(*part);
  if (df != NULL) {
    field_data_to_ioss(bulk_data, df, sides, &io, "distribution_factors", Ioss::Field::MESH);
  }

  const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
  std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const mesh::FieldBase *f = *I ; ++I ;
    const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
    if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
      stk::io::field_data_to_ioss(bulk_data, f, sides, &io, f->name(), Ioss::Field::ATTRIBUTE);
    }
  }
}

//----------------------------------------------------------------------
template <typename INT>
void output_node_block(Ioss::NodeBlock &nb,
                       stk::mesh::Part &part,
                       const stk::mesh::BulkData &bulk,
                       const stk::mesh::Selector *anded_selector, INT /*dummy*/)
{
  //----------------------------------
  // Exactly one node block to obtain the nodal coordinates and ids:
  // Note that the "ids" field of the nodes needs to be written
  // before any other bulk data that uses node ids since it sets up
  // the global->local mapping of nodes for the output database.
  // Similarly for the element "ids" field related to bulk data
  // using element ids.
  std::vector<mesh::Entity> nodes ;
  size_t num_nodes = get_entities(part, stk::mesh::MetaData::NODE_RANK,
				  bulk, nodes, true, anded_selector);

  std::vector<INT> node_ids; node_ids.reserve(num_nodes);
  for(size_t i=0; i<num_nodes; ++i) {
    const mesh::Entity node = nodes[i] ;
    node_ids.push_back(bulk.identifier(node));
  }

  size_t num_ids_written = nb.put_field_data("ids", node_ids);
  if ( num_nodes != num_ids_written) {
    std::ostringstream msg ;
    msg << " FAILED in Ioss::NodeBlock::put_field_data:" ;
    msg << " num_nodes = " << num_nodes ;
    msg << " , num_ids_written = " << num_ids_written ;
    throw std::runtime_error( msg.str() );
  }

  if (nb.get_database()->needs_shared_node_information()) {
    std::vector<INT> owning_processor; owning_processor.reserve(num_nodes);
    for(size_t i=0; i<num_nodes; ++i) {
      owning_processor.push_back(bulk.parallel_owner_rank(nodes[i]));
    }
    nb.put_field_data("owning_processor", owning_processor);
  }

  const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
  const mesh::FieldBase *coord_field = meta_data.coordinate_field();
  assert(coord_field != NULL);
  field_data_to_ioss(bulk, coord_field, nodes, &nb, "mesh_model_coordinates", Ioss::Field::MESH);

  const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
  std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const mesh::FieldBase *f = *I ; ++I ;
    if (stk::io::is_valid_part_field(f, part_primary_entity_rank(part), part,
				     Ioss::Field::ATTRIBUTE, false)) {
      stk::io::field_data_to_ioss(bulk, f, nodes, &nb, f->name(), Ioss::Field::ATTRIBUTE);
    }
  }
}

template <typename INT>
void output_element_block(Ioss::ElementBlock *block,
                          const stk::mesh::BulkData &bulk,
                          const stk::mesh::Selector *anded_selector, INT /*dummy*/)
{
  const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
  const std::string& name = block->name();
  mesh::Part* part = meta_data.get_part(name);

  assert(part != NULL);
  std::vector<mesh::Entity> elements;
  stk::mesh::EntityRank type = part_primary_entity_rank(*part);
  size_t num_elems = get_entities(*part, type, bulk, elements, false, anded_selector);

  stk::topology topo = part->topology();
  if (topo == stk::topology::INVALID_TOPOLOGY) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part->name() << " returned INVALID from get_topology()";
    throw std::runtime_error( msg.str() );
  }
  size_t nodes_per_elem = topo.num_nodes();

  std::vector<INT> elem_ids; elem_ids.reserve(num_elems);
  std::vector<INT> connectivity; connectivity.reserve(num_elems*nodes_per_elem);

  for (size_t i = 0; i < num_elems; ++i) {

    elem_ids.push_back(bulk.identifier(elements[i]));

    stk::mesh::Entity const * elem_nodes = bulk.begin_nodes(elements[i]);

    for (size_t j = 0; j < nodes_per_elem; ++j) {
      connectivity.push_back(bulk.identifier(elem_nodes[j]));
    }
  }

  const size_t num_ids_written = block->put_field_data("ids", elem_ids);
  const size_t num_con_written = block->put_field_data("connectivity", connectivity);

  if ( num_elems != num_ids_written || num_elems != num_con_written ) {
    std::ostringstream msg ;
    msg << " FAILED in Ioss::ElementBlock::put_field_data:" << std::endl ;
    msg << "  num_elems = " << num_elems << std::endl ;
    msg << "  num_ids_written = " << num_ids_written << std::endl ;
    msg << "  num_connectivity_written = " << num_con_written << std::endl ;
    throw std::runtime_error( msg.str() );
  }

  stk::mesh::EntityRank elem_rank = stk::mesh::MetaData::ELEMENT_RANK;
  const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
  std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const mesh::FieldBase *f = *I ; ++I ;
    const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
    if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
      const mesh::FieldBase::Restriction &res = f->restriction(elem_rank, *part);
      if (res.dimension() > 0) {
        stk::io::field_data_to_ioss(bulk, f, elements, block, f->name(), Ioss::Field::ATTRIBUTE);
      }
    }
  }
}

template <typename INT>
void output_node_set(Ioss::NodeSet *ns, const stk::mesh::BulkData &bulk,
                     const stk::mesh::Selector *anded_selector, INT /*dummy*/)
{
  const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
  const std::string& name = ns->name();
  stk::mesh::Part* part = meta_data.get_part(name);

  // If part is null, then it is possible that this nodeset is a "viz nodeset" which
  // means that it is a nodeset containing the nodes of an element block.
  // See if there is a property base_stk_part_name and if so, get the part with
  // that name.
  if (part == NULL) {
    if (ns->property_exists(base_stk_part_name)) {
      std::string base_name = ns->get_property(base_stk_part_name).get_string();
      part = meta_data.get_part(base_name);
    }
    if (part == NULL) {
      std::ostringstream msg ;
      msg << " FAILED in Ioss::NodeSet::output_node_set:"
        << " Could not find stk part corresponding to nodeset named '"
	  << name << "'";
      throw std::runtime_error( msg.str() );
    }
  }

  std::vector<stk::mesh::Entity> nodes ;
  size_t num_nodes = get_entities(*part, stk::mesh::MetaData::NODE_RANK, bulk, nodes, true, anded_selector);

  std::vector<INT> node_ids; node_ids.reserve(num_nodes);
  for(size_t i=0; i<num_nodes; ++i) {
    const stk::mesh::Entity node = nodes[i] ;
    node_ids.push_back(bulk.identifier(node));
  }

  size_t num_ids_written = ns->put_field_data("ids", node_ids);
  if ( num_nodes != num_ids_written ) {
    std::ostringstream msg ;
    msg << " FAILED in Ioss::NodeSet::output_node_set:"
        << " num_nodes = " << num_nodes
        << ", num_ids_written = " << num_ids_written;
    throw std::runtime_error( msg.str() );
  }

  stk::mesh::Field<double> *df_field = meta_data.get_field<stk::mesh::Field<double> >("distribution_factors");
  if (df_field != NULL) {
    const stk::mesh::FieldBase::Restriction &res = df_field->restriction(stk::mesh::MetaData::NODE_RANK, *part);
    if (res.dimension() > 0) {
      stk::io::field_data_to_ioss(bulk, df_field, nodes, ns, "distribution_factors", Ioss::Field::MESH);
    }
  }

  const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
  std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const mesh::FieldBase *f = *I ; ++I ;
    const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
    if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
      const mesh::FieldBase::Restriction &res = f->restriction(0, *part);
      if (res.dimension() > 0) {
        stk::io::field_data_to_ioss(bulk, f, nodes, ns, f->name(), Ioss::Field::ATTRIBUTE);
      }
    }
  }
}

template <typename INT>
void output_communication_maps(Ioss::Region &io_region,
                               const stk::mesh::BulkData &bulk,
                               const stk::mesh::Selector *anded_selector,
                               INT /*dummy*/)
{
  if (bulk.parallel_size() > 1) {
    const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);
    mesh::Selector selector = meta.globally_shared_part();
    if (anded_selector) selector &= *anded_selector;

    std::vector<mesh::Entity> entities;
    get_selected_entities(selector, bulk.buckets(stk::mesh::MetaData::NODE_RANK), entities);

    const std::string cs_name("node_symm_comm_spec");
    Ioss::CommSet * io_cs = io_region.get_commset(cs_name);
    ThrowRequire(io_cs != NULL);

    // Allocate data space to store <id, processor> pair
    assert(io_cs->field_exists("entity_processor"));
    int size = io_cs->get_field("entity_processor").raw_count();

    std::vector<INT> ep;
    ep.reserve(size*2);

    for (size_t i=0; i < entities.size(); i++) {
      for ( stk::mesh::PairIterEntityComm ec = bulk.entity_comm(bulk.entity_key(entities[i])); ! ec.empty() ; ++ec ) {
        if (ec->ghost_id == 0) {
          ep.push_back(bulk.identifier(entities[i]));
          ep.push_back(ec->proc);
        }
      }
    }
    io_cs->put_field_data("entity_processor", ep);
  }
}

template <typename INT>
void output_side_set(Ioss::SideSet *ss,
                     const stk::mesh::BulkData &bulk,
                     const stk::mesh::Selector *anded_selector, INT dummy)
{
  const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
  size_t block_count = ss->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *block = ss->get_block(i);
    if (stk::io::include_entity(block)) {
      stk::mesh::Part * const part = meta_data.get_part(block->name());
      stk::io::write_side_data_to_ioss(*block, part, bulk, anded_selector, dummy);
    }
  }
}

} // namespace <blank>

void write_output_db(Ioss::Region& io_region,
                     const stk::mesh::BulkData& bulk,
                     const stk::mesh::Selector *anded_selector)
{
  const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);

  bool ints64bit = db_api_int_size(&io_region) == 8;

  io_region.begin_mode( Ioss::STATE_MODEL );

  int64_t z64 = 0;
  int     z32 = 0;

  Ioss::NodeBlock & nb = *io_region.get_node_blocks()[0];

  if (ints64bit)
	output_node_block(nb, meta.universal_part(), bulk, anded_selector, z64);
  else
	output_node_block(nb, meta.universal_part(), bulk, anded_selector, z32);

  //----------------------------------
  const Ioss::ElementBlockContainer& elem_blocks = io_region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
	  it != elem_blocks.end(); ++it) {
	if (ints64bit)
	  output_element_block(*it, bulk, anded_selector, z64);
	else
	  output_element_block(*it, bulk, anded_selector, z32);
  }

  //----------------------------------
  const Ioss::NodeSetContainer& node_sets = io_region.get_nodesets();
  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
	  it != node_sets.end(); ++it) {
	if (ints64bit)
	  output_node_set(*it, bulk, anded_selector, z64);
	else
	  output_node_set(*it, bulk, anded_selector, z32);
  }

  //----------------------------------
  const Ioss::SideSetContainer& side_sets = io_region.get_sidesets();
  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
	  it != side_sets.end(); ++it) {
	if (ints64bit)
	  output_side_set(*it, bulk, anded_selector, z64);
	else
	  output_side_set(*it, bulk, anded_selector, z32);
  }

  if (ints64bit)
    output_communication_maps(io_region, bulk, anded_selector, z64);
  else
    output_communication_maps(io_region, bulk, anded_selector, z32);

  io_region.end_mode( Ioss::STATE_MODEL );
}

//----------------------------------------------------------------------
bool is_part_io_part(stk::mesh::Part &part)
{
  return NULL != part.attribute<Ioss::GroupingEntity>();
}

// TODO: NOTE: The use of "FieldBase" here basically eliminates the use of the attribute
// for any other fieldbase.  This is just being done now for a proof-of-concept for use
// in the framework-based stk-mesh handoff to a stk-based stk-mesh...  If this looks promising,
// then need to wrap it in a different classs...
const stk::mesh::FieldBase *get_distribution_factor_field(const stk::mesh::Part &p)
{
  return p.attribute<stk::mesh::FieldBase>();
}

void set_distribution_factor_field(stk::mesh::Part &p,
                                   const stk::mesh::FieldBase &df_field)
{
  stk::mesh::MetaData &m = mesh::MetaData::get(p);
  m.declare_attribute_no_delete(p,&df_field);
}

const Ioss::Field::RoleType* get_field_role(const stk::mesh::FieldBase &f)
{
  return f.attribute<Ioss::Field::RoleType>();
}

void set_field_role(stk::mesh::FieldBase &f, const Ioss::Field::RoleType &role)
{
  Ioss::Field::RoleType *my_role = new Ioss::Field::RoleType(role);
  stk::mesh::MetaData &m = mesh::MetaData::get(f);
  const Ioss::Field::RoleType *check = m.declare_attribute_with_delete(f, my_role);
  if ( check != my_role ) {
    if (*check != *my_role) {
      std::ostringstream msg ;
	  msg << " FAILED in IossBridge -- set_field_role:"
	      << " The role type had already been set to " << *check
	      << ", so it is not possible to change it to " << *my_role;
	  throw std::runtime_error( msg.str() );
    }
    delete my_role;
  }
}

}//namespace io
}//namespace stk

