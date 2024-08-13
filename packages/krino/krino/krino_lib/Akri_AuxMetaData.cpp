// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AuxMetaData.hpp>
#include "Akri_FieldRef.hpp"              // for FieldRef

#include "stk_mesh/base/FieldBase.hpp"  // for FieldState
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_io/IossBridge.hpp>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc

namespace krino
{
const FieldType FieldType::UNSIGNED_INTEGER_64("UNSIGNED INTEGER 64", typeid(uint64_t), 1);
const FieldType FieldType::UNSIGNED_INTEGER("UNSIGNED INTEGER", typeid(unsigned), 1);
const FieldType FieldType::INTEGER("INTEGER", typeid(int), 1);
const FieldType FieldType::REAL("REAL", typeid(double), 1);

const FieldType FieldType::VECTOR_2D("VECTOR_2D", typeid(double), 2);
const FieldType FieldType::VECTOR_3D("VECTOR_3D", typeid(double), 3);
const FieldType FieldType::MATRIX_22("MATRIX_22", typeid(double), 4);
const FieldType FieldType::MATRIX_33("MATRIX_33", typeid(double), 9);

AuxMetaData &
AuxMetaData::get(const stk::mesh::MetaData & stk_meta)
{
  AuxMetaData * aux_meta = const_cast<AuxMetaData*>(stk_meta.get_attribute<AuxMetaData>());
  STK_ThrowRequireMsg(nullptr != aux_meta, "AuxMetaData not found on MetaData.");
  return *aux_meta;
}

bool AuxMetaData::has(const stk::mesh::MetaData & stk_meta)
{
  AuxMetaData * aux_meta = const_cast<AuxMetaData*>(stk_meta.get_attribute<AuxMetaData>());
  return aux_meta != nullptr;
}

AuxMetaData &
AuxMetaData::create(stk::mesh::MetaData & stk_meta)
{
  AuxMetaData * aux_meta = const_cast<AuxMetaData*>(stk_meta.get_attribute<AuxMetaData>());
  STK_ThrowRequireMsg(nullptr == aux_meta, "AuxMetaData::create should be caled only once per MetaData.");
  if (nullptr == aux_meta)
  {
    aux_meta = new AuxMetaData(stk_meta);
    stk_meta.declare_attribute_with_delete<AuxMetaData>(aux_meta);
  }
  return *aux_meta;
}

AuxMetaData & 
AuxMetaData::get_or_create(stk::mesh::MetaData & stk_meta)
{
  if (AuxMetaData::has(stk_meta)) return AuxMetaData::get(stk_meta);

  return AuxMetaData::create(stk_meta);
}

AuxMetaData::AuxMetaData(stk::mesh::MetaData & stk_meta)
  : my_meta(stk_meta),
    is_fmwk(false),
    my_assert_32bit_flag(false),
    my_force_64bit_flag(true), // just to guarantee testing 64 bit as much as possible
    my_active_part(nullptr),
    my_exposed_boundary_part(nullptr),
    my_block_boundary_part(nullptr)
{
  stk::mesh::Part * existing_active_part = my_meta.get_part("ACTIVE_CONTEXT_BIT");
  if (nullptr != existing_active_part)
  {
    my_active_part = existing_active_part;
    my_exposed_boundary_part = my_meta.get_part("EXPOSED_BOUNDARY_CONTEXT_BIT");
    my_block_boundary_part = my_meta.get_part("BLOCK_BOUNDARY_CONTEXT_BIT");
  }
  else
  {
    STK_ThrowRequireMsg(my_meta.spatial_dimension() > 0, "For non-Fmwk usage, AuxMetaData cannot be created until after the spatial_dimension is set on the stk::mesh::MetaData.");
    STK_ThrowRequireMsg(!my_meta.is_commit(), "For non-Fmwk usage, AuxMetaData must be created before the stk::mesh::MetaData is committed.");
    my_active_part           = &my_meta.declare_part("ACTIVE_CONTEXT_BIT");
    my_exposed_boundary_part = &my_meta.declare_part("EXPOSED_BOUNDARY_CONTEXT_BIT", my_meta.side_rank());
    my_block_boundary_part   = &my_meta.declare_part("BLOCK_BOUNDARY_CONTEXT_BIT", my_meta.side_rank());
  }

  my_active_globally_shared_selector = stk::mesh::Selector(*my_active_part & my_meta.globally_shared_part());
  my_active_not_ghost_selector = stk::mesh::Selector(*my_active_part & (my_meta.locally_owned_part() | my_meta.globally_shared_part()));
  my_active_locally_owned_selector = stk::mesh::Selector(*my_active_part & my_meta.locally_owned_part());
}

void
AuxMetaData::set_fmwk_functions(
  const std::function<stk::mesh::Part *(const std::string & name)> & in_fmwk_get_iopart,
  const std::function<stk::mesh::Part &(const std::string & name, const stk::mesh::EntityRank entityRank)> & in_fmwk_iopart,
  const std::function<void(stk::mesh::Part & part, const std::string & alias)> & in_fmwk_define_iopart_alias,
  const std::function<stk::mesh::FieldBase &
    (const std::string & name,
     const std::string & field_type_name,
     const std::type_info & field_value_type,
     const unsigned field_value_dimension,
     const stk::mesh::EntityRank rank,
     const unsigned num_states,
     const unsigned dimension,
     const stk::mesh::Part & part,
     const void * value_type_init)> & in_fmwk_register_field,
  stk::mesh::Part * in_exposed_boundary_part,
  stk::mesh::Part * in_block_boundary_part)
{
  is_fmwk = true;
  clear_force_64bit_flag();
  fmwk_get_iopart = in_fmwk_get_iopart;
  fmwk_iopart = in_fmwk_iopart;
  fmwk_define_iopart_alias = in_fmwk_define_iopart_alias;
  fmwk_register_field = in_fmwk_register_field;
  my_exposed_boundary_part = in_exposed_boundary_part;
  my_block_boundary_part = in_block_boundary_part;
}

void
AuxMetaData::set_inducer_functions(
  const std::function<void(stk::mesh::Selector selector)> & in_inducer_induce_topology_nodesets,
  const std::function<stk::topology( const stk::mesh::FieldBase & field, const stk::mesh::Bucket & bucket )> & in_inducer_get_nodal_field_topology,
  const std::function<stk::mesh::Selector( const stk::mesh::FieldBase & field, const stk::mesh::EntityRank target_rank )> & in_inducer_selectField)
{
  fn_induce_topology_nodesets = in_inducer_induce_topology_nodesets;
  fn_get_nodal_field_topology = in_inducer_get_nodal_field_topology;
  fn_selectField = in_inducer_selectField;
}

stk::mesh::Part & AuxMetaData::block_boundary_part() const
{
  STK_ThrowAssert(nullptr != my_block_boundary_part); return *my_block_boundary_part;
}

stk::mesh::Part & AuxMetaData::exposed_boundary_part() const
{
  STK_ThrowAssert(nullptr != my_exposed_boundary_part); return *my_exposed_boundary_part;
}

bool
AuxMetaData::has_field( const stk::mesh::EntityRank obj_type, const std::string& name ) const
{
  return NULL != my_meta.get_field(obj_type, name);
}

bool
AuxMetaData::has_field( const stk::mesh::EntityRank obj_type, const std::string& name, stk::mesh::FieldState state ) const
{
  stk::mesh::FieldBase* field_ptr = my_meta.get_field(obj_type, name);
  if (NULL == field_ptr) return false;
  return NULL != field_ptr->field_state(state);
}

FieldRef
AuxMetaData::get_field( const stk::mesh::EntityRank obj_type, const std::string& name ) const
{
  stk::mesh::FieldBase* field_ptr = my_meta.get_field(obj_type, name);
  STK_ThrowRequireMsg(NULL != field_ptr, "Field \"" << name << "\" not found.");
  return FieldRef(field_ptr);
}


FieldRef
AuxMetaData::get_field( const stk::mesh::EntityRank obj_type, const std::string& name, stk::mesh::FieldState state ) const
{
  stk::mesh::FieldBase* field_ptr = my_meta.get_field(obj_type, name);
  if(field_ptr == nullptr)
  {
    std::ostringstream err_msg;
    err_msg << "Field \"" << name << "\" not found.\n";
    err_msg << "Registered fields are:\n";
    for(auto && field : my_meta.get_fields(obj_type))
    {
      err_msg << field->name() << "\n";
    }
    throw std::runtime_error(err_msg.str());
  }
  return FieldRef(field_ptr, state);
}

bool
AuxMetaData::has_part( const std::string& name ) const
{
  stk::mesh::Part * part = my_meta.get_part(name);
  // If the part is not found see if the name is actually a Fmwk alias
  if (!part && is_fmwk) part = fmwk_get_iopart(name);
  return NULL != part;
}

stk::mesh::Part&
AuxMetaData::get_part( const std::string& name ) const
{
  stk::mesh::Part * part = my_meta.get_part(name);
  // If the part is not found see if the name is actually a Fmwk alias
  if (!part && is_fmwk) part = fmwk_get_iopart(name);
  STK_ThrowRequireMsg(part, "Could not find part " << name;);
  return *part;
}

void
AuxMetaData::define_part_alias( stk::mesh::Part & part , const std::string & alias )
{
  if (is_fmwk)
  {
    fmwk_define_iopart_alias(part, alias);
    return;
  }
  // Silent no-op
}

void
AuxMetaData::assign_part_id(stk::mesh::Part& part)
{
  static int64_t max_existing_part_id = stk::mesh::Part::INVALID_ID;
  if (max_existing_part_id == stk::mesh::Part::INVALID_ID)
  {
    max_existing_part_id = 0;
    for (auto && existing_part : my_meta.get_parts())
    {
      max_existing_part_id = std::max(max_existing_part_id, existing_part->id());
    }
  }
  if (part.id() == stk::mesh::Part::INVALID_ID)
  {
    const int64_t part_id = ++max_existing_part_id;
    my_meta.set_part_id(part, part_id);
  }
}

stk::mesh::Part&
AuxMetaData::declare_io_part( const std::string& name, const stk::mesh::EntityRank obj_type, const bool restartOnlyIOPart )
{
  if (is_fmwk)
  {
    stk::mesh::Part & part = fmwk_iopart(name, obj_type);
    if (restartOnlyIOPart)
    {
      myRestartOnlyIOParts.push_back(&part);
    }
    assign_part_id(part);
    return part;
  }

  stk::mesh::Part & part = my_meta.declare_part(name, obj_type);

  if (restartOnlyIOPart)
  {
    myRestartOnlyIOParts.push_back(&part);
  }
  else
  {
    if (!stk::io::is_part_io_part(part))
    {
      stk::io::put_io_part_attribute(part);
    }
  }

  assign_part_id(part);

  return part;
}

stk::mesh::Part&
AuxMetaData::declare_io_part_with_topology( const std::string& name, const stk::topology topology, const bool restartOnlyIOPart )
{
  if (is_fmwk)
  {
    stk::mesh::Part & root_part = my_meta.get_topology_root_part(topology);
    stk::mesh::Part & part = fmwk_iopart(name, root_part.primary_entity_rank());
    my_meta.declare_part_subset(root_part, part);
    if (restartOnlyIOPart)
    {
      myRestartOnlyIOParts.push_back(&part);
    }
    assign_part_id(part);
    return part;
  }

  stk::mesh::Part & part = my_meta.declare_part_with_topology(name, topology);

  if (restartOnlyIOPart)
  {
    myRestartOnlyIOParts.push_back(&part);
  }
  else
  {
    if (!stk::io::is_part_io_part(part))
    {
      stk::io::put_io_part_attribute(part);
    }
  }

  assign_part_id(part);

  return part;
}

FieldRef
AuxMetaData::declare_field(
    const std::string & fld_name,
    const FieldType & field_type,
    const stk::mesh::EntityRank entity_rank,
    const unsigned num_states)
{
  stk::mesh::FieldBase * field = NULL;
  const std::type_info & value_type = field_type.type_info();
  if (value_type == typeid(int))
    field = &my_meta.declare_field<int>(entity_rank, fld_name, num_states);
  else if (value_type == typeid(double))
    field = &my_meta.declare_field<double>(entity_rank, fld_name, num_states);
  else if (value_type == typeid(unsigned))
    field = &my_meta.declare_field<unsigned>(entity_rank, fld_name, num_states);
  else if (value_type == typeid(int64_t))
    field = &my_meta.declare_field<int64_t>(entity_rank, fld_name, num_states);
  else if (value_type == typeid(uint64_t))
    field = &my_meta.declare_field<uint64_t>(entity_rank, fld_name, num_states);
  else {
    STK_ThrowRequireMsg(false, "Unhandled primitive type " << value_type.name());
  }

  return FieldRef(field);
}

FieldRef
AuxMetaData::register_field(
    const std::string & fld_name,
    const FieldType & field_type,
    const stk::mesh::EntityRank entity_rank,
    const unsigned num_states,
    const unsigned dimension,
    const stk::mesh::Part & part,
    const void * value_type_init)
{
  // don't rely on induction of non-ranked superset parts
  // Note that we are only checking nodal fields here, but this could also be a problem if we are counting on a non-ranked superset part to be induced on faces or edges (unlikely?).
  //ThrowRequire(entity_rank != stk::topology::NODE_RANK || !is_unranked_superset_part(part));

  if (is_fmwk)
  {
    return FieldRef(fmwk_register_field(fld_name, field_type.name(), field_type.type_info(), field_type.dimension(), entity_rank, num_states, dimension, part, value_type_init));
  }

  if (field_type.name() == FieldType::VECTOR_2D.name())
  {
    auto & field = my_meta.declare_field<double>(entity_rank, fld_name, num_states);
    stk::mesh::put_field_on_mesh(field, part, field_type.dimension(), dimension, nullptr);
    stk::io::set_field_output_type(field, stk::io::FieldOutputType::VECTOR_2D);
    return FieldRef(field);
  }
  else if (field_type.name() == FieldType::VECTOR_3D.name())
  {
    auto & field = my_meta.declare_field<double>(entity_rank, fld_name, num_states);
    stk::mesh::put_field_on_mesh(field, part, field_type.dimension(), dimension, nullptr);
    stk::io::set_field_output_type(field, stk::io::FieldOutputType::VECTOR_3D);
    return FieldRef(field);
  }

  FieldRef field = declare_field(fld_name, field_type, entity_rank, num_states);
  stk::mesh::put_field_on_mesh(field.field(), part, field_type.dimension(), dimension, value_type_init);
  return field;
}

void
AuxMetaData::induce_topology_nodesets(stk::mesh::Selector selector) const
{
  if (fn_induce_topology_nodesets)
  {
    fn_induce_topology_nodesets(selector);
  }
  // no-op if inducer is not set
  return;
}

stk::topology
AuxMetaData::get_nodal_field_topology( const stk::mesh::FieldBase & field, stk::mesh::Entity entity ) const
{
  return get_nodal_field_topology(field, my_meta.mesh_bulk_data().bucket(entity));
}

stk::topology AuxMetaData::get_nodal_field_topology( const stk::mesh::FieldBase & field, const stk::mesh::Bucket & bucket ) const
{
  if (fn_get_nodal_field_topology)
  {
    return fn_get_nodal_field_topology(field, bucket);
  }
  // return bucket topology if inducer is not set and field is defined on bucket
  const stk::mesh::Selector field_selector = stk::mesh::selectField(field);
  stk::topology field_topology = stk::topology::INVALID_TOPOLOGY;
  if (field_selector(bucket)) field_topology = bucket.topology();
  return field_topology;
}

stk::mesh::Selector AuxMetaData::selectField( const stk::mesh::FieldBase & field, const stk::mesh::EntityRank target_rank ) const
{
  if (fn_selectField && field.entity_rank() == stk::topology::NODE_RANK)
  {
    return fn_selectField( field, target_rank );
  }
  // return stk::mesh::selectField if inducer is not set
  return stk::mesh::selectField(field);
}

FieldRef AuxMetaData::get_current_coordinates() const
{
  if (!my_current_coordinates.valid())
  {
    const stk::mesh::FieldBase * meta_coords = my_meta.coordinate_field();
    STK_ThrowRequireMsg(nullptr != meta_coords, "Coordinates must be defined before calling AuxMetaData::get_current_coordinates().");
    my_current_coordinates = FieldRef(meta_coords);
  }
  return my_current_coordinates;
}

//----------------------------------------------------------------------

} // namespace krino
