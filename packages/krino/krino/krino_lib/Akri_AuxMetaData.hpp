// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_AuxMetaData_h
#define Akri_AuxMetaData_h

#include <functional>
#include <string>                       // for operator<<, basic_string, etc

#include "Akri_FieldRef.hpp"
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Entity.hpp"   // for Entity

#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 9)) || (__GNUC__ == 5)
// Looks like there is an issue with these compilers
// related to instatiating std::function with an incomplete 
// type leading to undefine behavior
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/FieldBase.hpp"  // for FieldState
#else
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class FieldBase; } }
#endif

namespace krino
{

class FieldType {
public:
  static const FieldType UNSIGNED_INTEGER_64;
  static const FieldType UNSIGNED_INTEGER;
  static const FieldType INTEGER;
  static const FieldType REAL;

  static const FieldType VECTOR_2D;
  static const FieldType VECTOR_3D;
  static const FieldType MATRIX_22;
  static const FieldType MATRIX_33;

  const std::string & name() const { return my_name; }
  const std::type_info & type_info() const { return my_type_info; }
  unsigned dimension() const { return my_dimension; }

  FieldType(const std::string & nme, const std::type_info & type, const unsigned dim) : my_name(nme), my_type_info(type), my_dimension(dim) {}

private:
  const std::string my_name;
  const std::type_info & my_type_info;
  const unsigned my_dimension;
};

class AuxMetaData final
{
public:
  static AuxMetaData & get(const stk::mesh::MetaData & stk_meta);
  static bool has(const stk::mesh::MetaData & stk_meta);
  static AuxMetaData & create(stk::mesh::MetaData & stk_meta); // must be called before calling get
  static AuxMetaData & get_or_create(stk::mesh::MetaData & stk_meta);

  AuxMetaData ( const AuxMetaData & ) = delete;
  AuxMetaData & operator= ( const AuxMetaData & ) = delete;

  bool using_fmwk() const { return is_fmwk; }

  stk::mesh::Part & active_part() const { return *my_active_part; }
  stk::mesh::Part & block_boundary_part() const;
  stk::mesh::Part & exposed_boundary_part() const;

  const stk::mesh::Selector & active_not_ghost_selector() const { return my_active_not_ghost_selector; }
  const stk::mesh::Selector & active_locally_owned_selector() const { return my_active_locally_owned_selector; }
  const stk::mesh::Selector & active_globally_shared_selector() const { return my_active_globally_shared_selector; }

  bool get_force_64bit_flag() const { return my_force_64bit_flag; }
  void clear_force_64bit_flag() { my_force_64bit_flag = false; }
  bool get_assert_32bit_flag() { return my_assert_32bit_flag; }
  void set_assert_32bit_flag() { my_force_64bit_flag = false; my_assert_32bit_flag = true; }

  FieldRef declare_field(
    const std::string & fld_name,
    const FieldType & field_type,
    const stk::mesh::EntityRank entity_rank,
    const unsigned num_states);

  FieldRef register_field(
    const std::string & fld_name,
    const FieldType & field_type,
    const stk::mesh::EntityRank entity_rank,
    const unsigned num_states,
    const unsigned dimension,
    const stk::mesh::Part & part,
    const void * value_type_init = nullptr);

  void assign_part_id(stk::mesh::Part& part);
  stk::mesh::Part & declare_io_part(const std::string & name, stk::mesh::EntityRank entityRank, const bool restartOnlyIOPart=false);
  stk::mesh::Part & declare_io_part_with_topology(const std::string & name, const stk::topology topology, const bool restartOnlyIOPart=false);

  const std::vector<stk::mesh::Part*> & get_restart_only_io_parts() const { return myRestartOnlyIOParts; }

  bool has_part( const std::string& name ) const;
  stk::mesh::Part& get_part( const std::string& name ) const;

  void define_part_alias( stk::mesh::Part & part, const std::string & alias );

  bool has_field( const stk::mesh::EntityRank obj_type, const std::string& name ) const;
  bool has_field( const stk::mesh::EntityRank obj_type, const std::string& name, stk::mesh::FieldState state ) const;

  FieldRef get_field( const stk::mesh::EntityRank obj_type, const std::string& name ) const;
  FieldRef get_field( const stk::mesh::EntityRank obj_type, const std::string& name, stk::mesh::FieldState state ) const;

  void induce_topology_nodesets(stk::mesh::Selector selector = stk::mesh::Selector()) const;
  stk::topology get_nodal_field_topology( const stk::mesh::FieldBase & field, stk::mesh::Entity entity ) const;
  stk::topology get_nodal_field_topology( const stk::mesh::FieldBase & field, const stk::mesh::Bucket & bucket ) const;
  stk::mesh::Selector selectField( const stk::mesh::FieldBase & field, const stk::mesh::EntityRank target_rank ) const;
  bool is_cell_edge(stk::mesh::Entity node0, stk::mesh::Entity node1) const { return fn_is_cell_edge(node0, node1); }

  void set_fmwk_functions(
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
      stk::mesh::Part * in_block_boundary_part);
  void set_inducer_functions(
      const std::function<void(stk::mesh::Selector selector)> & in_inducer_induce_topology_nodesets,
      const std::function<stk::topology( const stk::mesh::FieldBase & field, const stk::mesh::Bucket & bucket )> & in_inducer_get_nodal_field_topology,
      const std::function<stk::mesh::Selector( const stk::mesh::FieldBase & field, const stk::mesh::EntityRank target_rank )> & in_inducer_selectField);
  void set_is_cell_edge_function(const std::function<bool(stk::mesh::Entity node0, stk::mesh::Entity node1)> & in_is_cell_edge) { fn_is_cell_edge = in_is_cell_edge; }
  FieldRef get_current_coordinates() const;
  void set_current_coordinates(FieldRef current_coords) { my_current_coordinates = current_coords; }

private:
  explicit AuxMetaData(stk::mesh::MetaData & stk_meta);
private:
  stk::mesh::MetaData & my_meta;
  bool is_fmwk;
  bool my_assert_32bit_flag;
  bool my_force_64bit_flag;
  stk::mesh::Part * my_active_part;
  stk::mesh::Part * my_exposed_boundary_part;
  stk::mesh::Part * my_block_boundary_part;
  stk::mesh::Selector my_active_not_ghost_selector;
  stk::mesh::Selector my_active_locally_owned_selector;
  stk::mesh::Selector my_active_globally_shared_selector;
  std::function<stk::mesh::Part *(const std::string & name)> fmwk_get_iopart;
  std::function<stk::mesh::Part &(const std::string & name, const stk::mesh::EntityRank entityRank)> fmwk_iopart;
  std::function<void(stk::mesh::Part & part, const std::string & alias)> fmwk_define_iopart_alias;
  std::function<stk::mesh::FieldBase &
    (const std::string & name,
     const std::string & field_type_name,
     const std::type_info & field_value_type,
     const unsigned field_value_dimension,
     const stk::mesh::EntityRank rank,
     const unsigned num_states,
     const unsigned dimension,
     const stk::mesh::Part & part,
     const void * value_type_init)> fmwk_register_field;
  std::function<void(stk::mesh::Selector selector)> fn_induce_topology_nodesets;
  std::function<stk::topology( const stk::mesh::FieldBase & field, const stk::mesh::Bucket & bucket )> fn_get_nodal_field_topology;
  std::function<stk::mesh::Selector( const stk::mesh::FieldBase & field, const stk::mesh::EntityRank target_rank )> fn_selectField;
  std::function<bool(stk::mesh::Entity node0, stk::mesh::Entity node1)> fn_is_cell_edge;
  std::vector<stk::mesh::Part *> myRestartOnlyIOParts;
  mutable FieldRef my_current_coordinates;
};

} // namespace krino

#endif // Akri_AuxMetaData_h
