// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_FieldRef_h
#define Akri_FieldRef_h

#include <stddef.h>                     // for NULL
#include <string>                       // for string
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"  // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc
#include "stk_mesh/base/FieldState.hpp"  // for FieldState, etc
#include "stk_mesh/base/Types.hpp"      // for EntityRank
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, etc
namespace stk { namespace mesh { class Bucket; } }

namespace krino{

std::string
state_string(stk::mesh::FieldState state);

class FieldRef
{
public:
  // Only default constructed FieldRef (or copies of a default constructed FieldRef) should ever have my_field==NULL
  FieldRef() : my_field(NULL) {}

  FieldRef(const FieldRef &) = default;
  FieldRef & operator= (const FieldRef &) = default;

  FieldRef(const stk::mesh::FieldBase * in_field)
  : my_field(in_field)
  { 
    STK_ThrowAssertMsg(NULL != in_field, "Cannot set FieldRef with NULL field.");
  }
  FieldRef(const stk::mesh::FieldBase * in_field, const stk::mesh::FieldState in_state)
  : my_field(in_field)
  { 
    STK_ThrowAssertMsg(NULL != in_field, "Cannot set FieldRef with NULL field.");
    my_field = in_field->field_state(in_state);
    STK_ThrowAssertMsg(NULL != my_field, "Invalid state.");
  }
  FieldRef(const stk::mesh::FieldBase & in_fieldref)
  {
    my_field = &in_fieldref;
  }
  FieldRef(const FieldRef in_fieldref, const stk::mesh::FieldState in_state)
  { 
    STK_ThrowAssertMsg(NULL != in_fieldref.my_field, "Cannot set FieldRef with NULL field.");
    my_field = in_fieldref.my_field->field_state(in_state);
  }

  FieldRef field_state(stk::mesh::FieldState in_state) const
  {
    return FieldRef(my_field, in_state);
  }

  // assignment operators
  FieldRef& operator=( const stk::mesh::FieldBase * in_field )
  { 
    STK_ThrowAssertMsg(NULL != in_field, "Cannot set FieldRef with NULL field.");
    my_field = in_field;
    return *this;
  }
  FieldRef& operator=( const stk::mesh::FieldBase & in_field )
  {
    my_field = &in_field;
    return *this;
  }

  bool valid() const { return NULL != my_field; }
  bool valid_restriction_size() const { return NULL != my_field && 0 != my_field->restrictions().size(); }
  void assert_valid() const { STK_ThrowAssertMsg(NULL != my_field, "Attempt to access field of uninitialized FieldRef."); }

  // stk access
  const stk::mesh::FieldBase & field() const { assert_valid(); return *my_field; }
  stk::mesh::FieldBase & field() { assert_valid(); return const_cast<stk::mesh::FieldBase &>(*my_field); }
  operator const stk::mesh::FieldBase &() const { assert_valid(); return *my_field; }

  // pass through methods
  unsigned number_of_states() const { assert_valid(); return my_field->number_of_states(); }
  stk::mesh::FieldState state() const { assert_valid(); return my_field->state(); }
  stk::mesh::EntityRank entity_rank() const { assert_valid(); return my_field->entity_rank(); }

  template<class Type> bool type_is() const
  { return my_field->type_is<Type>(); }

  // name_with_state() includes any suffix indicating the state, ie "_STKFS_OLD"
  const std::string & name_with_state() const { assert_valid(); return my_field->name(); }
  // name() is the basic field name without any suffix for state, regardless of the state of the field
  const std::string & name() const { assert_valid(); return my_field->field_state(stk::mesh::StateNone)->name(); }

  unsigned length(const stk::mesh::Bucket& b) const
  {
    assert_valid();
    const unsigned type_stride = my_field->data_traits().stride_of;
    const unsigned fieldSize = stk::mesh::field_bytes_per_entity(*my_field, b) / type_stride;
    return fieldSize;
  }

  unsigned length(stk::mesh::Entity e) const
  {
    assert_valid();
    const unsigned type_stride = my_field->data_traits().stride_of;
    const unsigned fieldSize = stk::mesh::field_bytes_per_entity(*my_field, e) / type_stride;
    return fieldSize;
  }

  // Only safe if field has same length on all entities
  unsigned length() const
  {
    STK_ThrowAssert(valid());
    return my_field->max_size();
  }

  void sync_to_host() const
  {
    my_field->sync_to_host();
  }

  // testing
  bool operator ==(const FieldRef & rhs) const { return my_field == rhs.my_field; }
  bool operator !=(const FieldRef & rhs) const { return my_field != rhs.my_field; }
  bool operator < (const FieldRef & rhs) const;

private:
  const stk::mesh::FieldBase * my_field;
};

typedef std::set<FieldRef> FieldSet;

inline std::ostream & operator << (std::ostream & os, const FieldRef & field_ref)
{
  const bool valid = field_ref.valid();
  os << "FieldRef valid = " << std::boolalpha << valid;
  if(valid)
  {
    os << " field = " << field_ref.name_with_state();
  }
  return os;
}

template< class pod_type >
pod_type * field_data( const stk::mesh::FieldBase & field, stk::mesh::Entity entity )
{
  return stk::mesh::field_data(static_cast< const stk::mesh::Field<pod_type>& >(field), entity);
}

template< class pod_type >
pod_type * field_data( const stk::mesh::FieldBase & field, const stk::mesh::Bucket & bucket )
{
  return stk::mesh::field_data(static_cast< const stk::mesh::Field<pod_type>& >(field), bucket);
}

inline bool has_field_data( const FieldRef fieldref, stk::mesh::Entity entity )
{
  return stk::mesh::field_is_allocated_for_bucket(fieldref.field(), fieldref.field().get_mesh().bucket(entity));
}

inline bool has_field_data( const FieldRef fieldref, const stk::mesh::Bucket & bucket )
{
  return stk::mesh::field_is_allocated_for_bucket(fieldref.field(), bucket);
}

} // namespace krino

#endif // Akri_FieldRef_h
