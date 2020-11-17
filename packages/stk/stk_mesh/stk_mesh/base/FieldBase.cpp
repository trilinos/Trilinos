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

#include <stk_mesh/base/FieldBase.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <vector>                       // for vector, etc
#include "Shards_Array.hpp"             // for ArrayDimTag
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg


namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace mesh {

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  s << "Field<" ;
  s << field.data_traits().name ;
  for ( unsigned i = 0 ; i < field.field_array_rank() ; ++i ) {
    s << "," << field.dimension_tags()[i]->name();
  }
  s << ">" ;

  s << "[ name: \"" ;
  s << field.name() ;
  s << "\" , #states: " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  s << b << field << std::endl;
  std::string indent = b;
  indent += "  ";
  print_restrictions(s, indent.c_str(), field);
  return s ;
}

std::ostream & print_restrictions(std::ostream & s ,
                                  const char * const b ,
                                  const FieldBase & field )
{
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();

  for ( const FieldBase::Restriction& r : rMap ) {
    s << b;
    r.print( s, r.selector(), field.field_array_rank() );
    s << std::endl;
  }
  return s;
}

void FieldBase::set_mesh(stk::mesh::BulkData* bulk)
{
  if (m_mesh == NULL || bulk == NULL) {
    m_mesh = bulk;
  }
  else {
    ThrowRequireMsg(bulk == m_mesh, "Internal Error: Trying to use field " << name() << " on more than one bulk data");
  }
}

bool FieldBase::defined_on(const stk::mesh::Part& part) const
{
  return (length(part) > 0);
}

unsigned FieldBase::length(const stk::mesh::Part& part) const
{
  const stk::mesh::FieldRestriction& restriction =
    stk::mesh::find_restriction(*this, entity_rank(), part);
  return restriction.num_scalars_per_entity();
}

void FieldBase::rotate_multistate_data()
{
  const unsigned numStates = m_impl.number_of_states();
  if (numStates > 1 && StateNew == state()) {
    NgpFieldBase* ngpField = get_ngp_field();
    if (ngpField != nullptr) {
      ngpField->rotate_multistate_data();
    }
    for (unsigned s = 1; s < numStates; ++s) {
      FieldBase* sField = field_state(static_cast<FieldState>(s));
      m_field_meta_data.swap(sField->m_field_meta_data);
    }
  }
}

#ifdef STK_DEBUG_FIELD_SYNC
unsigned
FieldBase::get_num_components(const Entity & entity) const
{
  const MeshIndex & index = get_mesh().mesh_index(entity);
  const unsigned bytesPerScalar = data_traits().size_of;
  return m_field_meta_data[index.bucket->bucket_id()].m_bytes_per_entity/bytesPerScalar;
}

void
FieldBase::set_last_modification_view(const LastFieldModLocationType & lastModView) const
{
  m_debugFieldLastModification = lastModView;
}

bool
FieldBase::last_accessed_entity_value_has_changed(const Entity & entity, const uint8_t * lastEntityValues, unsigned component) const
{
  if (get_mesh().is_valid(entity)) {
    const MeshIndex & index = get_mesh().mesh_index(entity);
    const unsigned bytesPerScalar = data_traits().size_of;
    const FieldMetaData& field_meta_data = m_field_meta_data[index.bucket->bucket_id()];
    const unsigned entityOffset = field_meta_data.m_bytes_per_entity * index.bucket_ordinal;
    const unsigned componentOffset = bytesPerScalar * component;

    const uint8_t * currentValue = field_meta_data.m_data + entityOffset + componentOffset;
    const uint8_t * lastValue = lastEntityValues + componentOffset;

    return std::memcmp(currentValue, lastValue, bytesPerScalar);
  }
  else {
    return false;
  }
}

void
FieldBase::set_last_modification(const Entity & entity, unsigned component, LastModLocation location) const
{
  const MeshIndex & index = get_mesh().mesh_index(entity);
  if (get_mesh().in_modifiable_state()) {
    stk::mesh::FieldBase & lastModLocationField = get_last_mod_location_field();
    uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(stk::mesh::ngp_debug_field_data(lastModLocationField,
                                                                                           index.bucket->bucket_id(), index.bucket_ordinal));
    lastModLocation[component] = location;
  }
  else {
    m_debugFieldLastModification(index.bucket->bucket_id(), ORDER_INDICES(index.bucket_ordinal, component)) = location;
  }
}

void
FieldBase::detect_host_field_entity_modification() const
{
  if (!get_mesh().is_valid(m_lastFieldEntity)) return;

  for (unsigned component = 0; component < get_num_components(m_lastFieldEntity);  ++component) {
    if (last_accessed_entity_value_has_changed(m_lastFieldEntity, m_lastFieldValue.data(), component)) {
      set_last_modification(m_lastFieldEntity, component, LastModLocation::HOST);
    }
  }
}

bool
FieldBase::data_is_stale_on_host(const Entity & entity, unsigned component) const
{
  const MeshIndex & index = get_mesh().mesh_index(entity);
  if (get_mesh().in_modifiable_state()) {
    stk::mesh::FieldBase & lastModLocationField = get_last_mod_location_field();
    const uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(stk::mesh::ngp_debug_field_data(lastModLocationField,
                                                                                                 index.bucket->bucket_id(),
                                                                                                 index.bucket_ordinal));
    return !(lastModLocation[component] & LastModLocation::HOST);
  }
  else {
    return !(m_debugFieldLastModification(index.bucket->bucket_id(),
                                          ORDER_INDICES(index.bucket_ordinal, component)) & LastModLocation::HOST);
  }
}

void
FieldBase::store_last_entity_access_location(const Entity & entity) const
{
  const MeshIndex & index = get_mesh().mesh_index(entity);
  const FieldMetaData& field_meta_data = m_field_meta_data[index.bucket->bucket_id()];
  const uint8_t * data = field_meta_data.m_data + field_meta_data.m_bytes_per_entity * index.bucket_ordinal;

  m_lastFieldValue.resize(field_meta_data.m_bytes_per_entity);
  std::memcpy(m_lastFieldValue.data(), data, field_meta_data.m_bytes_per_entity);
  m_lastFieldEntity = entity;
}

std::string location_string(const char * fullPath, int lineNumber)
{
  if (lineNumber != -1) {
    std::string fileName(fullPath);
    std::size_t pathDelimeter = fileName.find_last_of("/");
    if (pathDelimeter < fileName.size()) {
      fileName = fileName.substr(pathDelimeter+1);
    }
    return fileName + ":" + std::to_string(lineNumber) + " ";
  }
  else {
    return "";
  }
}

void
FieldBase::check_stale_field_entity_access(const Entity & entity, const char * fileName, int lineNumber) const
{
  if (get_ngp_field()->lost_device_field_data()) {
    std::cout << location_string(fileName, lineNumber)
              << "*** WARNING: Lost Device values for Field " << name()
              << " due to a mesh modification before a sync to Host" << std::endl;
    return;
  }

  for (unsigned component = 0; component < get_num_components(entity);  ++component) {
    if (data_is_stale_on_host(entity, component)) {
      const MeshIndex & index = get_mesh().mesh_index(entity);
      const unsigned bytesPerScalar = data_traits().size_of;
      const FieldMetaData& field_meta_data = m_field_meta_data[index.bucket->bucket_id()];
      const unsigned entityOffset = field_meta_data.m_bytes_per_entity * index.bucket_ordinal;
      const unsigned componentOffset = bytesPerScalar * component;

      uint8_t * data = field_meta_data.m_data + entityOffset + componentOffset;

      const std::string baseMessage = location_string(fileName, lineNumber) +
                                      "*** WARNING: Accessing stale data on Host for Field " + name() +
                                      "[" + std::to_string(component) + "]=";
      if (data_traits().is_floating_point && data_traits().size_of == 8u) {
        std::cout << baseMessage << *reinterpret_cast<double*>(data) << std::endl;
      }
      else if (data_traits().is_floating_point && data_traits().size_of == 4u) {
        std::cout << baseMessage << *reinterpret_cast<float*>(data) << std::endl;
      }
      else if (data_traits().is_integral && data_traits().is_signed && data_traits().size_of == 4u) {
        std::cout << baseMessage << *reinterpret_cast<int32_t*>(data) << std::endl;
      }
      else if (data_traits().is_integral && data_traits().is_unsigned && data_traits().size_of == 4u) {
        std::cout << baseMessage << *reinterpret_cast<uint32_t*>(data) << std::endl;
      }
      else if (data_traits().is_integral && data_traits().is_signed && data_traits().size_of == 8u) {
        std::cout << baseMessage << *reinterpret_cast<int64_t*>(data) << std::endl;
      }
      else if (data_traits().is_integral && data_traits().is_unsigned && data_traits().size_of == 8u) {
        std::cout << baseMessage << *reinterpret_cast<uint64_t*>(data) << std::endl;
      }
      else if (data_traits().is_integral && data_traits().is_signed && data_traits().size_of == 2u) {
        std::cout << baseMessage << *reinterpret_cast<int16_t*>(data) << std::endl;
      }
      else if (data_traits().is_integral && data_traits().is_unsigned && data_traits().size_of == 2u) {
        std::cout << baseMessage << *reinterpret_cast<uint16_t*>(data) << std::endl;
      }
      else {
        std::cout << location_string(fileName, lineNumber)
                  << "*** WARNING: Accessing stale data on Host for Field " << name() << std::endl;
      }
    }
  }
}

void
FieldBase::store_last_bucket_access_location(const stk::mesh::Bucket & bucket) const
{
  const FieldMetaData& field_meta_data = m_field_meta_data[bucket.bucket_id()];
  m_lastFieldBucketValues.resize(field_meta_data.m_bytes_per_entity * bucket.capacity());

  std::memcpy(m_lastFieldBucketValues.data(), field_meta_data.m_data, field_meta_data.m_bytes_per_entity * bucket.size());
  m_lastFieldBucketEntities.clear();
  for (const stk::mesh::Entity & entity : bucket) {
    m_lastFieldBucketEntities.push_back(entity);
  }
}

void
FieldBase::detect_host_field_bucket_modification() const
{

  for (unsigned ordinal = 0; ordinal < m_lastFieldBucketEntities.size(); ++ordinal) {
    const stk::mesh::Entity & entity = m_lastFieldBucketEntities[ordinal];
    const MeshIndex & index = get_mesh().mesh_index(entity);
    const unsigned bytesPerEntity = m_field_meta_data[index.bucket->bucket_id()].m_bytes_per_entity;
    const uint8_t * lastEntityValues = m_lastFieldBucketValues.data() + ordinal * bytesPerEntity;

    for (unsigned component = 0; component < get_num_components(entity);  ++component) {
      if (last_accessed_entity_value_has_changed(entity, lastEntityValues, component)) {
        set_last_modification(entity, component, LastModLocation::HOST);
      }
    }
  }
}

void
FieldBase::check_stale_field_bucket_access(const stk::mesh::Bucket & bucket, const char * fileName, int lineNumber) const
{
  for (const Entity & entity : bucket) {
    check_stale_field_entity_access(entity, fileName, lineNumber);
  }
}

void
FieldBase::reset_debug_state()
{
  m_lastFieldEntity = stk::mesh::Entity();
  m_lastFieldBucketEntities.clear();
}

FieldBase &
FieldBase::get_last_mod_location_field() const
{
  if (m_lastModLocationField == nullptr) {
    ThrowRequire(get_debug_ngp_field() != nullptr);
    stk::mesh::MetaData & meta = get_mesh().mesh_meta_data();
    meta.enable_late_fields();
    stk::mesh::Field<uint8_t> & lastModLocationField =
        meta.declare_field<stk::mesh::Field<uint8_t>>(entity_rank(),
                                                      "DEBUG_lastFieldModLocation_"+name(),
                                                      number_of_states());

    const RestrictionVector & fieldRestrictions = restrictions();
    for (const Restriction & restriction : fieldRestrictions) {
      const stk::mesh::BucketVector & buckets = get_mesh().get_buckets(entity_rank(), restriction.selector());
      const unsigned numComponents = stk::mesh::field_scalars_per_entity(*this, *buckets.front());
      std::vector<uint8_t> initLastModLocation(numComponents, LastModLocation::HOST_OR_DEVICE);
      stk::mesh::put_field_on_mesh(lastModLocationField, restriction.selector(), numComponents, initLastModLocation.data());
    }

    m_lastModLocationField = &lastModLocationField;
  }
  return *m_lastModLocationField;
}

void
FieldBase::fill_last_mod_location_field_from_device() const
{
  FieldBase & lastModLocationField = get_last_mod_location_field();
  stk::mesh::NgpFieldBase & ngpField = *get_ngp_field();

  const stk::mesh::BucketVector & buckets = get_mesh().buckets(entity_rank());
  for (const stk::mesh::Bucket * bucket : buckets) {
    for (unsigned ordinal = 0; ordinal < bucket->size(); ++ordinal) {
      const stk::mesh::Entity & entity = (*bucket)[ordinal];
      const unsigned numComponents = stk::mesh::field_scalars_per_entity(lastModLocationField, entity);
      uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(stk::mesh::ngp_debug_field_data(lastModLocationField, entity));
      for (unsigned component = 0; component < numComponents; ++component) {
        const unsigned bucketOffset = ngpField.get_bucket_offset(bucket->bucket_id());
        lastModLocation[component] = m_debugFieldLastModification(bucketOffset, ORDER_INDICES(ordinal, component));
      }
    }
  }
}

void
FieldBase::fill_last_mod_location_view_from_host() const
{
  FieldBase & lastModLocationField = get_last_mod_location_field();
  stk::mesh::NgpFieldBase & ngpField = *get_ngp_field();

  const stk::mesh::BucketVector & buckets = get_mesh().buckets(entity_rank());
  for (const stk::mesh::Bucket * bucket : buckets) {
    for (unsigned ordinal = 0; ordinal < bucket->size(); ++ordinal) {
      const stk::mesh::Entity & entity = (*bucket)[ordinal];
      uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(stk::mesh::ngp_debug_field_data(lastModLocationField, entity));
      const unsigned bucketOffset = ngpField.get_bucket_offset(bucket->bucket_id());
      const unsigned numComponents = stk::mesh::field_scalars_per_entity(lastModLocationField, entity);

      for (unsigned component = 0; component < numComponents; ++component) {
        m_debugFieldLastModification(bucketOffset, ORDER_INDICES(ordinal, component)) =
            static_cast<LastModLocation>(lastModLocation[component]);
      }
    }
  }
}

#endif

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

