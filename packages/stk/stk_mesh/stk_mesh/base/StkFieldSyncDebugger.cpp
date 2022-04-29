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

#include "StkFieldSyncDebugger.hpp"
#include "FieldBase.hpp"
#include "BulkData.hpp"
#include "MetaData.hpp"
#include "FieldRestriction.hpp"
#include "stk_mesh/baseImpl/BucketRepository.hpp"

namespace stk {
namespace mesh {

StkFieldSyncDebugger::StkFieldSyncDebugger(const FieldBase* stkField)
  : m_stkField(stkField),
    m_isDataInitialized(false)
{
}

void
StkFieldSyncDebugger::host_stale_access_entity_check(const stk::mesh::Entity& entity, const char* fileName, int lineNumber)
{
  if (m_isDataInitialized) {
    impl::get_ngp_field(*m_stkField)->debug_detect_device_field_modification();
    detect_host_field_entity_modification();
    check_stale_field_entity_access(entity, fileName, lineNumber);
    store_last_entity_access_location(entity);
  }
}

void
StkFieldSyncDebugger::host_stale_access_entity_check(const unsigned & bucketId, const unsigned & bucketOrd,
                                                     const char * fileName, int lineNumber)
{
  BulkData & bulk = m_stkField->get_mesh();
  const Bucket * bucket = bulk.get_bucket_repository().get_bucket(m_stkField->entity_rank(), bucketId);
  ThrowRequire(bucket != nullptr);

  host_stale_access_entity_check((*bucket)[bucketOrd], fileName, lineNumber);
}

void
StkFieldSyncDebugger::host_stale_access_bucket_check(const stk::mesh::Bucket& bucket, const char* fileName, int lineNumber)
{
  if (m_isDataInitialized) {
    impl::get_ngp_field(*m_stkField)->debug_detect_device_field_modification();
    detect_host_field_bucket_modification();
    check_stale_field_bucket_access(bucket, fileName, lineNumber);
    store_last_bucket_access_location(bucket);
  }
}

void
StkFieldSyncDebugger::host_stale_access_bucket_check(const unsigned& bucketId, const char* fileName, int lineNumber)
{
  BulkData & bulk = m_stkField->get_mesh();
  const Bucket * bucket = bulk.get_bucket_repository().get_bucket(m_stkField->entity_rank(), bucketId);
  ThrowRequire(bucket != nullptr);

  host_stale_access_bucket_check(*bucket, fileName, lineNumber);
}

void
StkFieldSyncDebugger::clear_last_field_value()
{
  m_lastFieldEntity = stk::mesh::Entity();
  m_lastFieldBucketEntities.clear();
}

void
StkFieldSyncDebugger::fill_last_mod_location_view_from_host()
{
  FieldBase & lastModLocationField = get_last_mod_location_field();
  stk::mesh::NgpFieldBase & ngpField = *impl::get_ngp_field(*m_stkField);
  BulkData & bulk = m_stkField->get_mesh();

  const stk::mesh::BucketVector & buckets = bulk.buckets(m_stkField->entity_rank());
  for (const stk::mesh::Bucket * bucket : buckets) {
    for (unsigned ordinal = 0; ordinal < bucket->size(); ++ordinal) {
      const stk::mesh::Entity & entity = (*bucket)[ordinal];
      uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(field_data<FieldBase, EmptyStkFieldSyncDebugger>(lastModLocationField, entity));
      const unsigned bucketOffset = ngpField.debug_get_bucket_offset(bucket->bucket_id());
      const unsigned numComponents = stk::mesh::field_scalars_per_entity(lastModLocationField, entity);

      for (unsigned component = 0; component < numComponents; ++component) {
        m_debugFieldLastModification(bucketOffset, ORDER_INDICES(ordinal, component)) =
            static_cast<LastModLocation>(lastModLocation[component]);
      }
    }
  }
}

void
StkFieldSyncDebugger::fill_last_mod_location_field_from_device()
{
  FieldBase & lastModLocationField = get_last_mod_location_field();
  NgpFieldBase & ngpField = *impl::get_ngp_field(*m_stkField);
  BulkData & bulk = m_stkField->get_mesh();

  const BucketVector & buckets = bulk.buckets(m_stkField->entity_rank());
  for (const Bucket * bucket : buckets) {
    for (unsigned ordinal = 0; ordinal < bucket->size(); ++ordinal) {
      const Entity & entity = (*bucket)[ordinal];
      const unsigned numComponents = field_scalars_per_entity(lastModLocationField, entity);
      uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(field_data<FieldBase, EmptyStkFieldSyncDebugger>(lastModLocationField, entity));
      for (unsigned component = 0; component < numComponents; ++component) {
        const unsigned bucketOffset = ngpField.debug_get_bucket_offset(bucket->bucket_id());
        lastModLocation[component] = m_debugFieldLastModification(bucketOffset, ORDER_INDICES(ordinal, component));
      }
    }
  }
}

FieldBase &
StkFieldSyncDebugger::get_last_mod_location_field() const
{
  if (m_lastModLocationField == nullptr) {
    ThrowRequire(impl::get_ngp_field(*m_stkField) != nullptr);
    BulkData & bulk = m_stkField->get_mesh();
    MetaData & meta = bulk.mesh_meta_data();
    meta.enable_late_fields();
    FieldState state = m_stkField->state();
    FieldBase* fieldWithStateNew = m_stkField->field_state(stk::mesh::StateNew);
    Field<uint8_t> & lastModLocationField =
        meta.declare_field<Field<uint8_t>>(m_stkField->entity_rank(),
                                           "DEBUG_lastFieldModLocation_"+fieldWithStateNew->name(),
                                           m_stkField->number_of_states());

    meta.set_mesh_on_fields(&bulk);
    const FieldBase::RestrictionVector & fieldRestrictions = m_stkField->restrictions();
    if (not fieldRestrictions.empty()) {
      for (const FieldBase::Restriction & restriction : fieldRestrictions) {
        const unsigned numComponents = restriction.num_scalars_per_entity();
        std::vector<uint8_t> initLastModLocation(numComponents, LastModLocation::HOST_OR_DEVICE);
        put_field_on_mesh(lastModLocationField, restriction.selector(), numComponents, initLastModLocation.data());
      }
    }
    else {
      bulk.reallocate_field_data(lastModLocationField);
    }

    m_lastModLocationField = lastModLocationField.field_state(state);
  }
  return *m_lastModLocationField;
}

void
StkFieldSyncDebugger::set_last_modification_view(const LastFieldModLocationType & lastModView) const
{
  m_debugFieldLastModification = lastModView;
}

void
StkFieldSyncDebugger::set_lost_device_field_data_view(const ScalarUvmType<bool>& lostDeviceFieldDataView) const
{
  m_lostDeviceFieldData = lostDeviceFieldDataView;
}

void
StkFieldSyncDebugger::set_bucket_offset_view(const UnsignedViewType::HostMirror & hostSelectedBucketOffset) const
{
  m_hostSelectedBucketOffset = hostSelectedBucketOffset;
}

void
StkFieldSyncDebugger::mark_data_initialized() const
{
  m_isDataInitialized = true;
}

void
StkFieldSyncDebugger::detect_host_field_entity_modification() const
{
  BulkData & bulk = m_stkField->get_mesh();
  if (!bulk.is_valid(m_lastFieldEntity)) return;

  for (unsigned component = 0; component < get_num_components(m_lastFieldEntity);  ++component) {
    if (last_accessed_entity_value_has_changed(m_lastFieldEntity, m_lastFieldValue.data(), component)) {
      set_last_modification(m_lastFieldEntity, component, LastModLocation::HOST);
    }
  }
}

void
StkFieldSyncDebugger::detect_host_field_bucket_modification() const
{
  BulkData & bulk = m_stkField->get_mesh();
  for (unsigned ordinal = 0; ordinal < m_lastFieldBucketEntities.size(); ++ordinal) {
    const stk::mesh::Entity & entity = m_lastFieldBucketEntities[ordinal];
    const MeshIndex & index = bulk.mesh_index(entity);
    const unsigned bytesPerEntity = m_stkField->get_meta_data_for_field()[index.bucket->bucket_id()].m_bytes_per_entity;
    const uint8_t * lastEntityValues = m_lastFieldBucketValues.data() + ordinal * bytesPerEntity;

    for (unsigned component = 0; component < get_num_components(entity);  ++component) {
      if (last_accessed_entity_value_has_changed(entity, lastEntityValues, component)) {
        set_last_modification(entity, component, LastModLocation::HOST);
      }
    }
  }
}

unsigned
StkFieldSyncDebugger::get_num_components(const Entity & entity) const
{
  BulkData & bulk = m_stkField->get_mesh();
  const MeshIndex & index = bulk.mesh_index(entity);
  const unsigned bytesPerScalar = m_stkField->data_traits().size_of;
  return m_stkField->get_meta_data_for_field()[index.bucket->bucket_id()].m_bytes_per_entity/bytesPerScalar;
}

bool
StkFieldSyncDebugger::last_accessed_entity_value_has_changed(const Entity & entity,
                                                             const uint8_t * lastEntityValues,
                                                             unsigned component) const
{
  BulkData & bulk = m_stkField->get_mesh();
  if (bulk.is_valid(entity)) {
    const MeshIndex & index = bulk.mesh_index(entity);
    const unsigned bytesPerScalar = m_stkField->data_traits().size_of;
    const FieldMetaData& field_meta_data = m_stkField->get_meta_data_for_field()[index.bucket->bucket_id()];
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
StkFieldSyncDebugger::set_last_modification(const Entity & entity, unsigned component, LastModLocation location) const
{
  BulkData & bulk = m_stkField->get_mesh();
  if (bulk.in_modifiable_state()) {
    stk::mesh::FieldBase & lastModLocationField = get_last_mod_location_field();
    uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(field_data<FieldBase, EmptyStkFieldSyncDebugger>(lastModLocationField,
                                                                                                            entity));
    lastModLocation[component] = location;
  }
  else {
    const MeshIndex & index = bulk.mesh_index(entity);
    m_debugFieldLastModification(m_hostSelectedBucketOffset(index.bucket->bucket_id()),
                                 ORDER_INDICES(index.bucket_ordinal, component)) = location;
  }
}

std::string
location_string(const char * fullPath, int lineNumber)
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
StkFieldSyncDebugger::check_stale_field_entity_access(const Entity& entity, const char* fileName, int lineNumber) const
{
  if (m_lostDeviceFieldData()) {
    std::cout << location_string(fileName, lineNumber)
              << "*** WARNING: Lost Device values for Field " << m_stkField->name()
              << " due to a mesh modification before a sync to Host" << std::endl;
    return;
  }

  BulkData & bulk = m_stkField->get_mesh();
  for (unsigned component = 0; component < get_num_components(entity);  ++component) {
    if (data_is_stale_on_host(entity, component)) {
      const MeshIndex & index = bulk.mesh_index(entity);
      const unsigned bytesPerScalar = m_stkField->data_traits().size_of;
      const FieldMetaData& field_meta_data = m_stkField->get_meta_data_for_field()[index.bucket->bucket_id()];
      const unsigned entityOffset = field_meta_data.m_bytes_per_entity * index.bucket_ordinal;
      const unsigned componentOffset = bytesPerScalar * component;

      uint8_t * data = field_meta_data.m_data + entityOffset + componentOffset;

      const std::string baseMessage = location_string(fileName, lineNumber) +
                                      "*** WARNING: Accessing stale data on Host for Field " + m_stkField->name() +
                                      "[" + std::to_string(component) + "]=";
      const DataTraits & dataTraits = m_stkField->data_traits();
      if (dataTraits.is_floating_point && dataTraits.size_of == 8u) {
        std::cout << baseMessage << *reinterpret_cast<double*>(data) << std::endl;
      }
      else if (dataTraits.is_floating_point && dataTraits.size_of == 4u) {
        std::cout << baseMessage << *reinterpret_cast<float*>(data) << std::endl;
      }
      else if (dataTraits.is_integral && dataTraits.is_signed && dataTraits.size_of == 4u) {
        std::cout << baseMessage << *reinterpret_cast<int32_t*>(data) << std::endl;
      }
      else if (dataTraits.is_integral && dataTraits.is_unsigned && dataTraits.size_of == 4u) {
        std::cout << baseMessage << *reinterpret_cast<uint32_t*>(data) << std::endl;
      }
      else if (dataTraits.is_integral && dataTraits.is_signed && dataTraits.size_of == 8u) {
        std::cout << baseMessage << *reinterpret_cast<int64_t*>(data) << std::endl;
      }
      else if (dataTraits.is_integral && dataTraits.is_unsigned && dataTraits.size_of == 8u) {
        std::cout << baseMessage << *reinterpret_cast<uint64_t*>(data) << std::endl;
      }
      else if (dataTraits.is_integral && dataTraits.is_signed && dataTraits.size_of == 2u) {
        std::cout << baseMessage << *reinterpret_cast<int16_t*>(data) << std::endl;
      }
      else if (dataTraits.is_integral && dataTraits.is_unsigned && dataTraits.size_of == 2u) {
        std::cout << baseMessage << *reinterpret_cast<uint16_t*>(data) << std::endl;
      }
      else {
        std::cout << location_string(fileName, lineNumber)
                  << "*** WARNING: Accessing stale data on Host for Field " << m_stkField->name() << std::endl;
      }
    }
  }
}

void
StkFieldSyncDebugger::check_stale_field_bucket_access(const stk::mesh::Bucket& bucket, const char* fileName, int lineNumber) const
{
  for (const Entity & entity : bucket) {
    check_stale_field_entity_access(entity, fileName, lineNumber);
  }
}

bool
StkFieldSyncDebugger::data_is_stale_on_host(const Entity& entity, unsigned component) const
{
  BulkData & bulk = m_stkField->get_mesh();
  if (bulk.in_modifiable_state()) {
    stk::mesh::FieldBase & lastModLocationField = get_last_mod_location_field();
    const uint8_t * lastModLocation = reinterpret_cast<uint8_t*>(field_data<FieldBase, EmptyStkFieldSyncDebugger>(lastModLocationField,
                                                                                                                  entity));
    return !(lastModLocation[component] & LastModLocation::HOST);
  }
  else {
    const MeshIndex & index = bulk.mesh_index(entity);
    return !(m_debugFieldLastModification(m_hostSelectedBucketOffset(index.bucket->bucket_id()),
                                          ORDER_INDICES(index.bucket_ordinal, component)) & LastModLocation::HOST);
  }
}

void
StkFieldSyncDebugger::store_last_entity_access_location(const Entity & entity) const
{
  BulkData & bulk = m_stkField->get_mesh();
  const MeshIndex & index = bulk.mesh_index(entity);
  const FieldMetaData& field_meta_data = m_stkField->get_meta_data_for_field()[index.bucket->bucket_id()];
  const uint8_t * data = field_meta_data.m_data + field_meta_data.m_bytes_per_entity * index.bucket_ordinal;

  m_lastFieldValue.resize(field_meta_data.m_bytes_per_entity);
  std::memcpy(m_lastFieldValue.data(), data, field_meta_data.m_bytes_per_entity);
  m_lastFieldEntity = entity;
}

void
StkFieldSyncDebugger::store_last_bucket_access_location(const Bucket & bucket) const
{
  const FieldMetaData& field_meta_data = m_stkField->get_meta_data_for_field()[bucket.bucket_id()];
  m_lastFieldBucketValues.resize(field_meta_data.m_bytes_per_entity * bucket.capacity());

  std::memcpy(m_lastFieldBucketValues.data(), field_meta_data.m_data, field_meta_data.m_bytes_per_entity * bucket.size());
  m_lastFieldBucketEntities.clear();
  for (const stk::mesh::Entity & entity : bucket) {
    m_lastFieldBucketEntities.push_back(entity);
  }
}


}
}
