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
FieldBase::get_num_components(const FastMeshIndex & index) const
{
  const unsigned bytesPerScalar = data_traits().size_of;
  return m_field_meta_data[index.bucket_id].m_bytes_per_entity/bytesPerScalar;
}

void
FieldBase::set_last_modification_view(const LastFieldModLocationType & lastModView) const
{
  m_debugFieldLastModification = lastModView;
}

bool
FieldBase::last_accessed_entity_value_has_changed(const FastMeshIndex & index, unsigned component) const
{
  const unsigned bytesPerScalar = data_traits().size_of;
  const FieldMetaData& field_meta_data = m_field_meta_data[index.bucket_id];
  const unsigned entityOffset = field_meta_data.m_bytes_per_entity * index.bucket_ord;
  const unsigned componentOffset = bytesPerScalar * component;

  const uint8_t * data = field_meta_data.m_data + entityOffset + componentOffset;
  const uint8_t * lastValueData = m_lastFieldValue.data() + entityOffset + componentOffset;

  return std::memcmp(data, lastValueData, bytesPerScalar);
}

void
FieldBase::set_last_modification(const FastMeshIndex & index, unsigned component, LastModLocation location) const
{
  m_debugFieldLastModification(index.bucket_id, ORDER_INDICES(index.bucket_ord, component)) = location;
}

void
FieldBase::detect_host_field_entity_modification() const
{
  if (m_lastFieldEntityLocation.bucket_id == INVALID_BUCKET_ID) return;

  for (unsigned component = 0; component < get_num_components(m_lastFieldEntityLocation);  ++component) {
    if (last_accessed_entity_value_has_changed(m_lastFieldEntityLocation, component)) {
      set_last_modification(m_lastFieldEntityLocation, component, LastModLocation::HOST);
    }
  }
}

bool
FieldBase::last_entity_modification_not_on_host(const FastMeshIndex & index, unsigned component) const
{
  return !(m_debugFieldLastModification(index.bucket_id,
                                        ORDER_INDICES(index.bucket_ord, component)) & LastModLocation::HOST);
}

void
FieldBase::store_last_entity_access_location(const FastMeshIndex & index) const
{
  const FieldMetaData& field_meta_data = m_field_meta_data[index.bucket_id];
  const uint8_t * data = field_meta_data.m_data + field_meta_data.m_bytes_per_entity * index.bucket_ord;

  const unsigned bucketCapacity = get_mesh().buckets(entity_rank())[index.bucket_id]->capacity();
  m_lastFieldValue.resize(field_meta_data.m_bytes_per_entity * bucketCapacity);
  uint8_t * lastValueData = m_lastFieldValue.data() + field_meta_data.m_bytes_per_entity * index.bucket_ord;

  std::memcpy(lastValueData, data, field_meta_data.m_bytes_per_entity);
  m_lastFieldEntityLocation = index;
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
FieldBase::check_stale_field_entity_access(const FastMeshIndex & index, const char * fileName, int lineNumber) const
{
  if (get_ngp_field()->any_device_field_modification()) {
    std::cout << location_string(fileName, lineNumber)
              << "*** WARNING: Lost Device values for Field " << name()
              << " due to a mesh modification before a sync to Host" << std::endl;
    return;
  }

  for (unsigned component = 0; component < get_num_components(index);  ++component) {
    if (last_entity_modification_not_on_host(index, component)) {
      const unsigned bytesPerScalar = data_traits().size_of;
      const FieldMetaData& field_meta_data = m_field_meta_data[index.bucket_id];
      const unsigned entityOffset = field_meta_data.m_bytes_per_entity * index.bucket_ord;
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
  m_lastFieldValue.resize(field_meta_data.m_bytes_per_entity * bucket.capacity());

  std::memcpy(m_lastFieldValue.data(), field_meta_data.m_data, field_meta_data.m_bytes_per_entity * bucket.size());
  m_lastFieldBucketLocation = &const_cast<stk::mesh::Bucket&>(bucket);
}

void
FieldBase::detect_host_field_bucket_modification() const
{
  if (m_lastFieldBucketLocation == nullptr) return;

  for (unsigned ordinal = 0; ordinal < m_lastFieldBucketLocation->size(); ++ordinal) {
    const FastMeshIndex index {m_lastFieldBucketLocation->bucket_id(), ordinal};
    for (unsigned component = 0; component < get_num_components(index);  ++component) {
      if (last_accessed_entity_value_has_changed(index, component)) {
        set_last_modification(index, component, LastModLocation::HOST);
      }
    }
  }
}

void
FieldBase::check_stale_field_bucket_access(const stk::mesh::Bucket & bucket, const char * fileName, int lineNumber) const
{
  for (unsigned ordinal = 0; ordinal < bucket.size(); ++ordinal) {
    const FastMeshIndex index {bucket.bucket_id(), ordinal};
    check_stale_field_entity_access(index, fileName, lineNumber);
  }
}

void
FieldBase::reset_debug_state()
{
  m_lastFieldEntityLocation = {INVALID_BUCKET_ID, 0};
  m_lastFieldBucketLocation = nullptr;
}

#endif

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

