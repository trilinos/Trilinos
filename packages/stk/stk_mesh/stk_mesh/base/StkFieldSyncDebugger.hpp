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

#ifndef STKFIELDSYNCDEBUGGER_HPP
#define STKFIELDSYNCDEBUGGER_HPP

#include "FieldSyncDebugging.hpp"
#include "Entity.hpp"
#include "Bucket.hpp"
#include "NgpTypes.hpp"
#include <cstddef>
#include <vector>

namespace stk {
namespace mesh {

class Bucket;
class FieldBase;

class EmptyStkFieldSyncDebugger
{
public:
  EmptyStkFieldSyncDebugger(const FieldBase*) {}
  ~EmptyStkFieldSyncDebugger() = default;

  inline void host_stale_access_entity_check(const stk::mesh::Entity &, const char *, int) {}
  inline void host_stale_access_entity_check(const unsigned &, const unsigned &, const char *, int) {}
  inline void host_stale_access_bucket_check(const stk::mesh::Bucket &, const char *, int) {}
  inline void host_stale_access_bucket_check(const unsigned &, const char *, int) {}
  inline void clear_last_field_value() {}
  inline void fill_last_mod_location_view_from_host() {}
  inline void fill_last_mod_location_field_from_device() {}
};

class StkFieldSyncDebugger
{
public:
  StkFieldSyncDebugger(const FieldBase* stkField);
  ~StkFieldSyncDebugger() = default;

  void host_stale_access_entity_check(const stk::mesh::Entity & entity, const char * fileName, int lineNumber);
  void host_stale_access_entity_check(const unsigned & bucketId, const unsigned & bucketOrd, const char * fileName, int lineNumber);
  void host_stale_access_bucket_check(const stk::mesh::Bucket & bucket, const char * fileName, int lineNumber);
  void host_stale_access_bucket_check(const unsigned & bucketId, const char * fileName, int lineNumber);
  void clear_last_field_value();
  void fill_last_mod_location_view_from_host();
  void fill_last_mod_location_field_from_device();

  void set_last_modification_view(const LastFieldModLocationType & lastModView) const;
  void set_lost_device_field_data_view(const ScalarUvmType<bool> & lostDeviceFieldDataView) const;
  void set_bucket_offset_view(const UnsignedViewType::HostMirror & hostSelectedBucketOffset) const;
  void mark_data_initialized() const;

private:
  FieldBase & get_last_mod_location_field() const;
  void detect_host_field_entity_modification() const;
  void detect_host_field_bucket_modification() const;
  unsigned get_num_components(const Entity& entity) const;
  bool last_accessed_entity_value_has_changed(const Entity& entity, const uint8_t* lastEntityValues, unsigned component) const;
  void set_last_modification(const Entity& entity, unsigned component, LastModLocation location) const;
  void check_stale_field_entity_access(const Entity& entity, const char* fileName, int lineNumber) const;
  void check_stale_field_bucket_access(const stk::mesh::Bucket& bucket, const char* fileName, int lineNumber) const;
  bool data_is_stale_on_host(const Entity& entity, unsigned component) const;
  void store_last_entity_access_location(const Entity& entity) const;
  void store_last_bucket_access_location(const Bucket& bucket) const;

  mutable const FieldBase*             m_stkField;
  mutable LastFieldModLocationType     m_debugFieldLastModification;
  mutable ScalarUvmType<bool>          m_lostDeviceFieldData;
  mutable UnsignedViewType::HostMirror m_hostSelectedBucketOffset;
  mutable Entity                       m_lastFieldEntity = stk::mesh::Entity();
  mutable std::vector<uint8_t>         m_lastFieldValue;
  mutable std::vector<Entity>          m_lastFieldBucketEntities;
  mutable std::vector<uint8_t>         m_lastFieldBucketValues;
  mutable FieldBase*                   m_lastModLocationField = nullptr;
  mutable bool                         m_isDataInitialized;
};

}
}

#endif // STKFIELDSYNCDEBUGGER_HPP
