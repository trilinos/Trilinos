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

#ifndef STK_CONSTFIELDDATA_HPP
#define STK_CONSTFIELDDATA_HPP

#include "FieldDataBytes.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_mesh/base/EntityValues.hpp"
#include "stk_mesh/base/BucketValues.hpp"
#include "stk_mesh/baseImpl/NgpFieldAux.hpp"
#include <any>

namespace stk::mesh {

//==============================================================================
// Device ConstFieldData
//==============================================================================

template <typename T,
          typename Space = stk::ngp::HostSpace,
          Layout DataLayout = DefaultLayoutSelector<Space>::layout>
class ConstFieldData : public FieldDataBytes<Space>
{
public:
  using value_type = T;
  using space = Space;
  using exec_space = typename Space::exec_space;
  using mem_space = typename Space::mem_space;
  static constexpr Layout layout = DataLayout;

  KOKKOS_FUNCTION ConstFieldData();
  ConstFieldData(FieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes);
  KOKKOS_FUNCTION virtual ~ConstFieldData() override;

  ConstFieldData(const ConstFieldData& fieldData, FieldAccessTag accessTag);
  KOKKOS_FUNCTION ConstFieldData(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData(ConstFieldData&& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(ConstFieldData&& fieldData);

  KOKKOS_INLINE_FUNCTION
  EntityValues<const T, Space, DataLayout> entity_values(Entity entity,
                                                         const char* file = STK_DEVICE_FILE,
                                                         int line = STK_DEVICE_LINE) const;

  KOKKOS_INLINE_FUNCTION
  EntityValues<const T, Space, DataLayout> entity_values(const FastMeshIndex& fmi,
                                                         const char* file = STK_DEVICE_FILE,
                                                         int line = STK_DEVICE_LINE) const;

  KOKKOS_INLINE_FUNCTION
  BucketValues<const T, Space, DataLayout> bucket_values(int bucketId,
                                                         const char* file = STK_DEVICE_FILE,
                                                         int line = STK_DEVICE_LINE) const;

protected:
  template <typename _T, typename _MemSpace> friend class DeviceField;

  virtual void update_device_field_metadata() override;
  virtual void sync_to_host(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout) override;
  virtual void sync_to_device(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout) override;
  virtual void update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout, bool needsSync) override;
  virtual void fence(const stk::ngp::ExecSpace& execSpace) override;
};


//==============================================================================
// Host ConstFieldData
//==============================================================================

template <typename T>
class ConstFieldData<T, stk::ngp::HostSpace, Layout::Right> : public FieldDataBytes<stk::ngp::HostSpace>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Right;

  ConstFieldData();
  ConstFieldData(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName,
                 const DataTraits& dataTraits);
  KOKKOS_FUNCTION ~ConstFieldData() override;

  ConstFieldData(const ConstFieldData& fieldData, FieldAccessTag accessTag);
  KOKKOS_FUNCTION ConstFieldData(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData(ConstFieldData&& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(ConstFieldData&& fieldData);

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Right> entity_values(Entity entity,
                                                                          const char* file = STK_HOST_FILE,
                                                                          int line = STK_HOST_LINE) const;

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Right> entity_values(const MeshIndex& mi,
                                                                          const char* file = STK_HOST_FILE,
                                                                          int line = STK_HOST_LINE) const;

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Right> entity_values(const FastMeshIndex& fmi,
                                                                          const char* file = STK_HOST_FILE,
                                                                          int line = STK_HOST_LINE) const;

  inline
  BucketValues<const T, stk::ngp::HostSpace, Layout::Right> bucket_values(const Bucket& bucket,
                                                                          const char* file = STK_HOST_FILE,
                                                                          int line = STK_HOST_LINE) const;

  inline
  BucketValues<const T, stk::ngp::HostSpace, Layout::Right> bucket_values(int bucketId,
                                                                          const char* file = STK_HOST_FILE,
                                                                          int line = STK_HOST_LINE) const;

protected:
  template <typename _T, typename _MemSpace> friend class HostField;

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout, bool needsSync) override;
  virtual void fence(const stk::ngp::ExecSpace&) override {}
};

//------------------------------------------------------------------------------
template <typename T>
class ConstFieldData<T, stk::ngp::HostSpace, Layout::Left> : public FieldDataBytes<stk::ngp::HostSpace>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Left;

  ConstFieldData();
  ConstFieldData(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName,
                 const DataTraits& dataTraits);
  KOKKOS_FUNCTION ~ConstFieldData() override;

  ConstFieldData(const ConstFieldData& fieldData, FieldAccessTag accessTag);
  KOKKOS_FUNCTION ConstFieldData(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData(ConstFieldData&& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(ConstFieldData&& fieldData);

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Left> entity_values(Entity entity,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Left> entity_values(const MeshIndex& mi,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Left> entity_values(const FastMeshIndex& fmi,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  BucketValues<const T, stk::ngp::HostSpace, Layout::Left> bucket_values(const Bucket& bucket,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  BucketValues<const T, stk::ngp::HostSpace, Layout::Left> bucket_values(int bucketId,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

protected:
  template <typename _T, typename _MemSpace> friend class HostField;

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout, bool needsSync) override;
  virtual void fence(const stk::ngp::ExecSpace&) override {}
};

//------------------------------------------------------------------------------
template <typename T>
class ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto> : public FieldDataBytes<stk::ngp::HostSpace>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Auto;

  ConstFieldData();
  ConstFieldData(const FieldDataBytes<stk::ngp::HostSpace>& hostFieldBytes, FieldAccessTag accessTag);
  KOKKOS_FUNCTION ~ConstFieldData() override;

  KOKKOS_FUNCTION ConstFieldData(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData(ConstFieldData&& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(const ConstFieldData& fieldData);
  KOKKOS_FUNCTION ConstFieldData& operator=(ConstFieldData&& fieldData);

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Auto> entity_values(Entity entity,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Auto> entity_values(const MeshIndex& mi,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  EntityValues<const T, stk::ngp::HostSpace, Layout::Auto> entity_values(const FastMeshIndex& fmi,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  BucketValues<const T, stk::ngp::HostSpace, Layout::Auto> bucket_values(const Bucket& bucket,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

  inline
  BucketValues<const T, stk::ngp::HostSpace, Layout::Auto> bucket_values(int bucketId,
                                                                         const char* file = STK_HOST_FILE,
                                                                         int line = STK_HOST_LINE) const;

protected:
  template <typename _T, typename _MemSpace> friend class HostField;

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout, bool needsSync) override;
  virtual void fence(const stk::ngp::ExecSpace&) override {}
};


//==============================================================================
// Device ConstFieldData definitions
//==============================================================================

template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION
ConstFieldData<T, Space, DataLayout>::ConstFieldData()
  : FieldDataBytes<Space>()
{}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
ConstFieldData<T, Space, DataLayout>::ConstFieldData(FieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes)
  : FieldDataBytes<Space>(hostFieldBytes, DataLayout)
{}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION
ConstFieldData<T, Space, DataLayout>::~ConstFieldData()
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();
  )
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
ConstFieldData<T, Space, DataLayout>::ConstFieldData(const ConstFieldData& fieldData, FieldAccessTag accessTag)
  : FieldDataBytes<Space>(fieldData)
{
  this->track_copy(accessTag);
  this->update_field_meta_data_mod_count();
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION
ConstFieldData<T, Space, DataLayout>::ConstFieldData(const ConstFieldData& fieldData)
  : FieldDataBytes<Space>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION
ConstFieldData<T, Space, DataLayout>::ConstFieldData(ConstFieldData&& fieldData)
  : FieldDataBytes<Space>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION ConstFieldData<T, Space, DataLayout>&
ConstFieldData<T, Space, DataLayout>::operator=(const ConstFieldData& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<Space>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  KOKKOS_IF_ON_DEVICE(
    FieldDataBytes<Space>::operator=(fieldData);
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION ConstFieldData<T, Space, DataLayout>&
ConstFieldData<T, Space, DataLayout>::operator=(ConstFieldData&& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<Space>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  KOKKOS_IF_ON_DEVICE(
    FieldDataBytes<Space>::operator=(fieldData);
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_INLINE_FUNCTION EntityValues<const T, Space, DataLayout>
ConstFieldData<T, Space, DataLayout>::entity_values(Entity entity,
                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_entity_local_offset(entity.local_offset(), file, line);

  const FastMeshIndex& fmi = this->m_deviceFastMeshIndices[entity.local_offset()];

  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];

  return EntityValues<const T, Space, DataLayout>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}


//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_INLINE_FUNCTION EntityValues<const T, Space, DataLayout>
ConstFieldData<T, Space, DataLayout>::entity_values(const FastMeshIndex& fmi,
                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];

  return EntityValues<const T, Space, DataLayout>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_INLINE_FUNCTION BucketValues<const T, Space, DataLayout>
ConstFieldData<T, Space, DataLayout>::bucket_values(int bucketId,
                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[bucketId];

  return BucketValues<const T, Space, DataLayout>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
void
ConstFieldData<T, Space, DataLayout>::update_device_field_metadata()
{
  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<Space>(this->mesh());
  deviceFieldDataManager->set_device_field_meta_data(*this);
  this->m_deviceMeshSyncCountAtLastSync = deviceFieldDataManager->synchronized_count();
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
void
ConstFieldData<T, Space, DataLayout>::sync_to_host(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout)
{
  impl::transpose_to_pinned_and_mapped_memory<T>(execSpace, this->m_deviceFieldMetaData, hostDataLayout);
  execSpace.fence();
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
void
ConstFieldData<T, Space, DataLayout>::sync_to_device(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout)
{
  impl::transpose_from_pinned_and_mapped_memory<T>(execSpace, this->m_deviceFieldMetaData, hostDataLayout);
  execSpace.fence();
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
void
ConstFieldData<T, Space, DataLayout>::update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout,
                                             bool needsSync)
{
  ProfilingBlock prof("ConstFieldData::update()");

  this->m_fieldDataSynchronizedCount = this->mesh().synchronized_count();

  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<Space>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  if (not deviceFieldDataManager->update_all_bucket_allocations()) {
    deviceFieldDataManager->update_host_bucket_pointers(this->field_ordinal());
  }

  deviceFieldDataManager->set_device_field_meta_data(*this);


  if (not deviceFieldDataManager->has_unified_device_storage(this->field_ordinal())) {
    if (needsSync) {
      // Sync everything in one shot, in addition to just the modified Buckets
      impl::transpose_from_pinned_and_mapped_memory<T>(execSpace, this->m_deviceFieldMetaData, hostDataLayout);
    }
    else {
      // Just sync the modified Buckets
      int rankedOrdinal = -1;
      const auto deviceBucketsModified = std::any_cast<DeviceBucketsModifiedCollectionType<mem_space>>(
              deviceFieldDataManager->get_device_bucket_is_modified(this->field_ordinal(), rankedOrdinal));

      impl::transpose_modified_buckets_to_device<T>(execSpace, this->m_deviceFieldMetaData, rankedOrdinal,
                                                    deviceBucketsModified, hostDataLayout);
    }
  }
  execSpace.fence();

  deviceFieldDataManager->clear_bucket_is_modified(this->field_ordinal());

  this->m_deviceFastMeshIndices = this->mesh().template get_updated_fast_mesh_indices<mem_space>();
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
void
ConstFieldData<T, Space, DataLayout>::fence(const stk::ngp::ExecSpace& execSpace)
{
  execSpace.fence();
}


//==============================================================================
// Host ConstFieldData definitions: Layout::Right
//==============================================================================

template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::ConstFieldData()
  : FieldDataBytes<stk::ngp::HostSpace>()
{}

//------------------------------------------------------------------------------
template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::ConstFieldData(EntityRank entityRank, Ordinal fieldOrdinal,
                                                                      const std::string& fieldName,
                                                                      const DataTraits& dataTraits)
  : FieldDataBytes<stk::ngp::HostSpace>(entityRank, fieldOrdinal, fieldName, dataTraits, Layout::Right)
{}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::~ConstFieldData()
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();
  )
}

//------------------------------------------------------------------------------
template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::ConstFieldData(const ConstFieldData& fieldData,
                                                                      FieldAccessTag accessTag)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  this->track_copy(accessTag);
  this->update_field_meta_data_mod_count();
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::ConstFieldData(const ConstFieldData& fieldData)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::ConstFieldData(ConstFieldData&& fieldData)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>&
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::operator=(const ConstFieldData& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<stk::ngp::HostSpace>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>&
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::operator=(ConstFieldData&& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<stk::ngp::HostSpace>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Right>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::entity_values(Entity entity,
                                                                     const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_updated_field(file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<const T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Right>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::entity_values(const MeshIndex& mi,
                                                                     const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);
  this->check_bucket_ordinal(mi.bucket->bucket_id(), mi.bucket_ordinal, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<const T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Right>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::entity_values(const FastMeshIndex& fmi,
                                                                     const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityValues<const T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<const T, stk::ngp::HostSpace, Layout::Right>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::bucket_values(const Bucket& bucket,
                                                                     const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketValues<const T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<const T, stk::ngp::HostSpace, Layout::Right>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::bucket_values(int bucketId,
                                                                     const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketValues<const T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
void
ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>::update(const stk::ngp::ExecSpace&, Layout, bool)
{
  this->m_fieldDataSynchronizedCount = this->mesh().synchronized_count();
}


//==============================================================================
// Host ConstFieldData definitions: Layout::Left
//==============================================================================

template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::ConstFieldData()
  : FieldDataBytes<stk::ngp::HostSpace>()
{}

//------------------------------------------------------------------------------
template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::ConstFieldData(EntityRank entityRank, Ordinal fieldOrdinal,
                                                                     const std::string& fieldName,
                                                                     const DataTraits& dataTraits)
  : FieldDataBytes<stk::ngp::HostSpace>(entityRank, fieldOrdinal, fieldName, dataTraits, Layout::Left)
{}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::~ConstFieldData()
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();
  )
}

//------------------------------------------------------------------------------
template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::ConstFieldData(const ConstFieldData& fieldData,
                                                                     FieldAccessTag accessTag)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  this->track_copy(accessTag);
  this->update_field_meta_data_mod_count();
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::ConstFieldData(const ConstFieldData& fieldData)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::ConstFieldData(ConstFieldData&& fieldData)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>&
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::operator=(const ConstFieldData& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<stk::ngp::HostSpace>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>&
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::operator=(ConstFieldData&& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<stk::ngp::HostSpace>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Left>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::entity_values(Entity entity,
                                                                    const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_updated_field(file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<const T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Left>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::entity_values(const MeshIndex& mi,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);
  this->check_bucket_ordinal(mi.bucket->bucket_id(), mi.bucket_ordinal, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<const T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Left>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::entity_values(const FastMeshIndex& fmi,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityValues<const T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<const T, stk::ngp::HostSpace, Layout::Left>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::bucket_values(const Bucket& bucket,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketValues<const T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<const T, stk::ngp::HostSpace, Layout::Left>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::bucket_values(int bucketId,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketValues<const T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
void
ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>::update(const stk::ngp::ExecSpace&, Layout, bool)
{
  this->m_fieldDataSynchronizedCount = this->mesh().synchronized_count();
}


//==============================================================================
// Host ConstFieldData definitions: Layout::Auto
//==============================================================================

template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::ConstFieldData()
  : FieldDataBytes<stk::ngp::HostSpace>()
{}

//------------------------------------------------------------------------------
template <typename T>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::ConstFieldData(
    const FieldDataBytes<stk::ngp::HostSpace>& hostFieldBytes, FieldAccessTag accessTag)
  : FieldDataBytes<stk::ngp::HostSpace>(hostFieldBytes)
{
  this->track_copy(accessTag);
  this->update_field_meta_data_mod_count();
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::~ConstFieldData()
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::ConstFieldData(const ConstFieldData& fieldData)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::ConstFieldData(ConstFieldData&& fieldData)
  : FieldDataBytes<stk::ngp::HostSpace>(fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->track_copy(this->access_tag());
  )
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>&
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::operator=(const ConstFieldData& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<stk::ngp::HostSpace>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T>
KOKKOS_FUNCTION
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>&
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::operator=(ConstFieldData&& fieldData)
{
  KOKKOS_IF_ON_HOST(
    this->release_copy();  // Decrement first if becoming untracked
    FieldDataBytes<stk::ngp::HostSpace>::operator=(fieldData);
    this->track_copy(fieldData.access_tag());  // Increment after if becoming tracked
  )

  return *this;
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::entity_values(Entity entity,
                                                                    const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_updated_field(file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity, this->field_name());
  }
  else if (m_layout == Layout::Left) {
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout << ".  The actual run-time layout must be "
                                                                      "either Layout::Right or Layout::Left.");
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::entity_values(const MeshIndex& mi,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);
  this->check_bucket_ordinal(mi.bucket->bucket_id(), mi.bucket_ordinal, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity, this->field_name());
  }
  else if (m_layout == Layout::Left) {
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout << ".  The actual run-time layout must be "
                                                                      "either Layout::Right or Layout::Left.");
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::entity_values(const FastMeshIndex& fmi,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  if (m_layout == Layout::Right) {
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity, this->field_name());
  }
  else if (m_layout == Layout::Left) {
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout << ".  The actual run-time layout must be "
                                                                      "either Layout::Right or Layout::Left.");
    return EntityValues<const T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::bucket_values(const Bucket& bucket,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  if (m_layout == Layout::Right) {
    return BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize, this->field_name());
  }
  else if (m_layout == Layout::Left) {
    return BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout << ".  The actual run-time layout must be "
                                                                      "either Layout::Right or Layout::Left.");
    return BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::bucket_values(int bucketId,
                                                                    const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  if (m_layout == Layout::Right) {
    return BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize, this->field_name());
  }
  else if (m_layout == Layout::Left) {
    return BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout << ".  The actual run-time layout must be "
                                                                      "either Layout::Right or Layout::Left.");
    return BucketValues<const T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
void
ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>::update(const stk::ngp::ExecSpace&, Layout, bool)
{
  this->m_fieldDataSynchronizedCount = this->mesh().synchronized_count();
}


//==============================================================================
}

#endif // STK_CONSTFIELDDATA_HPP
