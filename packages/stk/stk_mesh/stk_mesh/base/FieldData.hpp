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

#ifndef STK_FIELDDATA_HPP
#define STK_FIELDDATA_HPP

#include "ConstFieldData.hpp"

namespace stk::mesh {

//==============================================================================
// Device FieldData
//==============================================================================

template <typename T,
          typename Space = stk::ngp::HostSpace,
          Layout DataLayout = DefaultLayoutSelector<Space>::layout>
class FieldData : public ConstFieldData<T, Space, DataLayout>
{
public:
  using value_type = T;
  using space = Space;
  using exec_space = typename Space::exec_space;
  using mem_space = typename Space::mem_space;
  static constexpr Layout layout = DataLayout;

  KOKKOS_FUNCTION FieldData();
  FieldData(FieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes);
  KOKKOS_FUNCTION virtual ~FieldData() override {}

  FieldData(const FieldData& fieldData, FieldAccessTag accessTag);
  KOKKOS_DEFAULTED_FUNCTION FieldData(const FieldData& fieldData) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData(FieldData&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(const FieldData&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(FieldData&&) = default;

  KOKKOS_INLINE_FUNCTION
  EntityValues<T, Space, DataLayout> entity_values(Entity entity,
                                                   const char* file = STK_DEVICE_FILE,
                                                   int line = STK_DEVICE_LINE) const;

  KOKKOS_INLINE_FUNCTION
  EntityValues<T, Space, DataLayout> entity_values(const FastMeshIndex& fmi,
                                                   const char* file = STK_DEVICE_FILE,
                                                   int line = STK_DEVICE_LINE) const;

  KOKKOS_INLINE_FUNCTION
  BucketValues<T, Space, DataLayout> bucket_values(int bucketId,
                                                   const char* file = STK_DEVICE_FILE,
                                                   int line = STK_DEVICE_LINE) const;
};


//==============================================================================
// Host FieldData
//==============================================================================

template <typename T>
class FieldData<T, stk::ngp::HostSpace, Layout::Right>
    : public ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Right;

  FieldData();
  FieldData(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName, const DataTraits& dataTraits);
  KOKKOS_DEFAULTED_FUNCTION virtual ~FieldData() override = default;

  FieldData(const FieldData& fieldData, FieldAccessTag accessTag);
  KOKKOS_DEFAULTED_FUNCTION FieldData(const FieldData& fieldData) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData(FieldData&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(const FieldData&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(FieldData&&) = default;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Right> entity_values(Entity entity,
                                                                    const char* file = STK_HOST_FILE,
                                                                    int line = STK_HOST_LINE) const;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Right> entity_values(const MeshIndex& mi,
                                                                    const char* file = STK_HOST_FILE,
                                                                    int line = STK_HOST_LINE) const;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Right> entity_values(const FastMeshIndex& fmi,
                                                                    const char* file = STK_HOST_FILE,
                                                                    int line = STK_HOST_LINE) const;

  inline
  BucketValues<T, stk::ngp::HostSpace, Layout::Right> bucket_values(const Bucket& bucket,
                                                                    const char* file = STK_HOST_FILE,
                                                                    int line = STK_HOST_LINE) const;

  inline
  BucketValues<T, stk::ngp::HostSpace, Layout::Right> bucket_values(int bucketId,
                                                                    const char* file = STK_HOST_FILE,
                                                                    int line = STK_HOST_LINE) const;
};

//------------------------------------------------------------------------------
template <typename T>
class FieldData<T, stk::ngp::HostSpace, Layout::Left>
    : public ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Left;

  FieldData();
  FieldData(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName, const DataTraits& dataTraits);
  KOKKOS_DEFAULTED_FUNCTION virtual ~FieldData() override = default;

  FieldData(const FieldData& fieldData, FieldAccessTag accessTag);
  KOKKOS_DEFAULTED_FUNCTION FieldData(const FieldData& fieldData) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData(FieldData&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(const FieldData&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(FieldData&&) = default;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Left> entity_values(Entity entity,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Left> entity_values(const MeshIndex& mi,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Left> entity_values(const FastMeshIndex& fmi,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  BucketValues<T, stk::ngp::HostSpace, Layout::Left> bucket_values(const Bucket& bucket,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  BucketValues<T, stk::ngp::HostSpace, Layout::Left> bucket_values(int bucketId,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;
};

//------------------------------------------------------------------------------
template <typename T>
class FieldData<T, stk::ngp::HostSpace, Layout::Auto>
    : public ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Auto;

  FieldData();
  FieldData(const FieldDataBytes<stk::ngp::HostSpace>& hostFieldBytes, FieldAccessTag accessTag);
  KOKKOS_DEFAULTED_FUNCTION virtual ~FieldData() override = default;

  KOKKOS_DEFAULTED_FUNCTION FieldData(const FieldData& fieldData) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData(FieldData&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(const FieldData&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldData& operator=(FieldData&&) = default;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Auto> entity_values(Entity entity,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Auto> entity_values(const MeshIndex& mi,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  EntityValues<T, stk::ngp::HostSpace, Layout::Auto> entity_values(const FastMeshIndex& fmi,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  BucketValues<T, stk::ngp::HostSpace, Layout::Auto> bucket_values(const Bucket& bucket,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;

  inline
  BucketValues<T, stk::ngp::HostSpace, Layout::Auto> bucket_values(int bucketId,
                                                                   const char* file = STK_HOST_FILE,
                                                                   int line = STK_HOST_LINE) const;
};


//==============================================================================
// Device FieldData definitions
//==============================================================================

template <typename T, typename Space, Layout DataLayout>
KOKKOS_FUNCTION
FieldData<T, Space, DataLayout>::FieldData()
  : ConstFieldData<T, Space, DataLayout>()
{}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
FieldData<T, Space, DataLayout>::FieldData(FieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes)
  : ConstFieldData<T, Space, DataLayout>(hostFieldBytes)
{
  static_assert(DataLayout == Layout::Left, "Only Layout::Left is supported for device data");
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
FieldData<T, Space, DataLayout>::FieldData(const FieldData& fieldData, FieldAccessTag accessTag)
  : ConstFieldData<T, Space, DataLayout>(fieldData, accessTag)
{}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_INLINE_FUNCTION EntityValues<T, Space, DataLayout>
FieldData<T, Space, DataLayout>::entity_values(Entity entity,
                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_entity_local_offset(entity.local_offset(), file, line);

  const FastMeshIndex& fmi = this->m_deviceFastMeshIndices[entity.local_offset()];

  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];

  return EntityValues<T, Space, DataLayout>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_INLINE_FUNCTION EntityValues<T, Space, DataLayout>
FieldData<T, Space, DataLayout>::entity_values(const FastMeshIndex& fmi,
                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];

  return EntityValues<T, Space, DataLayout>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T, typename Space, Layout DataLayout>
KOKKOS_INLINE_FUNCTION BucketValues<T, Space, DataLayout>
FieldData<T, Space, DataLayout>::bucket_values(int bucketId,
                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[bucketId];

  return BucketValues<T, Space, DataLayout>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity, this->field_name());
}


//==============================================================================
// Host FieldData definitions: Layout::Right
//==============================================================================

template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::FieldData()
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>()
{}

//------------------------------------------------------------------------------
template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::FieldData(EntityRank entityRank, Ordinal fieldOrdinal,
                                                            const std::string& fieldName,
                                                            const DataTraits& dataTraits)
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>(entityRank, fieldOrdinal, fieldName, dataTraits)
{}

//------------------------------------------------------------------------------
template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::FieldData(const FieldData& fieldData, FieldAccessTag accessTag)
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Right>(fieldData, accessTag)
{}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Right>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::entity_values(Entity entity,
                                                                const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_updated_field(file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Right>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::entity_values(const MeshIndex& mi,
                                                                const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);
  this->check_bucket_ordinal(mi.bucket->bucket_id(), mi.bucket_ordinal, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Right>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::entity_values(const FastMeshIndex& fmi,
                                                                const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityValues<T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<T, stk::ngp::HostSpace, Layout::Right>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::bucket_values(const Bucket& bucket,
                                                                const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketValues<T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<T, stk::ngp::HostSpace, Layout::Right>
FieldData<T, stk::ngp::HostSpace, Layout::Right>::bucket_values(int bucketId,
                                                                const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketValues<T, stk::ngp::HostSpace, Layout::Right>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize, this->field_name());
}


//==============================================================================
// Host FieldData definitions: Layout::Left
//==============================================================================

template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::FieldData()
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>()
{}

//------------------------------------------------------------------------------
template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::FieldData(EntityRank entityRank, Ordinal fieldOrdinal,
                                                           const std::string& fieldName,
                                                           const DataTraits& dataTraits)
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>(entityRank, fieldOrdinal, fieldName, dataTraits)
{}

//------------------------------------------------------------------------------
template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::FieldData(const FieldData& fieldData, FieldAccessTag accessTag)
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Left>(fieldData, accessTag)
{}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Left>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::entity_values(Entity entity,
                                                               const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_updated_field(file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Left>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::entity_values(const MeshIndex& mi,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);
  this->check_bucket_ordinal(mi.bucket->bucket_id(), mi.bucket_ordinal, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityValues<T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Left>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::entity_values(const FastMeshIndex& fmi,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityValues<T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<T, stk::ngp::HostSpace, Layout::Left>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::bucket_values(const Bucket& bucket,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketValues<T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity, this->field_name());
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<T, stk::ngp::HostSpace, Layout::Left>
FieldData<T, stk::ngp::HostSpace, Layout::Left>::bucket_values(int bucketId,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketValues<T, stk::ngp::HostSpace, Layout::Left>(
        reinterpret_cast<T*>(fieldMetaData.m_data),
        fieldMetaData.m_numComponentsPerEntity,
        fieldMetaData.m_numCopiesPerEntity,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity, this->field_name());
}


//==============================================================================
// Host FieldData definitions: Layout::Auto
//==============================================================================

template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::FieldData()
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>()
{}

//------------------------------------------------------------------------------
template <typename T>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::FieldData(
    const FieldDataBytes<stk::ngp::HostSpace>& hostFieldBytes, FieldAccessTag accessTag)
  : ConstFieldData<T, stk::ngp::HostSpace, Layout::Auto>(hostFieldBytes, accessTag)
{}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Auto>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::entity_values(Entity entity,
                                                               const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_updated_field(file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (this->m_layout == Layout::Right) {
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity, this->field_name());
  }
  else if (this->m_layout == Layout::Left) {
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << this->m_layout << ".  The actual run-time layout must be "
                                                                            "either Layout::Right or Layout::Left.");
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Auto>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::entity_values(const MeshIndex& mi,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);
  this->check_bucket_ordinal(mi.bucket->bucket_id(), mi.bucket_ordinal, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (this->m_layout == Layout::Right) {
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity, this->field_name());
  }
  else if (this->m_layout == Layout::Left) {
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data) + mi.bucket_ordinal,
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << this->m_layout << ".  The actual run-time layout must be "
                                                                            "either Layout::Right or Layout::Left.");
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline EntityValues<T, stk::ngp::HostSpace, Layout::Auto>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::entity_values(const FastMeshIndex& fmi,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);
  this->check_bucket_ordinal(fmi.bucket_id, fmi.bucket_ord, file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  if (this->m_layout == Layout::Right) {
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity, this->field_name());
  }
  else if (this->m_layout == Layout::Left) {
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data) + fmi.bucket_ord,
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << this->m_layout << ".  The actual run-time layout must be "
                                                                            "either Layout::Right or Layout::Left.");
    return EntityValues<T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<T, stk::ngp::HostSpace, Layout::Auto>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::bucket_values(const Bucket& bucket,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  if (this->m_layout == Layout::Right) {
    return BucketValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize, this->field_name());
  }
  else if (this->m_layout == Layout::Left) {
    return BucketValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << this->m_layout << ".  The actual run-time layout must be "
                                                                            "either Layout::Right or Layout::Left.");
    return BucketValues<T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, nullptr);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <typename T>
inline BucketValues<T, stk::ngp::HostSpace, Layout::Auto>
FieldData<T, stk::ngp::HostSpace, Layout::Auto>::bucket_values(int bucketId,
                                                               const char* file, int line) const
{
  this->check_updated_field(file, line);
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  if (this->m_layout == Layout::Right) {
    return BucketValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize, this->field_name());
  }
  else if (this->m_layout == Layout::Left) {
    return BucketValues<T, stk::ngp::HostSpace, Layout::Auto>(
          reinterpret_cast<T*>(fieldMetaData.m_data),
          fieldMetaData.m_numComponentsPerEntity,
          fieldMetaData.m_numCopiesPerEntity,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity, this->field_name());
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << this->m_layout << ".  The actual run-time layout must be "
                                                                            "either Layout::Right or Layout::Left.");
    return BucketValues<T, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, nullptr);  // Keep compiler happy
  }
}

//==============================================================================


}

#endif // STK_FIELDDATA_HPP
