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

#ifndef STK_FIELDDATABYTES_HPP
#define STK_FIELDDATABYTES_HPP

#include "ConstFieldDataBytes.hpp"

namespace stk::mesh {

//==============================================================================
// Device FieldDataBytes
//==============================================================================

template <typename Space = stk::ngp::HostSpace>
class FieldDataBytes : public ConstFieldDataBytes<Space>
{
public:
  KOKKOS_FUNCTION FieldDataBytes();
  FieldDataBytes(FieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes, Layout dataLayout);
  KOKKOS_FUNCTION virtual ~FieldDataBytes() override {}

  KOKKOS_DEFAULTED_FUNCTION FieldDataBytes(const FieldDataBytes& fieldData) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBytes(FieldDataBytes&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBytes& operator=(const FieldDataBytes&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBytes& operator=(FieldDataBytes&&) = default;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  EntityBytes<std::byte, Space, DataLayout> entity_bytes(Entity entity,
                                                         const char* file = STK_DEVICE_FILE,
                                                         int line = STK_DEVICE_LINE) const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  EntityBytes<std::byte, Space, DataLayout> entity_bytes(const FastMeshIndex& fmi,
                                                         const char* file = STK_DEVICE_FILE,
                                                         int line = STK_DEVICE_LINE) const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  BucketBytes<std::byte, Space, DataLayout> bucket_bytes(int bucketId,
                                                         const char* file = STK_DEVICE_FILE,
                                                         int line = STK_DEVICE_LINE) const;
};


//==============================================================================
// Host FieldDataBytes
//==============================================================================

template <>
class FieldDataBytes<stk::ngp::HostSpace> : public ConstFieldDataBytes<stk::ngp::HostSpace>
{
public:
  FieldDataBytes();
  FieldDataBytes(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName,
                 const DataTraits& dataTraits, Layout dataLayout);
  KOKKOS_FUNCTION virtual ~FieldDataBytes() override {}

  FieldDataBytes(const FieldDataBytes& fieldData) = default;
  FieldDataBytes(FieldDataBytes&&) = default;
  FieldDataBytes& operator=(const FieldDataBytes&) = default;
  FieldDataBytes& operator=(FieldDataBytes&&) = default;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<std::byte, stk::ngp::HostSpace, DataLayout> entity_bytes(Entity entity,
                                                                       const char* file = STK_HOST_FILE,
                                                                       int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<std::byte, stk::ngp::HostSpace, DataLayout> entity_bytes(const MeshIndex& mi,
                                                                       const char* file = STK_HOST_FILE,
                                                                       int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<std::byte, stk::ngp::HostSpace, DataLayout> entity_bytes(const FastMeshIndex& fmi,
                                                                       const char* file = STK_HOST_FILE,
                                                                       int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  BucketBytes<std::byte, stk::ngp::HostSpace, DataLayout> bucket_bytes(const Bucket& bucket,
                                                                       const char* file = STK_HOST_FILE,
                                                                       int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  BucketBytes<std::byte, stk::ngp::HostSpace, DataLayout> bucket_bytes(int bucketId,
                                                                       const char* file = STK_HOST_FILE,
                                                                       int line = STK_HOST_LINE) const;
};


//==============================================================================
// Device FieldDataBytes definitions
//==============================================================================

template <typename Space>
KOKKOS_FUNCTION
FieldDataBytes<Space>::FieldDataBytes()
  : ConstFieldDataBytes<Space>()
{}

//------------------------------------------------------------------------------
template <typename Space>
FieldDataBytes<Space>::FieldDataBytes(FieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes, Layout dataLayout)
  : ConstFieldDataBytes<Space>(hostFieldBytes, dataLayout)
{}

//------------------------------------------------------------------------------
template <typename Space>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
EntityBytes<std::byte, Space, DataLayout>
FieldDataBytes<Space>::entity_bytes(Entity entity,
                                    const char* file, int line) const
{
  this->check_entity_local_offset(entity.local_offset(), file, line);

  const FastMeshIndex& fmi = this->m_deviceFastMeshIndices[entity.local_offset()];

  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
      this->m_bytesPerScalar;

  return EntityBytes<std::byte, Space, DataLayout>(
        fieldMetaData.m_data + this->m_bytesPerScalar * fmi.bucket_ord,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename Space>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
EntityBytes<std::byte, Space, DataLayout>
FieldDataBytes<Space>::entity_bytes(const FastMeshIndex& fmi,
                                    const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
      this->m_bytesPerScalar;

  return EntityBytes<std::byte, Space, DataLayout>(
        fieldMetaData.m_data + this->m_bytesPerScalar * fmi.bucket_ord,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename Space>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
BucketBytes<std::byte, Space, DataLayout>
FieldDataBytes<Space>::bucket_bytes(int bucketId,
                                    const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[bucketId];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
      this->m_bytesPerScalar;

  return BucketBytes<std::byte, Space, DataLayout>(fieldMetaData.m_data,
                                                   bytesPerEntity,
                                                   this->m_bytesPerScalar,
                                                   fieldMetaData.m_bucketSize,
                                                   fieldMetaData.m_bucketCapacity);
}


//==============================================================================
// Host FieldDataBytes definitions
//==============================================================================

inline
FieldDataBytes<stk::ngp::HostSpace>::FieldDataBytes()
  : ConstFieldDataBytes<stk::ngp::HostSpace>()
{}

//------------------------------------------------------------------------------
inline
FieldDataBytes<stk::ngp::HostSpace>::FieldDataBytes(EntityRank entityRank, Ordinal fieldOrdinal,
                                                    const std::string& fieldName, const DataTraits& dataTraits,
                                                    Layout dataLayout)
  : ConstFieldDataBytes<stk::ngp::HostSpace>(entityRank, fieldOrdinal, fieldName, dataTraits, dataLayout)
{}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Auto>(Entity entity,
                                                                const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Auto>(const MeshIndex& mi,
                                                                const char* file, int line) const
{
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Auto>(const FastMeshIndex& fmi,
                                                                const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  if (m_layout == Layout::Right) {
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * fmi.bucket_ord,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>
FieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Auto>(const Bucket& bucket,
                                                                const char* file, int line) const
{
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  if (m_layout == Layout::Right) {
    return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>
FieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Auto>(int bucketId,
                                                                const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  if (m_layout == Layout::Right) {
    return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, 0);  // Keep compiler happy
  }
}


//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Left>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Left>(Entity entity,
                                                                const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Left>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Left>(const MeshIndex& mi,
                                                                const char* file, int line) const
{
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Left>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Left>(const FastMeshIndex& fmi,
                                                                const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * fmi.bucket_ord,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Left>
FieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Left>(const Bucket& bucket,
                                                                const char* file, int line) const
{
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Left>
FieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Left>(int bucketId,
                                                                const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}


//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Right>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Right>(Entity entity,
                                                                 const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Right>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Right>(const MeshIndex& mi,
                                                                 const char* file, int line) const
{
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Right>
FieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Right>(const FastMeshIndex& fmi,
                                                                 const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytes<std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Right>
FieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Right>(const Bucket& bucket,
                                                                 const char* file, int line) const
{
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Right>
FieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Right>(int bucketId,
                                                                 const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytes<std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize);
}

//==============================================================================

}

#endif // STK_FIELDDATABYTES_HPP
