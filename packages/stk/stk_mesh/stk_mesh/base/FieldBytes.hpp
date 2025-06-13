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

#ifndef STK_FIELDBYTES_HPP
#define STK_FIELDBYTES_HPP

#include "ConstFieldBytes.hpp"

namespace stk::mesh {

//==============================================================================
// Device FieldBytes
//==============================================================================

template <typename MemSpace = stk::ngp::HostMemSpace>
class FieldBytes : public ConstFieldBytes<MemSpace>
{
public:
  KOKKOS_FUNCTION FieldBytes();
  FieldBytes(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName, int bytesPerScalar,
             Layout dataLayout);
  KOKKOS_FUNCTION virtual ~FieldBytes() override {}

  KOKKOS_DEFAULTED_FUNCTION FieldBytes(const FieldBytes& fieldData) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldBytes(FieldBytes&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldBytes& operator=(const FieldBytes&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldBytes& operator=(FieldBytes&&) = default;

  KOKKOS_INLINE_FUNCTION
  EntityBytes<std::byte, MemSpace> entity_bytes(Entity entity,
                                                const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const;

  KOKKOS_INLINE_FUNCTION
  EntityBytes<std::byte, MemSpace> entity_bytes(const FastMeshIndex& fmi,
                                                const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const;

  KOKKOS_INLINE_FUNCTION
  BucketBytes<std::byte, MemSpace> bucket_bytes(int bucketId,
                                                const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const;
};


//==============================================================================
// Host FieldBytes
//==============================================================================

template <>
class FieldBytes<stk::ngp::HostMemSpace> : public ConstFieldBytes<stk::ngp::HostMemSpace>
{
public:
  FieldBytes();
  FieldBytes(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName, const DataTraits& dataTraits,
             Layout dataLayout);
  virtual ~FieldBytes() override = default;

  FieldBytes(const FieldBytes& fieldData) = default;
  FieldBytes(FieldBytes&&) = default;
  FieldBytes& operator=(const FieldBytes&) = default;
  FieldBytes& operator=(FieldBytes&&) = default;

  // These functions will adapt to Layout::Left or Layout::Right automatically, but they are slow
  inline
  EntityBytes<std::byte> entity_bytes(Entity entity,
                                      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  EntityBytes<std::byte> entity_bytes(const MeshIndex& mi,
                                      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  EntityBytes<std::byte> entity_bytes(const FastMeshIndex& fmi,
                                      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  BucketBytes<std::byte> bucket_bytes(const Bucket& bucket,
                                      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  BucketBytes<std::byte> bucket_bytes(int bucketId,
                                      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;


  // These functions will only work correctly if your data is Layout::Left, but they are fast
  inline
  EntityBytesLeft<std::byte> entity_bytes_left(Entity entity,
                                               const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  EntityBytesLeft<std::byte> entity_bytes_left(const MeshIndex& mi,
                                               const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  EntityBytesLeft<std::byte> entity_bytes_left(const FastMeshIndex& fmi,
                                               const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  BucketBytesLeft<std::byte> bucket_bytes_left(const Bucket& bucket,
                                               const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  BucketBytesLeft<std::byte> bucket_bytes_left(int bucketId,
                                               const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;


  // These functions will only work correctly if your data is Layout::Right, but they are fast
  inline
  EntityBytesRight<std::byte> entity_bytes_right(Entity entity,
                                                 const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  EntityBytesRight<std::byte> entity_bytes_right(const MeshIndex& mi,
                                                 const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  EntityBytesRight<std::byte> entity_bytes_right(const FastMeshIndex& fmi,
                                                 const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  BucketBytesRight<std::byte> bucket_bytes_right(const Bucket& bucket,
                                                 const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

  inline
  BucketBytesRight<std::byte> bucket_bytes_right(int bucketId,
                                                 const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const;

};


//==============================================================================
// Device FieldBytes definitions
//==============================================================================

template <typename MemSpace>
KOKKOS_FUNCTION
FieldBytes<MemSpace>::FieldBytes()
  : ConstFieldBytes<MemSpace>()
{}

//------------------------------------------------------------------------------
template <typename MemSpace>
FieldBytes<MemSpace>::FieldBytes(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName,
                                 int bytesPerScalar, Layout dataLayout)
  : ConstFieldBytes<MemSpace>(entityRank, fieldOrdinal, fieldName, bytesPerScalar, dataLayout)
{}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION
EntityBytes<std::byte, MemSpace>
FieldBytes<MemSpace>::entity_bytes(Entity entity,
                                   const char* file, int line) const
{
  this->check_entity_local_offset(entity.local_offset(), file, line);

  const FastMeshIndex& fmi = this->m_deviceFastMeshIndices[entity.local_offset()];

  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
                             this->m_bytesPerScalar;

  return EntityBytes<std::byte, MemSpace>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                       this->m_bytesPerScalar * fmi.bucket_ord),
                                          bytesPerEntity,
                                          this->m_bytesPerScalar,
                                          fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION
EntityBytes<std::byte, MemSpace>
FieldBytes<MemSpace>::entity_bytes(const FastMeshIndex& fmi,
                                   const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
                             this->m_bytesPerScalar;

  return EntityBytes<std::byte, MemSpace>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                       this->m_bytesPerScalar * fmi.bucket_ord),
                                          bytesPerEntity,
                                          this->m_bytesPerScalar,
                                          fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION
BucketBytes<std::byte, MemSpace>
FieldBytes<MemSpace>::bucket_bytes(int bucketId,
                                   const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[bucketId];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
                             this->m_bytesPerScalar;

  return BucketBytes<std::byte, MemSpace>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                          bytesPerEntity,
                                          this->m_bytesPerScalar,
                                          fieldMetaData.m_bucketSize,
                                          fieldMetaData.m_bucketCapacity);
}


//==============================================================================
// Host FieldBytes definitions
//==============================================================================

inline
FieldBytes<stk::ngp::HostMemSpace>::FieldBytes()
  : ConstFieldBytes<stk::ngp::HostMemSpace>()
{}

//------------------------------------------------------------------------------
inline
FieldBytes<stk::ngp::HostMemSpace>::FieldBytes(EntityRank entityRank, Ordinal fieldOrdinal,
                                               const std::string& fieldName, const DataTraits& dataTraits,
                                               Layout dataLayout)
  : ConstFieldBytes<stk::ngp::HostMemSpace>(entityRank, fieldOrdinal, fieldName, dataTraits, dataLayout)
{}

//------------------------------------------------------------------------------
inline
EntityBytes<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes(Entity entity,
                                                 const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                               fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
                                  fieldMetaData.m_bytesPerEntity);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                               this->m_dataTraits->alignment_of * mi.bucket_ordinal),
                                  fieldMetaData.m_bytesPerEntity,
                                  this->m_dataTraits->alignment_of,
                                  fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<std::byte>(nullptr, 0, 0, 0);  // Keep comiler happy
  }
}

//------------------------------------------------------------------------------
inline
EntityBytes<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes(const MeshIndex& mi,
                                                 const char* file, int line) const
{
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                               fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
                                  fieldMetaData.m_bytesPerEntity);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                               this->m_dataTraits->alignment_of * mi.bucket_ordinal),
                                  fieldMetaData.m_bytesPerEntity,
                                  this->m_dataTraits->alignment_of,
                                  fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<std::byte>(nullptr, 0, 0, 0);  // Keep comiler happy
  }
}

//------------------------------------------------------------------------------
inline
EntityBytes<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes(const FastMeshIndex& fmi,
                                                 const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  if (m_layout == Layout::Right) {
    return EntityBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                               fieldMetaData.m_bytesPerEntity * fmi.bucket_ord),
                                  fieldMetaData.m_bytesPerEntity);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                               this->m_dataTraits->alignment_of * fmi.bucket_ord),
                                  fieldMetaData.m_bytesPerEntity,
                                  this->m_dataTraits->alignment_of,
                                  fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<std::byte>(nullptr, 0, 0, 0);  // Keep comiler happy
  }
}

//------------------------------------------------------------------------------
inline
BucketBytes<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::bucket_bytes(const Bucket& bucket,
                                                 const char* file, int line) const
{
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  if (m_layout == Layout::Right) {
    return BucketBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                  fieldMetaData.m_bytesPerEntity,
                                  fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                  fieldMetaData.m_bytesPerEntity,
                                  this->m_dataTraits->alignment_of,
                                  fieldMetaData.m_bucketSize,
                                  fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<std::byte>(nullptr, 0, 0, 0, 0);  // Keep comiler happy
  }
}

//------------------------------------------------------------------------------
inline
BucketBytes<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::bucket_bytes(int bucketId,
                                                 const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  if (m_layout == Layout::Right) {
    return BucketBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                  fieldMetaData.m_bytesPerEntity,
                                  fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                  fieldMetaData.m_bytesPerEntity,
                                  this->m_dataTraits->alignment_of,
                                  fieldMetaData.m_bucketSize,
                                  fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<std::byte>(nullptr, 0, 0, 0, 0);  // Keep comiler happy
  }
}


//------------------------------------------------------------------------------
inline
EntityBytesLeft<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes_left(Entity entity,
                                                      const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytesLeft<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                 this->m_dataTraits->alignment_of * mi.bucket_ordinal),
                                    fieldMetaData.m_bytesPerEntity,
                                    this->m_dataTraits->alignment_of,
                                    fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
inline
EntityBytesLeft<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes_left(const MeshIndex& mi,
                                                      const char* file, int line) const
{
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytesLeft<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                 this->m_dataTraits->alignment_of * mi.bucket_ordinal),
                                    fieldMetaData.m_bytesPerEntity,
                                    this->m_dataTraits->alignment_of,
                                    fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
inline
EntityBytesLeft<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes_left(const FastMeshIndex& fmi,
                                                      const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytesLeft<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                 this->m_dataTraits->alignment_of * fmi.bucket_ord),
                                    fieldMetaData.m_bytesPerEntity,
                                    this->m_dataTraits->alignment_of,
                                    fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
inline
BucketBytesLeft<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::bucket_bytes_left(const Bucket& bucket,
                                                      const char* file, int line) const
{
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytesLeft<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                    fieldMetaData.m_bytesPerEntity,
                                    this->m_dataTraits->alignment_of,
                                    fieldMetaData.m_bucketSize,
                                    fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
inline
BucketBytesLeft<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::bucket_bytes_left(int bucketId,
                                                      const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytesLeft<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                    fieldMetaData.m_bytesPerEntity,
                                    this->m_dataTraits->alignment_of,
                                    fieldMetaData.m_bucketSize,
                                    fieldMetaData.m_bucketCapacity);
}


//------------------------------------------------------------------------------
inline
EntityBytesRight<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes_right(Entity entity,
                                                       const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytesRight<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                  fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
                                     fieldMetaData.m_bytesPerEntity);
}

//------------------------------------------------------------------------------
inline
EntityBytesRight<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes_right(const MeshIndex& mi,
                                                       const char* file, int line) const
{
  this->check_mesh(mi.bucket->mesh(), "Entity", file, line);
  this->check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytesRight<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                  fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal),
                                     fieldMetaData.m_bytesPerEntity);
}

//------------------------------------------------------------------------------
inline
EntityBytesRight<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::entity_bytes_right(const FastMeshIndex& fmi,
                                                       const char* file, int line) const
{
  this->check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytesRight<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data +
                                                                  fieldMetaData.m_bytesPerEntity * fmi.bucket_ord),
                                     fieldMetaData.m_bytesPerEntity);
}

//------------------------------------------------------------------------------
inline
BucketBytesRight<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::bucket_bytes_right(const Bucket& bucket,
                                                       const char* file, int line) const
{
  this->check_mesh(bucket.mesh(), "Bucket", file, line);
  this->check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytesRight<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                     fieldMetaData.m_bytesPerEntity,
                                     fieldMetaData.m_bucketSize);
}

//------------------------------------------------------------------------------
inline
BucketBytesRight<std::byte>
FieldBytes<stk::ngp::HostMemSpace>::bucket_bytes_right(int bucketId,
                                                       const char* file, int line) const
{
  this->check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytesRight<std::byte>(reinterpret_cast<std::byte*>(fieldMetaData.m_data),
                                     fieldMetaData.m_bytesPerEntity,
                                     fieldMetaData.m_bucketSize);
}

//==============================================================================

}

#endif // STK_FIELDBYTES_HPP
