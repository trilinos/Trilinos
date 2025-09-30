#ifndef STK_CONSTFIELDDATABYTES_HPP
#define STK_CONSTFIELDDATABYTES_HPP

#include "stk_mesh/base/FieldDataBase.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/DataTraits.hpp"
#include "stk_mesh/base/EntityBytes.hpp"
#include "stk_mesh/base/BucketBytes.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include <cstddef>
#include <string>

namespace stk::mesh {

//==============================================================================
// Device ConstFieldDataBytes
//==============================================================================

template <typename MemSpace = stk::ngp::HostMemSpace>
class ConstFieldDataBytes : public FieldDataBase
{
public:
  using mem_space = MemSpace;

  KOKKOS_FUNCTION ConstFieldDataBytes();
  ConstFieldDataBytes(ConstFieldDataBytes<stk::ngp::HostMemSpace>* hostFieldBytes, Layout dataLayout);
  KOKKOS_FUNCTION virtual ~ConstFieldDataBytes() override {}

  KOKKOS_INLINE_FUNCTION ConstFieldDataBytes(const ConstFieldDataBytes& other);
  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes(ConstFieldDataBytes&&) = default;
  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes& operator=(const ConstFieldDataBytes&) = default;
  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes& operator=(ConstFieldDataBytes&&) = default;

  KOKKOS_INLINE_FUNCTION Layout data_layout() const;
  KOKKOS_INLINE_FUNCTION EntityRank entity_rank() const;
  KOKKOS_INLINE_FUNCTION Ordinal field_ordinal() const;

  inline BulkData& mesh();
  inline const BulkData& mesh() const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  EntityBytes<const std::byte, MemSpace, DataLayout> entity_bytes(Entity entity,
                                                                  const char* file = STK_DEVICE_FILE,
                                                                  int line = STK_DEVICE_LINE) const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  EntityBytes<const std::byte, MemSpace, DataLayout> entity_bytes(const FastMeshIndex& fmi,
                                                                  const char* file = STK_DEVICE_FILE,
                                                                  int line = STK_DEVICE_LINE) const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  BucketBytes<const std::byte, MemSpace, DataLayout> bucket_bytes(int bucketId,
                                                                  const char* file = STK_DEVICE_FILE,
                                                                  int line = STK_DEVICE_LINE) const;

protected:
  template <typename MemSpace_> friend class DeviceFieldDataManager;
  friend sierra::Fmwk::Region;

  virtual void set_mesh(BulkData* bulkData) override;

  virtual bool needs_update() const override;
  virtual int field_data_synchronized_count() const override;
  virtual void swap_field_data(FieldDataBase& other) override;
  virtual void update_host_bucket_pointers() override;
  virtual void incomplete_swap_field_data(FieldDataBase& other) override;

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void fence(const stk::ngp::ExecSpace&) override {}

  KOKKOS_INLINE_FUNCTION const char* field_name() const;
  inline void modify_field_meta_data();
  inline void update_field_meta_data_mod_count();

#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
  KOKKOS_INLINE_FUNCTION void check_updated_field(const char* file, int line) const;
  KOKKOS_INLINE_FUNCTION void check_entity_local_offset(unsigned localOffset, const char* file, int line) const;
  KOKKOS_INLINE_FUNCTION void check_bucket_id(unsigned bucketId, const char* valuesType, const char* file,
                                              int line) const;
  KOKKOS_INLINE_FUNCTION void check_bucket_ordinal(unsigned bucketId, unsigned bucketOrd, const char* file,
                                                   int line) const;
#else
  KOKKOS_INLINE_FUNCTION void check_updated_field(const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_entity_local_offset(unsigned, const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_bucket_id(unsigned, const char*, const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_bucket_ordinal(unsigned, unsigned, const char*, int) const {}
#endif

  DeviceFieldMetaDataArrayType<MemSpace> m_deviceFieldMetaData;
  MeshIndexType<MemSpace> m_deviceFastMeshIndices;
  DeviceStringType m_fieldName;
  FieldMetaDataModCountType m_fieldMetaDataModCount;
  BulkData* m_hostBulk;
  Ordinal m_ordinal;
  int m_bytesPerScalar;
  int m_fieldDataSynchronizedCount;
  unsigned m_localFieldMetaDataModCount;
  EntityRank m_rank;
  Layout m_layout;
};


//==============================================================================
// Host ConstFieldDataBytes
//==============================================================================

template <>
class ConstFieldDataBytes<stk::ngp::HostMemSpace> : public FieldDataBase
{
public:
  using mem_space = stk::ngp::HostMemSpace;

  ConstFieldDataBytes();
  ConstFieldDataBytes(EntityRank entityRank, Ordinal fieldOrdinal, const std::string& fieldName,
                      const DataTraits& dataTraits, Layout dataLayout);
  virtual ~ConstFieldDataBytes() override = default;

  // The AMD ROCm compiler has an "undefined hidden symbol" link error when this copy constructor is defined
  // below, along with all of the other functions.  Not sure why this one is special.
  inline ConstFieldDataBytes(const ConstFieldDataBytes& other)
    : FieldDataBase(other),
      m_fieldMetaData(other.m_fieldMetaData),
      m_bulk(other.m_bulk),
      m_dataTraits(other.m_dataTraits),
    #if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
      m_fieldName(other.m_fieldName),
    #else
      m_fieldName(),
    #endif
      m_fieldMetaDataModCount(other.m_fieldMetaDataModCount),
      m_ordinal(other.m_ordinal),
      m_fieldDataSynchronizedCount(other.m_fieldDataSynchronizedCount),
      m_localFieldMetaDataModCount(other.m_localFieldMetaDataModCount),
      m_rank(other.m_rank),
      m_layout(other.m_layout)
  {}

  inline ConstFieldDataBytes(ConstFieldDataBytes&&) = default;
  inline ConstFieldDataBytes& operator=(const ConstFieldDataBytes&) = default;
  inline ConstFieldDataBytes& operator=(ConstFieldDataBytes&&) = default;

  inline Layout data_layout() const;
  inline EntityRank entity_rank() const;
  inline Ordinal field_ordinal() const;
  inline const DataTraits& data_traits() const;

  inline BulkData& mesh();
  inline const BulkData& mesh() const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<const std::byte, stk::ngp::HostMemSpace, DataLayout> entity_bytes(Entity entity,
                                                                                const char* file = STK_HOST_FILE,
                                                                                int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<const std::byte, stk::ngp::HostMemSpace, DataLayout> entity_bytes(const MeshIndex& mi,
                                                                                const char* file = STK_HOST_FILE,
                                                                                int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<const std::byte, stk::ngp::HostMemSpace, DataLayout> entity_bytes(const FastMeshIndex& fmi,
                                                                                const char* file = STK_HOST_FILE,
                                                                                int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  BucketBytes<const std::byte, stk::ngp::HostMemSpace, DataLayout> bucket_bytes(const Bucket& bucket,
                                                                                const char* file = STK_HOST_FILE,
                                                                                int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  BucketBytes<const std::byte, stk::ngp::HostMemSpace, DataLayout> bucket_bytes(int bucketId,
                                                                                const char* file = STK_HOST_FILE,
                                                                                int line = STK_HOST_LINE) const;

protected:
  friend FieldBase;
  friend sierra::Fmwk::Region;
  template <typename MemSpace_> friend class ConstFieldDataBytes;

  virtual void set_mesh(BulkData* bulkData) override;

  virtual bool needs_update() const override;
  virtual int field_data_synchronized_count() const override;
  virtual void swap_field_data(FieldDataBase& other) override;
  virtual void update_host_bucket_pointers() override {}
  virtual void incomplete_swap_field_data(FieldDataBase&) override {}

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void fence(const stk::ngp::ExecSpace&) override {}

  inline const char* field_name() const;
  inline void modify_field_meta_data();
  inline void update_field_meta_data_mod_count();

#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
  inline std::string location_string(const char* file, int line) const;
  inline void check_updated_field(const char* file, int line) const;
  inline void check_mesh(const stk::mesh::BulkData& bulk, const char* target, const char* file, int line) const;
  inline void check_rank(stk::mesh::EntityRank entityRank, const char* target, const char* file, int line) const;
  inline void check_bucket_id(unsigned bucketId, const char* valuesType, const char* file, int line) const;
  inline void check_bucket_ordinal(unsigned bucketId, unsigned bucketOrd, const char* file, int line) const;
#else
  inline void check_updated_field(const char*, int) const {}
  inline void check_mesh(const stk::mesh::BulkData&, const char*, const char*, int) const {}
  inline void check_rank(stk::mesh::EntityRank, const char*, const char*, int) const {}
  inline void check_bucket_id(unsigned, const char*, const char*, int) const {}
  inline void check_bucket_ordinal(unsigned, unsigned, const char*, int) const {}
#endif

  FieldMetaDataArrayType m_fieldMetaData;
  BulkData* m_bulk;
  const DataTraits* m_dataTraits;
  HostStringType m_fieldName;
  FieldMetaDataModCountType m_fieldMetaDataModCount;
  Ordinal m_ordinal;
  int m_fieldDataSynchronizedCount;
  unsigned m_localFieldMetaDataModCount;
  EntityRank m_rank;
  Layout m_layout;
};


//==============================================================================
// Device ConstFieldDataBytes definitions
//==============================================================================

template <typename MemSpace>
KOKKOS_FUNCTION
ConstFieldDataBytes<MemSpace>::ConstFieldDataBytes()
  : FieldDataBase(),
    m_hostBulk(nullptr),
    m_ordinal(InvalidOrdinal),
    m_bytesPerScalar(0),
    m_fieldDataSynchronizedCount(0),
    m_localFieldMetaDataModCount(0),
    m_rank(InvalidEntityRank),
    m_layout(Layout::Left)
{
}

//------------------------------------------------------------------------------
template <typename MemSpace>
ConstFieldDataBytes<MemSpace>::ConstFieldDataBytes(ConstFieldDataBytes<stk::ngp::HostMemSpace>* hostFieldBytes,
                                                   Layout dataLayout)
  : FieldDataBase(true),
    m_hostBulk(nullptr),
    m_ordinal(hostFieldBytes->field_ordinal()),
    m_bytesPerScalar(hostFieldBytes->data_traits().alignment_of),
    m_fieldDataSynchronizedCount(0),
    m_localFieldMetaDataModCount(0),
    m_rank(hostFieldBytes->entity_rank()),
    m_layout(dataLayout)
{
  const std::string fieldName(hostFieldBytes->field_name());
  m_fieldName = DeviceStringType(Kokkos::view_alloc(Kokkos::WithoutInitializing, fieldName), fieldName.size()+1);
  std::strcpy(m_fieldName.data(), fieldName.c_str());
  m_fieldMetaDataModCount = hostFieldBytes->m_fieldMetaDataModCount;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION
ConstFieldDataBytes<MemSpace>::ConstFieldDataBytes(const ConstFieldDataBytes& other)
  : FieldDataBase(other),
    m_deviceFieldMetaData(other.m_deviceFieldMetaData),
    m_deviceFastMeshIndices(other.m_deviceFastMeshIndices),
    #if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
    m_fieldName(other.m_fieldName),
    #else
    m_fieldName(),
    #endif
    m_fieldMetaDataModCount(other.m_fieldMetaDataModCount),
    m_hostBulk(other.m_hostBulk),
    m_ordinal(other.m_ordinal),
    m_bytesPerScalar(other.m_bytesPerScalar),
    m_fieldDataSynchronizedCount(other.m_fieldDataSynchronizedCount),
    m_localFieldMetaDataModCount(other.m_localFieldMetaDataModCount),
    m_rank(other.m_rank),
    m_layout(other.m_layout)
{
}

//------------------------------------------------------------------------------
template <typename MemSpace>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
EntityBytes<const std::byte, MemSpace, DataLayout>
ConstFieldDataBytes<MemSpace>::entity_bytes(Entity entity,
                                            const char* file, int line) const
{
  check_entity_local_offset(entity.local_offset(), file, line);

  const FastMeshIndex& fmi = this->m_deviceFastMeshIndices[entity.local_offset()];

  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
                             this->m_bytesPerScalar;

  return EntityBytes<const std::byte, MemSpace, DataLayout>(
        fieldMetaData.m_data + this->m_bytesPerScalar * fmi.bucket_ord,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename MemSpace>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
EntityBytes<const std::byte, MemSpace, DataLayout>
ConstFieldDataBytes<MemSpace>::entity_bytes(const FastMeshIndex& fmi,
                                            const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
                             this->m_bytesPerScalar;

  return EntityBytes<const std::byte, MemSpace, DataLayout>(
        fieldMetaData.m_data + this->m_bytesPerScalar * fmi.bucket_ord,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename MemSpace>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
BucketBytes<const std::byte, MemSpace, DataLayout>
ConstFieldDataBytes<MemSpace>::bucket_bytes(int bucketId,
                                            const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[bucketId];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
                             this->m_bytesPerScalar;

  return BucketBytes<const std::byte, MemSpace, DataLayout>(
        fieldMetaData.m_data,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}


//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION Layout
ConstFieldDataBytes<MemSpace>::data_layout() const
{
  return m_layout;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION EntityRank
ConstFieldDataBytes<MemSpace>::entity_rank() const
{
  return m_rank;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION Ordinal
ConstFieldDataBytes<MemSpace>::field_ordinal() const
{
  return m_ordinal;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION const char*
ConstFieldDataBytes<MemSpace>::field_name() const
{
#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
  return m_fieldName.data();
#else
  return "";
#endif
}

//------------------------------------------------------------------------------
template <typename MemSpace>
void
ConstFieldDataBytes<MemSpace>::modify_field_meta_data()
{
  ++m_fieldMetaDataModCount();
}

//------------------------------------------------------------------------------
template <typename MemSpace>
void
ConstFieldDataBytes<MemSpace>::update_field_meta_data_mod_count()
{
#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
  m_localFieldMetaDataModCount = m_fieldMetaDataModCount();
#endif
}

//------------------------------------------------------------------------------
template <typename MemSpace>
BulkData&
ConstFieldDataBytes<MemSpace>::mesh()
{
  STK_ThrowAssert(m_hostBulk != nullptr);
  return *m_hostBulk;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
const BulkData&
ConstFieldDataBytes<MemSpace>::mesh() const
{
  STK_ThrowAssert(m_hostBulk != nullptr);
  return *m_hostBulk;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
void
ConstFieldDataBytes<MemSpace>::set_mesh(BulkData* bulkData)
{
  m_hostBulk = bulkData;
  m_fieldDataSynchronizedCount = 0;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
bool
ConstFieldDataBytes<MemSpace>::needs_update() const
{
#ifndef NDEBUG
  const int maxValidSyncCount = static_cast<int>(mesh().synchronized_count()+1);
  STK_ThrowAssertMsg(m_fieldDataSynchronizedCount <= maxValidSyncCount,
                     "Invalid sync state detected for Field: " << field_name() << ": field-sync-count (" <<
                     m_fieldDataSynchronizedCount << ") shouldn't be greater than mesh-sync-count (" <<
                     mesh().synchronized_count() << ")");
#endif
  return m_fieldDataSynchronizedCount != static_cast<int>(mesh().synchronized_count());
}

//------------------------------------------------------------------------------
template <typename MemSpace>
int
ConstFieldDataBytes<MemSpace>::field_data_synchronized_count() const
{
  return m_fieldDataSynchronizedCount;
}

//------------------------------------------------------------------------------
template <typename MemSpace>
void
ConstFieldDataBytes<MemSpace>::swap_field_data(FieldDataBase& other)
{
  ConstFieldDataBytes<MemSpace>* otherFieldBytes = dynamic_cast<ConstFieldDataBytes<MemSpace>*>(&other);
  STK_ThrowRequireMsg(otherFieldBytes != nullptr,
                      "ConstFieldDataBytes::swap_field_data() called with an imcompatible ConstFieldDataBytes object.");

  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<MemSpace>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  deviceFieldDataManager->swap_field_data(this->field_ordinal(), otherFieldBytes->field_ordinal());

  deviceFieldDataManager->set_device_field_meta_data(*this);
  deviceFieldDataManager->set_device_field_meta_data(other);
}

//------------------------------------------------------------------------------
template <typename MemSpace>
void
ConstFieldDataBytes<MemSpace>::update_host_bucket_pointers()
{
  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<MemSpace>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  deviceFieldDataManager->update_host_bucket_pointers(this->field_ordinal());
}

//------------------------------------------------------------------------------
template <typename MemSpace>
void
ConstFieldDataBytes<MemSpace>::incomplete_swap_field_data(FieldDataBase& other)
{
  ConstFieldDataBytes<MemSpace>* otherFieldBytes = dynamic_cast<ConstFieldDataBytes<MemSpace>*>(&other);
  STK_ThrowRequireMsg(otherFieldBytes != nullptr,
                      "ConstFieldDataBytes::incomplete_swap_field_data() called with an imcompatible "
                      "ConstFieldDataBytes object.");

  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<MemSpace>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  deviceFieldDataManager->swap_field_data(this->field_ordinal(), otherFieldBytes->field_ordinal());

  // Reset the host bucket pointers after the rotation, to mimic the incomplete behavior
  // of only rotating the device field data and leaving the host bucket pointers alone.
  // This is only called (in Sierra) after neglecting to rotate the host-side pointers.
  deviceFieldDataManager->update_host_bucket_pointers(this->field_ordinal());
  deviceFieldDataManager->update_host_bucket_pointers( otherFieldBytes->field_ordinal());

  deviceFieldDataManager->set_device_field_meta_data(*this);
  deviceFieldDataManager->set_device_field_meta_data(other);
}

#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<MemSpace>::check_updated_field(const char* file, int line) const
{
  if (this->m_localFieldMetaDataModCount != this->m_fieldMetaDataModCount()) {
    if (line == -1) {
      printf("Error: Accessing out-of-date FieldData after a mesh modification for Field '%s'.  "
             "please re-acquire this FieldData instance.", this->field_name());
    }
    else {
      printf("Error: %s:%i: Accessing out-of-date FieldData after a mesh modification for Field '%s'.  "
             "please re-acquire this FieldData instance.", file, line, this->field_name());
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<MemSpace>::check_entity_local_offset(unsigned localOffset, const char* file, int line) const
{
  if (localOffset >= this->m_deviceFastMeshIndices.extent(0)) {
    if (line == -1) {
      printf("Error: Called FieldData::entity_values() for Field '%s' with an out-of-bounds Entity local "
             "offset (%u) for a FastMeshIndex array with extent %lu.\n", this->field_name(), localOffset,
             this->m_deviceFastMeshIndices.extent(0));
    }
    else {
      printf("Error: %s:%i: Called FieldData::entity_values() for Field '%s' with an out-of-bounds Entity local "
             "offset (%u) for a FastMeshIndex array with extent %lu.\n", file, line, this->field_name(), localOffset,
             this->m_deviceFastMeshIndices.extent(0));
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<MemSpace>::check_bucket_id(unsigned bucketId, const char* valuesType, const char* file,
                                               int line) const
{
  if (bucketId >= this->m_deviceFieldMetaData.extent(0)) {
    if (line == -1) {
      printf("Error: Called FieldData::%s_values() for Field '%s' with an out-of-bounds Bucket ID (%u) for a "
             "DeviceFieldMetaData array with extent %lu.\n", valuesType, this->field_name(), bucketId,
             this->m_deviceFieldMetaData.extent(0));
    }
    else {
      printf("Error: %s:%i: Called FieldData::%s_values() for Field '%s' with an out-of-bounds Bucket ID (%u) for a "
             "DeviceFieldMetaData array with extent %lu.\n", file, line, valuesType, this->field_name(), bucketId,
             this->m_deviceFieldMetaData.extent(0));
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}

//------------------------------------------------------------------------------
template <typename MemSpace>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<MemSpace>::check_bucket_ordinal(unsigned bucketId, unsigned bucketOrd, const char* file,
                                                    int line) const
{
  // Only trip if we're referencing an out-of-bounds Entity in a Bucket where the Field is registered,
  // because accessing an EntityValues where the Field isn't valid should be allowed.  This allows users
  // to query EntityValues::is_field_defined() inside a loop.
  if ((bucketOrd >= static_cast<unsigned>(this->m_deviceFieldMetaData[bucketId].m_bucketSize)) &&
      (this->m_deviceFieldMetaData[bucketId].m_bucketSize > 0)) {
    if (line == -1) {
      printf("Error: Called FieldData::entity_values() for Field '%s' with an out-of-bounds Bucket ordinal (%u) "
             "for Bucket %u with size %i.\n", this->field_name(), bucketOrd, bucketId,
             this->m_deviceFieldMetaData[bucketId].m_bucketSize);
    }
    else {
      printf("Error: %s:%i: Called FieldData::entity_values() for Field '%s' with an out-of-bounds Bucket ordinal (%u) "
             "for Bucket %u with size %i.\n", file, line, this->field_name(), bucketOrd, bucketId,
             this->m_deviceFieldMetaData[bucketId].m_bucketSize);
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}
#endif  // !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)


//==============================================================================
// Host ConstFieldDataBytes definitions
//==============================================================================

inline
ConstFieldDataBytes<stk::ngp::HostMemSpace>::ConstFieldDataBytes()
  : FieldDataBase(),
    m_bulk(nullptr),
    m_dataTraits(&stk::mesh::data_traits<void>()),
    m_ordinal(InvalidOrdinal),
    m_fieldDataSynchronizedCount(0),
    m_localFieldMetaDataModCount(0),
    m_rank(InvalidEntityRank),
    m_layout(Layout::Right)
{
}

//------------------------------------------------------------------------------
inline
ConstFieldDataBytes<stk::ngp::HostMemSpace>::ConstFieldDataBytes(EntityRank entityRank, Ordinal fieldOrdinal,
                                                                 [[maybe_unused]] const std::string& fieldName,
                                                                 const DataTraits& dataTraits, Layout dataLayout)
  : FieldDataBase(true),
    m_bulk(nullptr),
    m_dataTraits(&dataTraits),
    m_ordinal(fieldOrdinal),
    m_fieldDataSynchronizedCount(0),
    m_localFieldMetaDataModCount(0),
    m_rank(entityRank),
    m_layout(dataLayout)
{
  m_fieldName = HostStringType(fieldName, fieldName.size()+1);
  std::strcpy(m_fieldName.data(), fieldName.c_str());
  m_fieldMetaDataModCount = FieldMetaDataModCountType("FieldMetaDataModCount");
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Auto>(Entity entity,
                                                                        const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Auto>(const MeshIndex& mi,
                                                                        const char* file, int line) const
{
  check_mesh(mi.bucket->mesh(), "Entity", file, line);
  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Auto>(const FastMeshIndex& fmi,
                                                                        const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  if (m_layout == Layout::Right) {
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * fmi.bucket_ord,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::bucket_bytes<Layout::Auto>(const Bucket& bucket,
                                                                        const char* file, int line) const
{
  check_mesh(bucket.mesh(), "Bucket", file, line);
  check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  if (m_layout == Layout::Right) {
    return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(nullptr, 0, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::bucket_bytes<Layout::Auto>(int bucketId,
                                                                        const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  if (m_layout == Layout::Right) {
    return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Auto>(nullptr, 0, 0, 0, 0);  // Keep compiler happy
  }
}


//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Left>(Entity entity,
                                                                        const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Left>(const MeshIndex& mi,
                                                                        const char* file, int line) const
{
  check_mesh(mi.bucket->mesh(), "Entity", file, line);
  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Left>(const FastMeshIndex& fmi,
                                                                        const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * fmi.bucket_ord,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::bucket_bytes<Layout::Left>(const Bucket& bucket,
                                                                        const char* file, int line) const
{
  check_mesh(bucket.mesh(), "Bucket", file, line);
  check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::bucket_bytes<Layout::Left>(int bucketId,
                                                                        const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Left>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}


//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Right>(Entity entity,
                                                                         const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Right>(const MeshIndex& mi,
                                                                         const char* file, int line) const
{
  check_mesh(mi.bucket->mesh(), "Entity", file, line);
  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_bytes<Layout::Right>(const FastMeshIndex& fmi,
                                                                         const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::bucket_bytes<Layout::Right>(const Bucket& bucket,
                                                                         const char* file, int line) const
{
  check_mesh(bucket.mesh(), "Bucket", file, line);
  check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostMemSpace>::bucket_bytes<Layout::Right>(int bucketId,
                                                                         const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytes<const std::byte, stk::ngp::HostMemSpace, Layout::Right>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize);
}


//------------------------------------------------------------------------------
inline Layout
ConstFieldDataBytes<stk::ngp::HostMemSpace>::data_layout() const
{
  return m_layout;
}

//------------------------------------------------------------------------------
inline EntityRank
ConstFieldDataBytes<stk::ngp::HostMemSpace>::entity_rank() const
{
  return m_rank;
}

//------------------------------------------------------------------------------
inline Ordinal
ConstFieldDataBytes<stk::ngp::HostMemSpace>::field_ordinal() const
{
  return m_ordinal;
}

//------------------------------------------------------------------------------
inline const char*
ConstFieldDataBytes<stk::ngp::HostMemSpace>::field_name() const
{
#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
  return m_fieldName.data();
#else
  return "";
#endif
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::modify_field_meta_data()
{
  ++m_fieldMetaDataModCount();
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::update_field_meta_data_mod_count()
{
#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
  m_localFieldMetaDataModCount = m_fieldMetaDataModCount();
#endif
}


//------------------------------------------------------------------------------
inline const DataTraits&
ConstFieldDataBytes<stk::ngp::HostMemSpace>::data_traits() const
{
  return *m_dataTraits;
}

//------------------------------------------------------------------------------
inline BulkData&
ConstFieldDataBytes<stk::ngp::HostMemSpace>::mesh()
{
  STK_ThrowAssert(m_bulk != nullptr);
  return *m_bulk;
}

//------------------------------------------------------------------------------
inline const BulkData&
ConstFieldDataBytes<stk::ngp::HostMemSpace>::mesh() const
{
  STK_ThrowAssert(m_bulk != nullptr);
  return *m_bulk;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::set_mesh(BulkData* bulkData)
{
  m_bulk = bulkData;
  m_fieldDataSynchronizedCount = 0;
}

//------------------------------------------------------------------------------
inline bool
ConstFieldDataBytes<stk::ngp::HostMemSpace>::needs_update() const
{
#ifndef NDEBUG
  const int maxValidSyncCount = static_cast<int>(mesh().synchronized_count()+1);
  STK_ThrowAssertMsg(m_fieldDataSynchronizedCount <= maxValidSyncCount,
                     "Invalid sync state detected for Field: " << field_name()
                     << ": field-sync-count (" << m_fieldDataSynchronizedCount
                     << ") shouldn't be greater than mesh-sync-count ("
                     << mesh().synchronized_count() << ")");
#endif
  return m_fieldDataSynchronizedCount != static_cast<int>(mesh().synchronized_count());
}

//------------------------------------------------------------------------------
inline int
ConstFieldDataBytes<stk::ngp::HostMemSpace>::field_data_synchronized_count() const
{
  return m_fieldDataSynchronizedCount;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::swap_field_data(FieldDataBase& other)
{
  ConstFieldDataBytes<stk::ngp::HostMemSpace>* otherFieldBytes =
      dynamic_cast<ConstFieldDataBytes<stk::ngp::HostMemSpace>*>(&other);
  STK_ThrowRequireMsg(otherFieldBytes != nullptr,
                      "ConstFieldDataBytes::swap_field_data() called with an imcompatible ConstFieldDataBytes object.");

  std::swap(this->m_fieldMetaData, otherFieldBytes->m_fieldMetaData);
}

#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)

//------------------------------------------------------------------------------
inline std::string
ConstFieldDataBytes<stk::ngp::HostMemSpace>::location_string(const char* file, int line) const
{
  if (line != -1) {
    std::string fileName(file);
    std::size_t pathDelimeter = fileName.find_last_of("/");
    if (pathDelimeter < fileName.size()) {
      fileName = fileName.substr(pathDelimeter+1);
    }
    return fileName + ":" + std::to_string(line) + ": ";
  }
  else {
    return "";
  }
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::check_updated_field(const char* file, int line) const
{
  STK_ThrowRequireMsg(m_localFieldMetaDataModCount == m_fieldMetaDataModCount(),
                      location_string(file, line) << "Accessing out-of-date FieldData after a mesh modification "
                      "for Field '" << field_name() << "'.  Please re-acquire this FieldData instance.");
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::check_mesh(const stk::mesh::BulkData& bulk, const char* target,
                                                        const char* file, int line) const
{
  STK_ThrowRequireMsg(&bulk == &mesh(),
                      location_string(file, line) << "Accessing " << target << " from a different mesh for Field '" <<
                      field_name() << "'.");
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::check_rank(stk::mesh::EntityRank targetRank, const char* target,
                                                        const char* file, int line) const
{
  STK_ThrowRequireMsg(entity_rank() == targetRank,
                      location_string(file, line) << "Accessing " << target << " with rank " << targetRank <<
                      " for Field '" << field_name() << "' with rank " << entity_rank() <<
                      ".  Are the Field and " << target << " from the same mesh?");
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::check_bucket_id(unsigned bucketId, const char* valuesType,
                                                             const char* file, int line) const
{
  STK_ThrowRequireMsg(bucketId < m_fieldMetaData.extent(0),
                      location_string(file, line) << "Called FieldData::" << valuesType << "_values() for Field '" <<
                      field_name() << "' with an out-of-bounds Bucket ID (" << bucketId <<
                      ") for a FieldMetaData array with extent " << m_fieldMetaData.extent(0));
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostMemSpace>::check_bucket_ordinal(unsigned bucketId, unsigned bucketOrd,
                                                                  const char* file, int line) const
{
  // Only trip if we're referencing an out-of-bounds Entity in a Bucket where the Field is registered,
  // because accessing an EntityValues where the Field isn't valid should be allowed.  This allows users
  // to query EntityValues::is_field_defined() inside a loop.
  STK_ThrowRequireMsg((bucketOrd < static_cast<unsigned>(m_fieldMetaData[bucketId].m_bucketSize)) ||
                      (m_fieldMetaData[bucketId].m_bucketSize == 0),
                      location_string(file, line) << "Called FieldData::entity_values() for Field '" <<
                      field_name() << "' with an out-of-bounds Bucket ordinal (" << bucketOrd << ") for Bucket " <<
                      bucketId << " with size " << m_fieldMetaData[bucketId].m_bucketSize);
}

#endif  // !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK

//==============================================================================

}

#endif // STK_CONSTFIELDDATABYTES_HPP
