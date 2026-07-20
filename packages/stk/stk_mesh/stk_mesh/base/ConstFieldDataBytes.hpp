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

template <typename Space = stk::ngp::HostSpace>
class ConstFieldDataBytes : public FieldDataBase
{
public:
  using space = Space;
  using exec_space = typename Space::exec_space;
  using mem_space = typename Space::mem_space;

  KOKKOS_FUNCTION ConstFieldDataBytes();
  ConstFieldDataBytes(ConstFieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes, Layout dataLayout,
                      FieldDataCopyTracking* copyTracking);
  KOKKOS_FUNCTION virtual ~ConstFieldDataBytes() override {}

  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes(const ConstFieldDataBytes&) = default;
  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes(ConstFieldDataBytes&&) = default;
  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes& operator=(const ConstFieldDataBytes&) = default;
  KOKKOS_DEFAULTED_FUNCTION ConstFieldDataBytes& operator=(ConstFieldDataBytes&&) = default;

  inline BulkData& mesh();
  inline const BulkData& mesh() const;
  KOKKOS_INLINE_FUNCTION const char* field_name() const;
  KOKKOS_INLINE_FUNCTION int num_buckets() const;
  virtual bool needs_update() const override;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  EntityBytes<const std::byte, Space, DataLayout> entity_bytes(Entity entity,
                                                               const char* file = STK_DEVICE_FILE,
                                                               int line = STK_DEVICE_LINE) const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  EntityBytes<const std::byte, Space, DataLayout> entity_bytes(const FastMeshIndex& fmi,
                                                               const char* file = STK_DEVICE_FILE,
                                                               int line = STK_DEVICE_LINE) const;

  template <Layout DataLayout = Layout::Left>
  KOKKOS_INLINE_FUNCTION
  BucketBytes<const std::byte, Space, DataLayout> bucket_bytes(int bucketId,
                                                               const char* file = STK_DEVICE_FILE,
                                                               int line = STK_DEVICE_LINE) const;

protected:
  friend FieldBase;
  template <typename MemSpace_> friend class DeviceFieldDataManager;
  template <typename MemSpace_> friend class impl::DeviceBucketRepository;
  friend sierra::Fmwk::Region;

  virtual void set_mesh(BulkData* bulkData) override;

  virtual void swap_field_data(FieldDataBase& other) override;
  virtual void update_host_bucket_pointers() override;
  virtual void incomplete_swap_field_data(FieldDataBase& other) override;

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace&, Layout, bool) override {}
  virtual void fence(const stk::ngp::ExecSpace&) override {}

  inline void set_up_to_date();
  inline void set_fast_mesh_indices(FastMeshIndex* fastMeshIndicesPtr, unsigned fastMeshIndicesSize);

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

  DeviceFieldMetaData* m_deviceFieldMetaData;         //  8 : pointer
  FastMeshIndex* m_deviceFastMeshIndices;             //  8 : pointer
  const char* m_fieldName;                            //  8 : pointer
  unsigned* m_fieldMetaDataModCount;                  //  8 : pointer
  BulkData* m_hostBulk;                               //  8 : pointer
  unsigned m_localOffsetExtent;                       //  4 : unsigned
  int m_numBuckets;                                   //  4 : int
  int m_bytesPerScalar;                               //  4 : int
  unsigned m_localFieldMetaDataModCount;              //  4 : unsigned
};


//==============================================================================
// Host ConstFieldDataBytes
//==============================================================================

template <>
class ConstFieldDataBytes<stk::ngp::HostSpace> : public FieldDataBase
{
public:
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;

  ConstFieldDataBytes();
  ConstFieldDataBytes(EntityRank entityRank, Ordinal fieldOrdinal, const DataTraits& dataTraits, Layout dataLayout);
  KOKKOS_FUNCTION virtual ~ConstFieldDataBytes() override {}

  inline ConstFieldDataBytes(const ConstFieldDataBytes&) = default;
  inline ConstFieldDataBytes(ConstFieldDataBytes&&) = default;
  inline ConstFieldDataBytes& operator=(const ConstFieldDataBytes&) = default;
  inline ConstFieldDataBytes& operator=(ConstFieldDataBytes&&) = default;

  inline BulkData& mesh();
  inline const BulkData& mesh() const;
  inline const DataTraits& data_traits() const;
  inline const char* field_name() const;
  inline int num_buckets() const;
  virtual bool needs_update() const override;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<const std::byte, stk::ngp::HostSpace, DataLayout> entity_bytes(Entity entity,
                                                                             const char* file = STK_HOST_FILE,
                                                                             int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<const std::byte, stk::ngp::HostSpace, DataLayout> entity_bytes(const MeshIndex& mi,
                                                                             const char* file = STK_HOST_FILE,
                                                                             int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  EntityBytes<const std::byte, stk::ngp::HostSpace, DataLayout> entity_bytes(const FastMeshIndex& fmi,
                                                                             const char* file = STK_HOST_FILE,
                                                                             int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  BucketBytes<const std::byte, stk::ngp::HostSpace, DataLayout> bucket_bytes(const Bucket& bucket,
                                                                             const char* file = STK_HOST_FILE,
                                                                             int line = STK_HOST_LINE) const;

  template <Layout DataLayout = Layout::Auto>
  inline
  BucketBytes<const std::byte, stk::ngp::HostSpace, DataLayout> bucket_bytes(int bucketId,
                                                                             const char* file = STK_HOST_FILE,
                                                                             int line = STK_HOST_LINE) const;

protected:
  friend FieldBase;
  friend sierra::Fmwk::Region;
  friend FieldDataManager;
  template <typename MemSpace_> friend class ConstFieldDataBytes;
  template <typename MemSpace_> friend class DeviceMeshT;
  template <typename MemSpace_> friend class DeviceFieldDataManager;
  template <typename MemSpace_> friend class impl::DeviceBucketRepository;

  virtual void set_mesh(BulkData* bulkData) override;

  virtual void swap_field_data(FieldDataBase& other) override;
  virtual void update_host_bucket_pointers() override {}
  virtual void incomplete_swap_field_data(FieldDataBase&) override {}

  virtual void sync_to_host(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void sync_to_device(const stk::ngp::ExecSpace&, Layout) override {}
  virtual void update(const stk::ngp::ExecSpace&, Layout, bool) override {}
  virtual void fence(const stk::ngp::ExecSpace&) override {}

  inline void set_up_to_date();
  inline void set_field_name(const char* fieldName);
  inline void set_field_meta_data_mod_count_pointer(unsigned* fieldMetaDataModCountPtr);
  inline void set_field_meta_data_mod_count(unsigned fieldMetaDataModCount);
  inline unsigned field_meta_data_mod_count();
  inline void set_fast_mesh_indices(FastMeshIndex*, unsigned) {}

#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
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

  FieldMetaData* m_fieldMetaData;                      //  8 : pointer
  BulkData* m_bulk;                                    //  8 : pointer
  const DataTraits* m_dataTraits;                      //  8 : pointer
  const char* m_fieldName;                             //  8 : pointer
  unsigned* m_fieldMetaDataModCount;                   //  8 : pointer
  int m_numBuckets;                                    //  4 : int
  unsigned m_localFieldMetaDataModCount;               //  4 : unsigned
};


//==============================================================================
// Device ConstFieldDataBytes definitions
//==============================================================================

template <typename Space>
KOKKOS_FUNCTION
ConstFieldDataBytes<Space>::ConstFieldDataBytes()
  : FieldDataBase(),
    m_deviceFieldMetaData(nullptr),
    m_deviceFastMeshIndices(nullptr),
    m_fieldName(nullptr),
    m_fieldMetaDataModCount(nullptr),
    m_hostBulk(nullptr),
    m_localOffsetExtent(0),
    m_numBuckets(0),
    m_bytesPerScalar(0),
    m_localFieldMetaDataModCount(0)
{
}

//------------------------------------------------------------------------------
template <typename Space>
ConstFieldDataBytes<Space>::ConstFieldDataBytes(ConstFieldDataBytes<stk::ngp::HostSpace>* hostFieldBytes,
                                                Layout dataLayout, FieldDataCopyTracking* copyTracking)
  : FieldDataBase(copyTracking, hostFieldBytes->field_ordinal(), hostFieldBytes->entity_rank(), dataLayout),
    m_deviceFieldMetaData(nullptr),
    m_deviceFastMeshIndices(nullptr),
    m_fieldName(hostFieldBytes->m_fieldName),
    m_fieldMetaDataModCount(hostFieldBytes->m_fieldMetaDataModCount),
    m_hostBulk(nullptr),
    m_localOffsetExtent(0),
    m_numBuckets(0),
    m_bytesPerScalar(hostFieldBytes->data_traits().alignment_of),
    m_localFieldMetaDataModCount(0)
{
}

//------------------------------------------------------------------------------
template <typename Space>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
EntityBytes<const std::byte, Space, DataLayout>
ConstFieldDataBytes<Space>::entity_bytes(Entity entity,
                                         const char* file, int line) const
{
  check_entity_local_offset(entity.local_offset(), file, line);

  const FastMeshIndex& fmi = this->m_deviceFastMeshIndices[entity.local_offset()];

  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
      this->m_bytesPerScalar;

  return EntityBytes<const std::byte, Space, DataLayout>(
        fieldMetaData.m_data + this->m_bytesPerScalar * fmi.bucket_ord,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename Space>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
EntityBytes<const std::byte, Space, DataLayout>
ConstFieldDataBytes<Space>::entity_bytes(const FastMeshIndex& fmi,
                                         const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[fmi.bucket_id];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
      this->m_bytesPerScalar;

  return EntityBytes<const std::byte, Space, DataLayout>(
        fieldMetaData.m_data + this->m_bytesPerScalar * fmi.bucket_ord,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <typename Space>
template <Layout DataLayout>
KOKKOS_INLINE_FUNCTION
BucketBytes<const std::byte, Space, DataLayout>
ConstFieldDataBytes<Space>::bucket_bytes(int bucketId,
                                         const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const DeviceFieldMetaData& fieldMetaData = this->m_deviceFieldMetaData[bucketId];
  const int bytesPerEntity = fieldMetaData.m_numComponentsPerEntity * fieldMetaData.m_numCopiesPerEntity *
      this->m_bytesPerScalar;

  return BucketBytes<const std::byte, Space, DataLayout>(
        fieldMetaData.m_data,
        bytesPerEntity,
        this->m_bytesPerScalar,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}


//------------------------------------------------------------------------------
template <typename Space>
KOKKOS_INLINE_FUNCTION const char*
ConstFieldDataBytes<Space>::field_name() const
{
  return m_fieldName;
}

template <typename Space>
KOKKOS_INLINE_FUNCTION int
ConstFieldDataBytes<Space>::num_buckets() const
{
  return m_numBuckets;
}

//------------------------------------------------------------------------------
template <typename Space>
void
ConstFieldDataBytes<Space>::set_up_to_date()
{
  m_localFieldMetaDataModCount = *m_fieldMetaDataModCount;
}

//------------------------------------------------------------------------------
template <typename Space>
void
ConstFieldDataBytes<Space>::set_fast_mesh_indices(FastMeshIndex* fastMeshIndicesPtr, unsigned fastMeshIndicesSize)
{
  m_deviceFastMeshIndices = fastMeshIndicesPtr;
  m_localOffsetExtent = fastMeshIndicesSize;
}

//------------------------------------------------------------------------------
template <typename Space>
BulkData&
ConstFieldDataBytes<Space>::mesh()
{
  STK_ThrowAssert(m_hostBulk != nullptr);
  return *m_hostBulk;
}

//------------------------------------------------------------------------------
template <typename Space>
const BulkData&
ConstFieldDataBytes<Space>::mesh() const
{
  STK_ThrowAssert(m_hostBulk != nullptr);
  return *m_hostBulk;
}

//------------------------------------------------------------------------------
template <typename Space>
void
ConstFieldDataBytes<Space>::set_mesh(BulkData* bulkData)
{
  m_hostBulk = bulkData;
  m_localFieldMetaDataModCount = 0;
}

//------------------------------------------------------------------------------
template <typename Space>
bool
ConstFieldDataBytes<Space>::needs_update() const
{
  return (m_localFieldMetaDataModCount != *m_fieldMetaDataModCount || mesh().in_modifiable_state());
}

//------------------------------------------------------------------------------
[[maybe_unused]]
static
void swap_all_meta_data_pointers_on_device(DeviceFieldMetaData* deviceFieldMetaDataA,
                                           DeviceFieldMetaData* deviceFieldMetaDataB, unsigned numBuckets)
{
  Kokkos::parallel_for(numBuckets,
                       KOKKOS_LAMBDA(unsigned bucketIdx) {
                         std::byte* tmpData = deviceFieldMetaDataA[bucketIdx].m_data;
                         deviceFieldMetaDataA[bucketIdx].m_data = deviceFieldMetaDataB[bucketIdx].m_data;
                         deviceFieldMetaDataB[bucketIdx].m_data = tmpData;

                         std::byte* tmpHostData = deviceFieldMetaDataA[bucketIdx].m_hostData;
                         deviceFieldMetaDataA[bucketIdx].m_hostData = deviceFieldMetaDataB[bucketIdx].m_hostData;
                         deviceFieldMetaDataB[bucketIdx].m_hostData = tmpHostData;
                       }
                      );
}

//------------------------------------------------------------------------------
[[maybe_unused]]
static
void swap_device_meta_data_pointers_on_device(DeviceFieldMetaData* deviceFieldMetaDataA,
                                              DeviceFieldMetaData* deviceFieldMetaDataB, unsigned numBuckets)
{
  Kokkos::parallel_for(numBuckets,
                       KOKKOS_LAMBDA(unsigned bucketIdx) {
                         std::byte* tmpData = deviceFieldMetaDataA[bucketIdx].m_data;
                         deviceFieldMetaDataA[bucketIdx].m_data = deviceFieldMetaDataB[bucketIdx].m_data;
                         deviceFieldMetaDataB[bucketIdx].m_data = tmpData;
                       }
                      );
}

//------------------------------------------------------------------------------
template <typename Space>
void
ConstFieldDataBytes<Space>::swap_field_data(FieldDataBase& other)
{
  ConstFieldDataBytes<Space>* otherFieldBytes = dynamic_cast<ConstFieldDataBytes<Space>*>(&other);
  STK_ThrowRequireMsg(otherFieldBytes != nullptr,
                      "ConstFieldDataBytes::swap_field_data() called with an imcompatible ConstFieldDataBytes object.");

  swap_all_meta_data_pointers_on_device(this->m_deviceFieldMetaData, otherFieldBytes->m_deviceFieldMetaData,
                                        this->m_numBuckets);

  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<Space>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  deviceFieldDataManager->swap_host_cache_all_meta_data_pointers(this->field_ordinal(),
                                                                 otherFieldBytes->field_ordinal());

  std::swap(this->m_fieldMetaDataModCount, otherFieldBytes->m_fieldMetaDataModCount);
  std::swap(this->m_localFieldMetaDataModCount, otherFieldBytes->m_localFieldMetaDataModCount);
}

//------------------------------------------------------------------------------
template <typename Space>
void
ConstFieldDataBytes<Space>::update_host_bucket_pointers()
{
  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<Space>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  deviceFieldDataManager->update_host_bucket_pointers(this->field_ordinal());
}

//------------------------------------------------------------------------------
template <typename Space>
void
ConstFieldDataBytes<Space>::incomplete_swap_field_data(FieldDataBase& other)
{
  ConstFieldDataBytes<Space>* otherFieldBytes = dynamic_cast<ConstFieldDataBytes<Space>*>(&other);
  STK_ThrowRequireMsg(otherFieldBytes != nullptr,
                      "ConstFieldDataBytes::incomplete_swap_field_data() called with an imcompatible "
                      "ConstFieldDataBytes object.");

  swap_device_meta_data_pointers_on_device(this->m_deviceFieldMetaData, otherFieldBytes->m_deviceFieldMetaData,
                                           this->m_numBuckets);

  DeviceFieldDataManagerBase* deviceFieldDataManager = impl::get_device_field_data_manager<Space>(this->mesh());
  STK_ThrowRequire(deviceFieldDataManager != nullptr);

  deviceFieldDataManager->swap_host_cache_device_meta_data_pointers(this->field_ordinal(),
                                                                    otherFieldBytes->field_ordinal());

  std::swap(this->m_fieldMetaDataModCount, otherFieldBytes->m_fieldMetaDataModCount);
  std::swap(this->m_localFieldMetaDataModCount, otherFieldBytes->m_localFieldMetaDataModCount);
}


#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)
//------------------------------------------------------------------------------
template <typename Space>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<Space>::check_updated_field(const char* file, int line) const
{
  if (this->m_localFieldMetaDataModCount != *this->m_fieldMetaDataModCount) {
    if (line == 0) {
      printf("Error: Accessing out-of-date FieldData after a host or device mesh modification for Field '%s'.  "
             "Please re-acquire this FieldData instance.", this->field_name());
    }
    else {
      printf("Error: %s:%i: Accessing out-of-date FieldData after a host or device mesh modification for Field '%s'.  "
             "Please re-acquire this FieldData instance.", file, line, this->field_name());
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}

//------------------------------------------------------------------------------
template <typename Space>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<Space>::check_entity_local_offset(unsigned localOffset, const char* file, int line) const
{
  if (localOffset >= this->m_localOffsetExtent) {
    if (line == 0) {
      printf("Error: Called FieldData::entity_values() for Field '%s' with an out-of-bounds Entity local "
             "offset (%u) for a FastMeshIndex array with extent %u.\n", this->field_name(), localOffset,
             this->m_localOffsetExtent);
    }
    else {
      printf("Error: %s:%i: Called FieldData::entity_values() for Field '%s' with an out-of-bounds Entity local "
             "offset (%u) for a FastMeshIndex array with extent %u.\n", file, line, this->field_name(), localOffset,
             this->m_localOffsetExtent);
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}

//------------------------------------------------------------------------------
template <typename Space>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<Space>::check_bucket_id(unsigned bucketId, const char* valuesType, const char* file, int line) const
{
  if (bucketId >= static_cast<unsigned>(this->m_numBuckets)) {
    if (line == 0) {
      printf("Error: Called FieldData::%s_values() for Field '%s' with an out-of-bounds Bucket ID (%u) for a "
             "DeviceFieldMetaData array with extent %i.\n", valuesType, this->field_name(), bucketId,
             this->m_numBuckets);
    }
    else {
      printf("Error: %s:%i: Called FieldData::%s_values() for Field '%s' with an out-of-bounds Bucket ID (%u) for a "
             "DeviceFieldMetaData array with extent %i.\n", file, line, valuesType, this->field_name(), bucketId,
             this->m_numBuckets);
    }
    STK_NGP_ThrowErrorMsg("Field consistency error.");
  }
}

//------------------------------------------------------------------------------
template <typename Space>
KOKKOS_INLINE_FUNCTION void
ConstFieldDataBytes<Space>::check_bucket_ordinal(unsigned bucketId, unsigned bucketOrd, const char* file, int line) const
{
  // Only trip if we're referencing an out-of-bounds Entity in a Bucket where the Field is registered,
  // because accessing an EntityValues where the Field isn't valid should be allowed.  This allows users
  // to query EntityValues::is_field_defined() inside a loop.
  if ((bucketOrd >= static_cast<unsigned>(this->m_deviceFieldMetaData[bucketId].m_bucketSize)) &&
      (this->m_deviceFieldMetaData[bucketId].m_bucketSize > 0)) {
    if (line == 0) {
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
ConstFieldDataBytes<stk::ngp::HostSpace>::ConstFieldDataBytes()
  : FieldDataBase(),
    m_fieldMetaData(nullptr),
    m_bulk(nullptr),
    m_dataTraits(&stk::mesh::data_traits<void>()),
    m_fieldName(nullptr),
    m_fieldMetaDataModCount(nullptr),
    m_numBuckets(0),
    m_localFieldMetaDataModCount(0)
{
}

//------------------------------------------------------------------------------
inline
ConstFieldDataBytes<stk::ngp::HostSpace>::ConstFieldDataBytes(EntityRank entityRank, Ordinal fieldOrdinal,
                                                              const DataTraits& dataTraits, Layout dataLayout)
  : FieldDataBase(fieldOrdinal, entityRank, dataLayout),
    m_fieldMetaData(nullptr),
    m_bulk(nullptr),
    m_dataTraits(&dataTraits),
    m_fieldName(nullptr),
    m_fieldMetaDataModCount(nullptr),
    m_numBuckets(0),
    m_localFieldMetaDataModCount(0)
{
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Auto>(Entity entity,
                                                                     const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Auto>(const MeshIndex& mi,
                                                                     const char* file, int line) const
{
  check_mesh(mi.bucket->mesh(), "Entity", file, line);
  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  if (m_layout == Layout::Right) {
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Auto>(const FastMeshIndex& fmi,
                                                                     const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  if (m_layout == Layout::Right) {
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of);
  }
  else if (m_layout == Layout::Left) {
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data + this->m_dataTraits->alignment_of * fmi.bucket_ord,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Auto>(const Bucket& bucket,
                                                                     const char* file, int line) const
{
  check_mesh(bucket.mesh(), "Bucket", file, line);
  check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  if (m_layout == Layout::Right) {
    return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, 0);  // Keep compiler happy
  }
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>
ConstFieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Auto>(int bucketId,
                                                                     const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  if (m_layout == Layout::Right) {
    return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize);
  }
  else if (m_layout == Layout::Left) {
    return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(
          fieldMetaData.m_data,
          fieldMetaData.m_bytesPerEntity,
          this->m_dataTraits->alignment_of,
          fieldMetaData.m_bucketSize,
          fieldMetaData.m_bucketCapacity);
  }
  else {
    STK_ThrowErrorMsg("Unsupported host data layout: " << m_layout);
    return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Auto>(nullptr, 0, 0, 0, 0);  // Keep compiler happy
  }
}


//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Left>(Entity entity,
                                                                     const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Left>(const MeshIndex& mi,
                                                                     const char* file, int line) const
{
  check_mesh(mi.bucket->mesh(), "Entity", file, line);
  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Left>(const FastMeshIndex& fmi,
                                                                     const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data + this->m_dataTraits->alignment_of * fmi.bucket_ord,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Left>(const Bucket& bucket,
                                                                     const char* file, int line) const
{
  check_mesh(bucket.mesh(), "Bucket", file, line);
  check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>
ConstFieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Left>(int bucketId,
                                                                     const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Left>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize,
        fieldMetaData.m_bucketCapacity);
}


//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Right>(Entity entity,
                                                                      const char* file, int line) const
{
  const MeshIndex& mi = this->mesh().mesh_index(entity);

  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Right>(const MeshIndex& mi,
                                                                      const char* file, int line) const
{
  check_mesh(mi.bucket->mesh(), "Entity", file, line);
  check_rank(mi.bucket->entity_rank(), "Entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[mi.bucket->bucket_id()];

  return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * mi.bucket_ordinal,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostSpace>::entity_bytes<Layout::Right>(const FastMeshIndex& fmi,
                                                                      const char* file, int line) const
{
  check_bucket_id(fmi.bucket_id, "entity", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[fmi.bucket_id];

  return EntityBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data + fieldMetaData.m_bytesPerEntity * fmi.bucket_ord,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Right>(const Bucket& bucket,
                                                                      const char* file, int line) const
{
  check_mesh(bucket.mesh(), "Bucket", file, line);
  check_rank(bucket.entity_rank(), "Bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucket.bucket_id()];

  return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize);
}

//------------------------------------------------------------------------------
template <>
inline
BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>
ConstFieldDataBytes<stk::ngp::HostSpace>::bucket_bytes<Layout::Right>(int bucketId,
                                                                      const char* file, int line) const
{
  check_bucket_id(bucketId, "bucket", file, line);

  const FieldMetaData& fieldMetaData = this->m_fieldMetaData[bucketId];

  return BucketBytes<const std::byte, stk::ngp::HostSpace, Layout::Right>(
        fieldMetaData.m_data,
        fieldMetaData.m_bytesPerEntity,
        this->m_dataTraits->alignment_of,
        fieldMetaData.m_bucketSize);
}


//------------------------------------------------------------------------------
inline const char*
ConstFieldDataBytes<stk::ngp::HostSpace>::field_name() const
{
  return m_fieldName;
}

//------------------------------------------------------------------------------
inline int
ConstFieldDataBytes<stk::ngp::HostSpace>::num_buckets() const
{
  return m_numBuckets;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::set_up_to_date()
{
  m_localFieldMetaDataModCount = *m_fieldMetaDataModCount;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::set_field_name(const char* fieldName)
{
  m_fieldName = fieldName;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::set_field_meta_data_mod_count_pointer(unsigned* fieldMetaDataModCountPtr)
{
  m_fieldMetaDataModCount = fieldMetaDataModCountPtr;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::set_field_meta_data_mod_count(unsigned fieldMetaDataModCount)
{
  *m_fieldMetaDataModCount = fieldMetaDataModCount;
}

//------------------------------------------------------------------------------
inline unsigned
ConstFieldDataBytes<stk::ngp::HostSpace>::field_meta_data_mod_count()
{
  return *m_fieldMetaDataModCount;
}

//------------------------------------------------------------------------------
inline const DataTraits&
ConstFieldDataBytes<stk::ngp::HostSpace>::data_traits() const
{
  return *m_dataTraits;
}

//------------------------------------------------------------------------------
inline BulkData&
ConstFieldDataBytes<stk::ngp::HostSpace>::mesh()
{
  STK_ThrowAssert(m_bulk != nullptr);
  return *m_bulk;
}

//------------------------------------------------------------------------------
inline const BulkData&
ConstFieldDataBytes<stk::ngp::HostSpace>::mesh() const
{
  STK_ThrowAssert(m_bulk != nullptr);
  return *m_bulk;
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::set_mesh(BulkData* bulkData)
{
  m_bulk = bulkData;
  m_localFieldMetaDataModCount = 0;
}

//------------------------------------------------------------------------------
inline bool
ConstFieldDataBytes<stk::ngp::HostSpace>::needs_update() const
{
  return (m_localFieldMetaDataModCount != *m_fieldMetaDataModCount);
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::swap_field_data(FieldDataBase& other)
{
  ConstFieldDataBytes<stk::ngp::HostSpace>* otherFieldBytes =
      dynamic_cast<ConstFieldDataBytes<stk::ngp::HostSpace>*>(&other);
  STK_ThrowRequireMsg(otherFieldBytes != nullptr,
                      "ConstFieldDataBytes::swap_field_data() called with an imcompatible ConstFieldDataBytes object.");
  STK_ThrowRequireMsg(this->m_numBuckets == otherFieldBytes->m_numBuckets,
                      "Multiple states of the same Field (" << this->m_fieldName <<
                      ") must be registered on identical subsets of the mesh.");

  const unsigned numBuckets = this->m_numBuckets;

  for (unsigned bucketIdx = 0; bucketIdx < numBuckets; ++bucketIdx) {
    std::swap(this->m_fieldMetaData[bucketIdx].m_data, otherFieldBytes->m_fieldMetaData[bucketIdx].m_data);
  }

  std::swap(this->m_fieldMetaDataModCount, otherFieldBytes->m_fieldMetaDataModCount);
  std::swap(this->m_localFieldMetaDataModCount, otherFieldBytes->m_localFieldMetaDataModCount);
}

#if !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK)

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::check_updated_field(const char* file, int line) const
{
  STK_ThrowRequireMsg(m_localFieldMetaDataModCount == *m_fieldMetaDataModCount,
                      source_location_string(file, line) << "Accessing out-of-date FieldData after a host or device "
                      "mesh modification for Field '" << field_name() << "'.  Please re-acquire this FieldData "
                      "instance, potentially after calling NgpMesh::update_bulk_data() if you have modified the "
                      "mesh on device.");
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::check_mesh(const stk::mesh::BulkData& bulk, const char* target,
                                                     const char* file, int line) const
{
  STK_ThrowRequireMsg(&bulk == &mesh(),
                      source_location_string(file, line) << "Accessing " << target <<
                      " from a different mesh for Field '" << field_name() << "'.");
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::check_rank(stk::mesh::EntityRank targetRank, const char* target,
                                                     const char* file, int line) const
{
  STK_ThrowRequireMsg(entity_rank() == targetRank,
                      source_location_string(file, line) << "Accessing " << target << " with rank " << targetRank <<
                      " for Field '" << field_name() << "' with rank " << entity_rank() <<
                      ".  Are the Field and " << target << " from the same mesh?");
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::check_bucket_id(unsigned bucketId, const char* valuesType,
                                                          const char* file, int line) const
{
  STK_ThrowRequireMsg(bucketId < static_cast<unsigned>(m_numBuckets),
                      source_location_string(file, line) << "Called FieldData::" << valuesType <<
                      "_values() for Field '" << field_name() << "' with an out-of-bounds Bucket ID (" << bucketId <<
                      ") for a FieldMetaData array with size " << m_numBuckets);
}

//------------------------------------------------------------------------------
inline void
ConstFieldDataBytes<stk::ngp::HostSpace>::check_bucket_ordinal(unsigned bucketId, unsigned bucketOrd,
                                                               const char* file, int line) const
{
  // Only trip if we're referencing an out-of-bounds Entity in a Bucket where the Field is registered,
  // because accessing an EntityValues where the Field isn't valid should be allowed.  This allows users
  // to query EntityValues::is_field_defined() inside a loop.
  STK_ThrowRequireMsg((bucketOrd < static_cast<unsigned>(m_fieldMetaData[bucketId].m_bucketSize)) ||
                      (m_fieldMetaData[bucketId].m_bucketSize == 0),
                      source_location_string(file, line) << "Called FieldData::entity_values() for Field '" <<
                      field_name() << "' with an out-of-bounds Bucket ordinal (" << bucketOrd << ") for Bucket " <<
                      bucketId << " with size " << m_fieldMetaData[bucketId].m_bucketSize);
}

#endif  // !defined(NDEBUG) || defined(STK_FIELD_BOUNDS_CHECK

//==============================================================================

}

#endif // STK_CONSTFIELDDATABYTES_HPP
