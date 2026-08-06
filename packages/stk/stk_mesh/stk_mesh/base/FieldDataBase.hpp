#ifndef STK_FIELDDATABASE_HPP
#define STK_FIELDDATABASE_HPP

#include "stk_mesh/base/Types.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace sierra::Fmwk {
  class Region;
}

namespace stk::mesh {

template <typename MemSpace_> class DeviceFieldDataManager;
template <typename T, typename MemSpace_> class HostField;
template <typename T, typename MemSpace_> class DeviceField;
template <typename T, typename MemSpace_, Layout DataLayout_> class ConstFieldData;

namespace impl {

class MeshModification;
template <typename MemSpace_> class DeviceBucketRepository;

KOKKOS_FORCEINLINE_FUNCTION
constexpr bool is_called_on_host() {
  KOKKOS_IF_ON_HOST(return true;);
  KOKKOS_IF_ON_DEVICE(return false;);
}

}

struct FieldDataCopyTracking {
  const char* m_file = nullptr;
  int m_line = 0;
  int m_copyCount = 0;
};

class FieldDataBase
{
public:
  KOKKOS_FUNCTION FieldDataBase()
    : m_copyTracking(nullptr),
      m_ordinal(InvalidOrdinal),
      m_rank(InvalidEntityRank),
      m_layout(Layout::Auto),
      m_accessTag(InvalidAccess)
  {}

  FieldDataBase(Ordinal ordinal, EntityRank rank, Layout layout)
    : m_copyTracking(nullptr),
      m_ordinal(ordinal),
      m_rank(rank),
      m_layout(layout),
      m_accessTag(InvalidAccess)
  {}

  FieldDataBase(FieldDataCopyTracking* copyTracking, Ordinal ordinal, EntityRank rank, Layout layout)
    : m_copyTracking(copyTracking),
      m_ordinal(ordinal),
      m_rank(rank),
      m_layout(layout),
      m_accessTag(InvalidAccess)
  {}

  KOKKOS_FUNCTION virtual ~FieldDataBase() {}

  KOKKOS_DEFAULTED_FUNCTION FieldDataBase(const FieldDataBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBase(FieldDataBase&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBase& operator=(const FieldDataBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBase& operator=(FieldDataBase&&) = default;

  KOKKOS_INLINE_FUNCTION Layout data_layout() const {
    return m_layout;
  }

  KOKKOS_INLINE_FUNCTION EntityRank entity_rank() const {
    return m_rank;
  }

  KOKKOS_INLINE_FUNCTION Ordinal field_ordinal() const {
    return m_ordinal;
  }

  KOKKOS_INLINE_FUNCTION FieldAccessTag access_tag() const {
    return m_accessTag;
  }

  virtual bool needs_update() const = 0;

protected:
  friend FieldBase;
  friend BulkData;
  friend impl::MeshModification;
  friend sierra::Fmwk::Region;

  void track_copy(FieldAccessTag tag) {
    if (tag < NumTrackedFieldAccessTags) {
      ++m_copyTracking[tag].m_copyCount;
    }
    m_accessTag = tag;
  }

  void track_copy(FieldAccessTag tag, [[maybe_unused]] const char* file, [[maybe_unused]] int line) {
    if (tag < NumTrackedFieldAccessTags) {
      ++m_copyTracking[tag].m_copyCount;
      m_copyTracking[tag].m_file = file;
      m_copyTracking[tag].m_line = line;
    }
    m_accessTag = tag;
  }

  void release_copy() {
    if (m_accessTag < NumTrackedFieldAccessTags) {
      --m_copyTracking[m_accessTag].m_copyCount;
    }
  }

  bool has_copies(FieldAccessTag tag) const {
    return (m_copyTracking[tag].m_copyCount > 0);
  }

  const FieldDataCopyTracking& copy_tracking(FieldAccessTag tag) const {
    return m_copyTracking[tag];
  }

  void set_copy_tracking(FieldDataCopyTracking* copyTracking) {
    m_copyTracking = copyTracking;
  }

  virtual void set_mesh(BulkData* bulkData) = 0;

  virtual void swap_field_data(FieldDataBase& other) = 0;
  virtual void update_host_bucket_pointers() = 0;
  virtual void incomplete_swap_field_data(FieldDataBase& other) = 0;

  virtual void sync_to_host(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout) = 0;
  virtual void sync_to_device(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout) = 0;
  virtual void update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout, bool needsSync) = 0;
  virtual void fence(const stk::ngp::ExecSpace& execSpace) = 0;

  FieldDataCopyTracking* m_copyTracking;  //  8 : pointer
  Ordinal m_ordinal;                      //  4 : unsigned
  EntityRank m_rank;                      //  1 : int8_t
  Layout m_layout;                        //  1 : uint8_t
  FieldAccessTag m_accessTag;             //  1 : uint8_t
};

}

#endif // STK_FIELDDATABASE_HPP
