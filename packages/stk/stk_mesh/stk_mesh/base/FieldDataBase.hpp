#ifndef STK_FIELDDATABASE_HPP
#define STK_FIELDDATABASE_HPP

#include "stk_mesh/base/Types.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace sierra::Fmwk {
  class Region;
}

namespace stk::mesh {

template <typename _MemSpace> class DeviceFieldDataManager;
template <typename _T, typename _MemSpace> class HostField;
template <typename _T, typename _MemSpace> class DeviceField;
template <typename _T, typename _MemSpace, Layout _DataLayout> class ConstFieldData;

namespace impl {

class MeshModification;

KOKKOS_FORCEINLINE_FUNCTION
constexpr bool is_called_on_host() {
  KOKKOS_IF_ON_HOST(return true;);
  KOKKOS_IF_ON_DEVICE(return false;);
}

}

class FieldDataBase
{
public:
  KOKKOS_FUNCTION FieldDataBase()
    : m_accessTag(InvalidAccess)
  {}

  explicit FieldDataBase(bool)
    : m_copyCount("CopyCount"),
      m_accessTag(InvalidAccess)
  {}

  KOKKOS_FUNCTION virtual ~FieldDataBase() {}

  KOKKOS_DEFAULTED_FUNCTION FieldDataBase(const FieldDataBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBase(FieldDataBase&&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBase& operator=(const FieldDataBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION FieldDataBase& operator=(FieldDataBase&&) = default;

protected:
  friend FieldBase;
  friend BulkData;
  friend impl::MeshModification;
  friend sierra::Fmwk::Region;

  void track_copy(FieldAccessTag tag) {
    if (tag < NumTrackedFieldAccessTags) {
      ++m_copyCount[tag];
    }
    m_accessTag = tag;
  }

  void release_copy() {
    if (m_accessTag < NumTrackedFieldAccessTags) {
      --m_copyCount[m_accessTag];
    }
  }

  bool has_copies(FieldAccessTag tag) const {
    return (m_copyCount[tag] > 0);
  }

  FieldAccessTag access_tag() const {
    return m_accessTag;
  }

  virtual void set_mesh(BulkData* bulkData) = 0;

  virtual bool needs_update() const = 0;
  virtual int field_data_synchronized_count() const = 0;
  virtual void swap_field_data(FieldDataBase& other) = 0;
  virtual void update_host_bucket_pointers() = 0;
  virtual void incomplete_swap_field_data(FieldDataBase& other) = 0;

  virtual bool need_device_metadata_update() = 0;
  virtual void update_device_field_metadata() = 0;

  virtual void sync_to_host(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout) = 0;
  virtual void sync_to_device(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout) = 0;
  virtual void update(const stk::ngp::ExecSpace& execSpace, Layout hostDataLayout, bool needsSync) = 0;
  virtual void fence(const stk::ngp::ExecSpace& execSpace) = 0;

  Kokkos::View<int[NumTrackedFieldAccessTags], stk::ngp::HostMemSpace> m_copyCount;
  FieldAccessTag m_accessTag;
};

}

#endif // STK_FIELDDATABASE_HPP
