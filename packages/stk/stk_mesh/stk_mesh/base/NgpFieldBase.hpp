#ifndef NGPFIELDBASE_HPP
#define NGPFIELDBASE_HPP

#include "stk_util/stk_config.h"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/NgpSpaces.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include <stddef.h>

namespace stk {
namespace mesh {

class NgpFieldBase
{
public:
  KOKKOS_DEFAULTED_FUNCTION NgpFieldBase() = default;
  KOKKOS_DEFAULTED_FUNCTION NgpFieldBase(const NgpFieldBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION NgpFieldBase(NgpFieldBase&&) = default;
  KOKKOS_FUNCTION NgpFieldBase& operator=(const NgpFieldBase&) { return *this; }
  KOKKOS_FUNCTION NgpFieldBase& operator=(NgpFieldBase&&) { return *this; }
  KOKKOS_FUNCTION virtual ~NgpFieldBase() {}
  virtual void update_field(bool needToSyncAllDataToDevice = false) = 0;
  virtual void rotate_multistate_data() = 0;
  virtual void modify_on_host() = 0;
  virtual void modify_on_host(const Selector& selector) = 0;
  virtual void modify_on_device() = 0;
  virtual void modify_on_device(const Selector& selector) = 0;
  virtual void clear_sync_state() = 0;
  virtual void clear_host_sync_state() = 0;
  virtual void clear_device_sync_state() = 0;
  virtual void sync_to_host() = 0;
  virtual void sync_to_device() = 0;
  virtual size_t synchronized_count() const = 0;
  virtual size_t num_syncs_to_host() const = 0;
  virtual size_t num_syncs_to_device() const = 0;

  virtual void debug_modification_begin() = 0;
  virtual void debug_modification_end(size_t synchronizationCount) = 0;
  virtual void debug_detect_device_field_modification() = 0;
  virtual unsigned debug_get_bucket_offset(unsigned bucketOrdinal) const = 0;
};

}
}

#endif // NGPFIELDBASE_HPP
