#ifndef NGPFIELDBASE_HPP
#define NGPFIELDBASE_HPP

#include "stk_util/stk_config.h"
#include <stddef.h>

namespace stk {
namespace mesh {

class NgpFieldBase
{
public:
  KOKKOS_DEFAULTED_FUNCTION NgpFieldBase() = default;
  KOKKOS_FUNCTION virtual ~NgpFieldBase() {}
  virtual void update_field(bool needToSyncAllDataToDevice = false) = 0;
  virtual void rotate_multistate_data() = 0;
  virtual void modify_on_host() = 0;
  virtual void modify_on_device() = 0;
  virtual void clear_sync_state() = 0;
  virtual void sync_to_host() = 0;
  virtual void sync_to_device() = 0;
  virtual size_t synchronized_count() const = 0;
  virtual size_t num_syncs_to_host() const = 0;
  virtual size_t num_syncs_to_device() const = 0;
#ifdef STK_DEBUG_FIELD_SYNC
  virtual void detect_device_field_modification() = 0;
  virtual void update_debug_storage(size_t hostSynchronizedCount) = 0;
  virtual bool any_device_field_modification() const = 0;
#endif
};

#ifdef KOKKOS_ENABLE_CUDA
#define ORDER_INDICES(i,j) j,i
#else
#define ORDER_INDICES(i,j) i,j
#endif

}
}

#endif // NGPFIELDBASE_HPP
