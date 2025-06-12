#ifndef NGPFIELDBASE_HPP
#define NGPFIELDBASE_HPP

#include "stk_util/stk_config.h"
#include "stk_mesh/base/Selector.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
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
  virtual void modify_on_host() = 0;
  virtual void modify_on_host(const Selector& selector) = 0;
  virtual void modify_on_device() = 0;
  virtual void modify_on_device(const Selector& selector) = 0;
  virtual bool need_sync_to_host() const = 0;
  virtual bool need_sync_to_device() const = 0;
  virtual void clear_sync_state() = 0;
  virtual void clear_host_sync_state() = 0;
  virtual void clear_device_sync_state() = 0;
  virtual void sync_to_host() = 0;
  virtual void sync_to_host(const stk::ngp::ExecSpace& newExecSpace) = 0;
  virtual void sync_to_host(stk::ngp::ExecSpace&& newExecSpace) = 0;
  virtual void sync_to_device() = 0;
  virtual void sync_to_device(const stk::ngp::ExecSpace& newExecSpace) = 0;
  virtual void sync_to_device(stk::ngp::ExecSpace&& newExecSpace) = 0;
  virtual size_t synchronized_count() const = 0;
  virtual size_t num_syncs_to_host() const = 0;
  virtual size_t num_syncs_to_device() const = 0;
  virtual void fence() = 0;

  virtual void update_field(const stk::ngp::ExecSpace& newExecSpace) = 0;
  virtual void update_field(stk::ngp::ExecSpace&& newExecSpace) = 0;
};

}
}

#endif // NGPFIELDBASE_HPP
