#ifndef NGPMESHBASE_HPP
#define NGPMESHBASE_HPP

#include "stk_util/stk_config.h"
#include <stddef.h>

namespace stk {
namespace mesh {

class NgpMeshBase
{
public:
  KOKKOS_DEFAULTED_FUNCTION NgpMeshBase() = default;
  KOKKOS_FUNCTION virtual ~NgpMeshBase() {}

  KOKKOS_DEFAULTED_FUNCTION NgpMeshBase(const NgpMeshBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION NgpMeshBase(NgpMeshBase &&) = default;
  KOKKOS_DEFAULTED_FUNCTION NgpMeshBase& operator=(const NgpMeshBase&) = default;
  KOKKOS_DEFAULTED_FUNCTION NgpMeshBase& operator=(NgpMeshBase&&) = default;

  virtual void update_mesh() = 0;
  virtual void update_bulk_data() = 0;
#ifndef STK_HIDE_DEPRECATED_CODE
  STK_DEPRECATED_MSG("Use need_update_bulk_data() instead.")
  virtual bool need_sync_to_host() const = 0;
#endif
  virtual bool need_update_bulk_data() const = 0;
  virtual unsigned synchronized_count() const = 0;
};

}}

#endif // NGPMESHBASE_HPP
