#ifndef FIELDSYNCDEBUGGING_HPP
#define FIELDSYNCDEBUGGING_HPP

#ifdef STK_DEBUG_FIELD_SYNC

#include "Kokkos_Core.hpp"
#include "stk_mesh/base/NgpSpaces.hpp"

#ifndef __has_builtin
  #define __has_builtin(x) 0  // Compatibility with non-clang compilers.
#endif

#if !defined(__INTEL_COMPILER) && !(defined(__clang__) && !(__has_builtin(__builtin_LINE)))
#define HOST_USE_LOCATION_BUILTINS
#define HOST_DEBUG_FILE_NAME __builtin_FILE()
#define HOST_DEBUG_LINE_NUMBER __builtin_LINE()
#else
#define HOST_DEBUG_FILE_NAME ""
#define HOST_DEBUG_LINE_NUMBER -1
#endif

#if !defined(__INTEL_COMPILER) && !(defined(__clang__) && !(__has_builtin(__builtin_LINE))) && !defined(KOKKOS_ENABLE_CUDA)
#define DEVICE_USE_LOCATION_BUILTINS
#define DEVICE_DEBUG_FILE_NAME __builtin_FILE()
#define DEVICE_DEBUG_LINE_NUMBER __builtin_LINE()
#else
#define DEVICE_DEBUG_FILE_NAME ""
#define DEVICE_DEBUG_LINE_NUMBER -1
#endif

namespace stk {
namespace mesh {

enum LastModLocation : uint8_t {
  NONE   = 0x00,
  HOST   = 0x01,
  DEVICE = 0x02
};

using LastFieldModLocationType = Kokkos::View<LastModLocation***, Kokkos::LayoutRight, UVMMemSpace>;

struct DummyOverload
{
  DummyOverload() = default;
  ~DummyOverload() = default;
};

}
}

#endif

#endif // FIELDSYNCDEBUGGING_HPP
