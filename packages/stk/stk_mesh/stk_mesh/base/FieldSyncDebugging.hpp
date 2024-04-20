#ifndef FIELDSYNCDEBUGGING_HPP
#define FIELDSYNCDEBUGGING_HPP

#include "Kokkos_Core.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

#ifndef __has_builtin
  #define __has_builtin(x) 0  // Compatibility with non-clang compilers.
#endif

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
  #define HAS_GCC_LOCATION_BUILTINS
#endif

#if defined(__clang__) && __has_builtin(__builtin_LINE)
  #define HAS_CLANG_LOCATION_BUILTINS
#endif

#if defined(HAS_GCC_LOCATION_BUILTINS) || defined(HAS_CLANG_LOCATION_BUILTINS)
  #define HAS_LOCATION_BUILTINS
#endif

#if defined(HAS_LOCATION_BUILTINS)
  #define HOST_USE_LOCATION_BUILTINS
  #define HOST_DEBUG_FILE_NAME __builtin_FILE()
  #define HOST_DEBUG_LINE_NUMBER __builtin_LINE()
#else
  #define HOST_DEBUG_FILE_NAME ""
  #define HOST_DEBUG_LINE_NUMBER -1
#endif

#if defined(HAS_LOCATION_BUILTINS) && !defined(STK_ENABLE_GPU)
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
  NONE           = 0x00,
  HOST           = 0x01,
  DEVICE         = 0x02,
  HOST_OR_DEVICE = HOST | DEVICE
};

using LastFieldModLocationType = Kokkos::View<LastModLocation***, Kokkos::LayoutRight, stk::ngp::UVMMemSpace>;

template <typename T>
using ScalarUvmType = Kokkos::View<T, stk::ngp::UVMMemSpace>;

struct DummyOverload
{
  DummyOverload() = default;
  ~DummyOverload() = default;
};

}
}

#endif // FIELDSYNCDEBUGGING_HPP
