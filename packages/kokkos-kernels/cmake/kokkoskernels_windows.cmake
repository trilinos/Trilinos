# Detect MSVC ABI once; result persists for the second include (below /bigobj).
# cl.exe, clang-cl, Intel icl/icx on Windows all use the MSVC C++ ABI.
# MinGW and Cygwin use the GCC/POSIX model – excluded intentionally.
if(WIN32 AND NOT CYGWIN AND NOT DEFINED KK_WINDOWS_MSVC_ABI_CXX)
  set(KK_WINDOWS_MSVC_ABI_CXX FALSE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(KK_WINDOWS_MSVC_ABI_CXX TRUE)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
         "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
    set(KK_WINDOWS_MSVC_ABI_CXX TRUE)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(KK_WINDOWS_MSVC_ABI_CXX TRUE)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(KK_WINDOWS_MSVC_ABI_CXX TRUE)
  endif()
endif()

# x64 -------------------------------------------------------------------------
# 16-byte host atomics need CMPXCHG16B / _InterlockedCompareExchange128 (x64 only).
# clang-cl: -mcx16 below.  MSVC: Compare_Exchange_MSVC.hpp (no -mcx16).
# Do not key off -A x64 alone: Ninja + Hostx64/x64/cl.exe is valid without it.
if(KK_WINDOWS_MSVC_ABI_CXX AND NOT DEFINED KK_WINDOWS_X64_CHECK_DONE)
  if(DEFINED CMAKE_SIZEOF_VOID_P AND NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    if(CMAKE_GENERATOR MATCHES "Visual Studio")
      message(FATAL_ERROR
        "KokkosKernels requires a 64-bit Windows (x64) build. "
        "Reconfigure with -A x64, for example:\n"
        "  cmake -G \"Visual Studio 17 2022\" -A x64 -S <src> -B <build>")
    else()
      message(FATAL_ERROR
        "KokkosKernels requires a 64-bit Windows build but detected "
        "${CMAKE_SIZEOF_VOID_P}-byte pointers. For Ninja or single-config generators, "
        "point CMAKE_CXX_COMPILER (and CUDA host compiler if used) to the x64 toolset, "
        "e.g. .../VC/Tools/MSVC/<ver>/bin/Hostx64/x64/cl.exe")
    endif()
  endif()
  set(KK_WINDOWS_X64_CHECK_DONE TRUE)
endif()

# /EHsc -----------------------------------------------------------------------
# Required for stack-unwinding semantics in throw/catch on the MSVC ABI.
# Without it the compiler generates non-unwind tables and <chrono>/<stdexcept>
# emit C4530; behaviour on actual throws is undefined.
# https://learn.microsoft.com/en-us/cpp/build/reference/eh-exception-handling-model
# https://learn.microsoft.com/en-us/cpp/error-messages/compiler-warnings/compiler-warning-level-3-c4530
if(KK_WINDOWS_MSVC_ABI_CXX AND NOT DEFINED KK_WINDOWS_EHSC_DONE)
  if(CMAKE_CXX_FLAGS MATCHES "(^| )/EHa( |$)")
    message(WARNING
      "CMAKE_CXX_FLAGS contains /EHa; KokkosKernels will not add /EHsc "
      "(exception model left to the user).")
  elseif(CMAKE_CXX_FLAGS MATCHES "(^| )/EHsc-( |$)")
    message(WARNING
      "CMAKE_CXX_FLAGS contains /EHsc-; KokkosKernels will not add /EHsc "
      "(C++ exceptions explicitly disabled by the user).")
  elseif(NOT CMAKE_CXX_FLAGS MATCHES "(^| )/EHsc( |$)")
    message(STATUS
      "Adding '/EHsc' for C++ exception handling on Windows "
      "(${CMAKE_CXX_COMPILER_ID})")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/EHsc>)
  else()
    # /EHsc already in CMAKE_CXX_FLAGS (e.g. MSVC generator default).
    # Pin it via add_compile_options so a later -DCMAKE_CXX_FLAGS= override
    # cannot silently drop unwind semantics for targets in this tree.
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/EHsc>)
  endif()
  set(KK_WINDOWS_EHSC_DONE TRUE)
endif()

# /bigobj (CUDA only) ---------------------------------------------------------
# C1128 on the host pass of nvcc (.cudafe1.cpp): template-heavy ETI TUs exceed
# the default COFF section limit. With COMPILE_AS_CMAKE_LANGUAGE=CUDA, library
# sources compile as CUDA; cl.exe needs /bigobj via -Xcompiler=/bigobj.
# https://learn.microsoft.com/en-us/cpp/error-messages/compiler-errors-1/fatal-error-c1128
# Included a second time (after kokkos_backends.cmake) so KOKKOS_ENABLE_CUDA is set.
if(KK_WINDOWS_MSVC_ABI_CXX AND KOKKOS_ENABLE_CUDA AND
   NOT DEFINED KK_WINDOWS_BIGOBJ_DONE)
  message(STATUS
    "Adding '-Xcompiler=/bigobj' for CUDA language on Windows "
    "(${CMAKE_CXX_COMPILER_ID})")
  add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/bigobj>)
  set(KK_WINDOWS_BIGOBJ_DONE TRUE)
endif()

# -mcx16 -----------------------------------------------------------------------
# LNK2001 __atomic_compare_exchange_16: clang/clang-cl on Windows uses GCC-style
# __atomic_* builtins (desul atomics). Without -mcx16, 16-byte compare-exchange
# (e.g. Kokkos::complex<double>, Dummy16ByteValue reductions) lowers to libatomic
# which is not linked in the MSVC link environment.
# MSVC (cl.exe) uses _InterlockedCompareExchange128 instead; -mcx16 is clang-only.
# Kokkos adds the same flag in cmake/kokkos_arch.cmake for Clang+WIN32.
if(KK_WINDOWS_MSVC_ABI_CXX AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
   NOT DEFINED KK_WINDOWS_MCX16_DONE)
  message(STATUS
    "Adding '-mcx16' for 16-byte atomics on Windows "
    "(${CMAKE_CXX_COMPILER_ID})")
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mcx16>)
  set(KK_WINDOWS_MCX16_DONE TRUE)
endif()
