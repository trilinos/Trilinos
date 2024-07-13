// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_TYPEDEFS_HPP_
#define _COMPADRE_TYPEDEFS_HPP_

#include "Compadre_Config.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <type_traits>
#include <vector>
#include <sstream>
#include <cstddef>
#include <functional>
#include <string>

/*!
 
  Data types in Compadre Toolkit:

    - Intention is to do local work, i.e. on a single node, so the default ordinal is local_index_type
    - When doing pointer arithmetic, it is possible to overflow local_index_type, so use global_index_type

*/

// Indices and data types
typedef double      scalar_type;
typedef int         local_index_type;
typedef std::size_t global_index_type;

// helper function when doing pointer arithmetic
#define TO_GLOBAL(variable) ((global_index_type)variable)

// KOKKOS TYPEDEFS

// execution spaces
typedef Kokkos::DefaultHostExecutionSpace host_execution_space;
typedef Kokkos::DefaultExecutionSpace device_execution_space;

// memory spaces
typedef typename host_execution_space::memory_space host_memory_space;
#ifdef COMPADRE_USE_CUDA
    typedef typename Kokkos::CudaSpace device_memory_space;
#else
    typedef typename device_execution_space::memory_space device_memory_space;
#endif

// team policies
typedef typename Kokkos::TeamPolicy<device_execution_space> team_policy;
typedef typename team_policy::member_type  member_type;

typedef typename Kokkos::TeamPolicy<host_execution_space> host_team_policy;
typedef typename host_team_policy::member_type  host_member_type;

// layout types
typedef Kokkos::LayoutRight layout_right;
typedef Kokkos::LayoutLeft layout_left;

// unmanaged data wrappers
typedef Kokkos::View<double**, layout_right, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            scratch_matrix_right_type;
typedef Kokkos::View<double**, layout_left, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            scratch_matrix_left_type;
typedef Kokkos::View<double*, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            scratch_vector_type;
typedef Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            scratch_local_index_type;

// host equivalents
typedef Kokkos::View<double**, layout_right, host_execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            host_scratch_matrix_right_type;
typedef Kokkos::View<double**, layout_left, host_execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            host_scratch_matrix_left_type;
typedef Kokkos::View<double*, host_execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            host_scratch_vector_type;
typedef Kokkos::View<int*, host_execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
            host_scratch_local_index_type;

// managed device data views
typedef Kokkos::View<double**, layout_right, device_memory_space>
            device_managed_matrix_right_type;
typedef Kokkos::View<double**, layout_left, device_memory_space>
            device_managed_matrix_left_type;
typedef Kokkos::View<double*, device_memory_space>
            device_managed_vector_type;
typedef Kokkos::View<int*, device_memory_space>
            device_managed_local_index_type;

// managed host data views
typedef Kokkos::View<double**, layout_right, host_execution_space>
            host_managed_matrix_right_type;
typedef Kokkos::View<double**, layout_left, host_execution_space>
            host_managed_matrix_left_type;
typedef Kokkos::View<double*, host_execution_space>
            host_managed_vector_type;
typedef Kokkos::View<int*, host_execution_space>
            host_managed_local_index_type;

// random number generator
typedef Kokkos::Random_XorShift64_Pool<> pool_type;
typedef typename pool_type::generator_type generator_type;

// KOKKOS_VERSION % 100 is the patch level
// KOKKOS_VERSION / 100 % 100 is the minor version
// KOKKOS_VERSION / 10000 is the major version
#ifdef KOKKOS_VERSION
  #define COMPADRE_KOKKOS_VERSION_MAJOR KOKKOS_VERSION / 10000
  #define COMPADRE_KOKKOS_VERSION_MINOR KOKKOS_VERSION / 100 % 100 
  #if COMPADRE_KOKKOS_VERSION_MAJOR < 4
    #if COMPADRE_KOKKOS_VERSION_MINOR >= 7
        using KokkosInitArguments = Kokkos::InitializationSettings;
        #define COMPADRE_KOKKOS_GREATEREQUAL_3_7
        constexpr char KOKKOS_THREADS_ARG[] = "--kokkos-num-threads";
    #elif COMPADRE_KOKKOS_VERSION_MINOR < 7
        using KokkosInitArguments = Kokkos::InitArguments;
        constexpr char KOKKOS_THREADS_ARG[] = "--kokkos-threads";
    #endif
  #elif COMPADRE_KOKKOS_VERSION_MAJOR >= 4
      using KokkosInitArguments = Kokkos::InitializationSettings;
      #define COMPADRE_KOKKOS_GREATEREQUAL_3_7
      constexpr char KOKKOS_THREADS_ARG[] = "--kokkos-num-threads";
  #endif
#else // older version
  using KokkosInitArguments = Kokkos::InitArguments;
  constexpr char KOKKOS_THREADS_ARG[] = "--kokkos-threads";
#endif

template< bool B, class T = void >
using enable_if_t = typename std::enable_if<B,T>::type;

template<typename T>
typename std::enable_if<1==T::rank,T>::type createView(std::string str, int dim_0, int dim_1)
{ return T(str, dim_0); }

template<typename T>
typename std::enable_if<2==T::rank,T>::type createView(std::string str, int dim_0, int dim_1)
{ return T(str, dim_0, dim_1); }

//void compadre_rethrow_exception(std::exception &e, const std::string &extra_message) {
//    std::cout << extra_message + "\n\n" + e.what() << std::endl;
//}

//! compadre_assert_release is used for assertions that should always be checked, but generally 
//! are not expensive to verify or are not called frequently. 
# define compadre_assert_release(condition) do {                                \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

//! compadre_kernel_assert_release is similar to compadre_assert_release, but is a call on the device, 
//! namely inside of a function marked KOKKOS_INLINE_FUNCTION
# define compadre_kernel_assert_release(condition) do { \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)

//! compadre_assert_debug is used for assertions that are checked in loops, as these significantly
//! impact performance. When NDEBUG is set, these conditions are not checked.
#ifdef COMPADRE_DEBUG
# define compadre_assert_debug(condition) do {                                \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
# define compadre_kernel_assert_debug(condition) do { \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#else
#  define compadre_assert_debug(condition)
#  define compadre_kernel_assert_debug(condition)
#endif
//! compadre_kernel_assert_debug is similar to compadre_assert_debug, but is a call on the device, 
//! namely inside of a function marked KOKKOS_INLINE_FUNCTION

#ifdef COMPADRE_EXTREME_DEBUG
# define compadre_assert_extreme_debug(condition) do {                                \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
# define compadre_kernel_assert_extreme_debug(condition) do { \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#else
#  define compadre_assert_extreme_debug(condition)
#  define compadre_kernel_assert_extreme_debug(condition)
#endif
//! compadre_kernel_assert_extreme_debug is similar to compadre_assert_debug, but is a call on the device, 
//! namely inside of a function marked KOKKOS_INLINE_FUNCTION

#endif
