// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_NUMERIC_TOOLS_FACTORY_HPP__
#define __TACHO_NUMERIC_TOOLS_FACTORY_HPP__

/// \file Tacho_NumericTools_Serial.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_NumericTools_Base.hpp"
#include "Tacho_NumericTools_LevelSet.hpp"
#include "Tacho_NumericTools_Serial.hpp"

namespace Tacho {

///
///
///
template <typename ValueType, typename DeviceType> class NumericToolsFactory;

/// partial specialization
#define TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING                                                                         \
  using ordinal_type_array = typename numeric_tools_base_type::ordinal_type_array;                                     \
  using size_type_array = typename numeric_tools_base_type::size_type_array;                                           \
  using ordinal_type_array_host = typename numeric_tools_base_type::ordinal_type_array_host;                           \
  using size_type_array_host = typename numeric_tools_base_type::size_type_array_host

#define TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER                                                                        \
  ordinal_type _method;                                                                                                \
  ordinal_type _m;                                                                                                     \
  size_type_array _ap;                                                                                                 \
  ordinal_type_array _aj;                                                                                              \
  ordinal_type_array _perm;                                                                                            \
  ordinal_type_array _peri;                                                                                            \
  ordinal_type _nsupernodes;                                                                                           \
  ordinal_type_array _supernodes;                                                                                      \
  size_type_array _gid_ptr;                                                                                            \
  ordinal_type_array _gid_colidx;                                                                                      \
  size_type_array _sid_ptr;                                                                                            \
  ordinal_type_array _sid_colidx;                                                                                      \
  ordinal_type_array _blk_colidx;                                                                                      \
  ordinal_type_array _stree_parent;                                                                                    \
  size_type_array _stree_ptr;                                                                                          \
  ordinal_type_array _stree_children;                                                                                  \
  ordinal_type_array_host _stree_level;                                                                                \
  ordinal_type_array_host _stree_roots;                                                                                \
  ordinal_type _verbose

#define TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER                                                                    \
  do {                                                                                                                 \
    _method = method;                                                                                                  \
    _m = m;                                                                                                            \
    _ap = ap;                                                                                                          \
    _aj = aj;                                                                                                          \
    _perm = perm;                                                                                                      \
    _peri = peri;                                                                                                      \
    _nsupernodes = nsupernodes;                                                                                        \
    _supernodes = supernodes;                                                                                          \
    _gid_ptr = gid_ptr;                                                                                                \
    _gid_colidx = gid_colidx;                                                                                          \
    _sid_ptr = sid_ptr;                                                                                                \
    _sid_colidx = sid_colidx;                                                                                          \
    _blk_colidx = blk_colidx;                                                                                          \
    _stree_parent = stree_parent;                                                                                      \
    _stree_ptr = stree_ptr;                                                                                            \
    _stree_children = stree_children;                                                                                  \
    _stree_level = stree_level;                                                                                        \
    _stree_roots = stree_roots;                                                                                        \
    _verbose = verbose;                                                                                                \
  } while (false)

#define TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER                                                                    \
  ordinal_type _variant;                                                                                               \
  ordinal_type _device_level_cut;                                                                                      \
  ordinal_type _device_factor_thres;                                                                                   \
  ordinal_type _device_solve_thres;                                                                                    \
  ordinal_type _nstreams

#define TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER                                                                \
  do {                                                                                                                 \
    _variant = variant;                                                                                                \
    _device_level_cut = device_level_cut;                                                                              \
    _device_factor_thres = device_factor_thres;                                                                        \
    _device_solve_thres = device_solve_thres;                                                                          \
    _nstreams = nstreams;                                                                                              \
  } while (false)

#define TACHO_NUMERIC_TOOLS_FACTORY_SERIAL_BODY                                                                        \
  do {                                                                                                                 \
    if (object == nullptr)                                                                                             \
      object = (numeric_tools_base_type *)::operator new(sizeof(numeric_tools_serial_type));                           \
                                                                                                                       \
    new (object) numeric_tools_serial_type(_method, _m, _ap, _aj, _perm, _peri, _nsupernodes, _supernodes, _gid_ptr,   \
                                           _gid_colidx, _sid_ptr, _sid_colidx, _blk_colidx, _stree_parent, _stree_ptr, \
                                           _stree_children, _stree_level, _stree_roots);                               \
  } while (false)

#define TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_name)                                         \
  do {                                                                                                                 \
    if (object == nullptr)                                                                                             \
      object = (numeric_tools_base_type *)::operator new(sizeof(numeric_tools_levelset_name));                         \
    new (object) numeric_tools_levelset_name(_method, _m, _ap, _aj, _perm, _peri, _nsupernodes, _supernodes, _gid_ptr, \
                                             _gid_colidx, _sid_ptr, _sid_colidx, _blk_colidx, _stree_parent,           \
                                             _stree_ptr, _stree_children, _stree_level, _stree_roots);                 \
    numeric_tools_levelset_name *N = dynamic_cast<numeric_tools_levelset_name *>(object);                              \
    N->initialize(_device_level_cut, _device_factor_thres, _device_solve_thres, _verbose);                             \
    N->createStream(_nstreams, _verbose);                                                                              \
  } while (false)

///
/// Serial construction
///
#if defined(KOKKOS_ENABLE_SERIAL)
template <typename ValueType> class NumericToolsFactory<ValueType, typename UseThisDevice<Kokkos::Serial>::type> {
public:
  using value_type = ValueType;
  using device_type = typename UseThisDevice<Kokkos::Serial>::type;
  using numeric_tools_base_type = NumericToolsBase<value_type, device_type>;
  using numeric_tools_serial_type = NumericToolsSerial<value_type, device_type>;
  using numeric_tools_levelset_var0_type = NumericToolsLevelSet<value_type, device_type, 0>;
  using numeric_tools_levelset_var1_type = NumericToolsLevelSet<value_type, device_type, 1>;
  using numeric_tools_levelset_var2_type = NumericToolsLevelSet<value_type, device_type, 2>;
  using numeric_tools_levelset_var3_type = NumericToolsLevelSet<value_type, device_type, 3>;

  TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING;
  TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER;
  // TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER;

  void setBaseMember(const ordinal_type method,
                     // input matrix A
                     const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                     // input permutation
                     const ordinal_type_array &perm, const ordinal_type_array &peri,
                     // supernodes
                     const ordinal_type nsupernodes, const ordinal_type_array &supernodes,
                     const size_type_array &gid_ptr, const ordinal_type_array &gid_colidx,
                     const size_type_array &sid_ptr, const ordinal_type_array &sid_colidx,
                     const ordinal_type_array &blk_colidx, const ordinal_type_array &stree_parent,
                     const size_type_array &stree_ptr, const ordinal_type_array &stree_children,
                     const ordinal_type_array_host &stree_level, const ordinal_type_array_host &stree_roots,
                     const ordinal_type verbose) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER;
  }

  void setLevelSetMember(const ordinal_type variant, const ordinal_type device_level_cut,
                         const ordinal_type device_factor_thres, const ordinal_type device_solve_thres,
                         const ordinal_type nstreams) {
    // TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER;
  }

  void createObject(numeric_tools_base_type *&object) {
    KOKKOS_IF_ON_HOST((TACHO_NUMERIC_TOOLS_FACTORY_SERIAL_BODY;))
  }
};
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
template <typename ValueType> class NumericToolsFactory<ValueType, typename UseThisDevice<Kokkos::OpenMP>::type> {
public:
  using value_type = ValueType;
  using device_type = typename UseThisDevice<Kokkos::OpenMP>::type;
  using numeric_tools_base_type = NumericToolsBase<value_type, device_type>;
  using numeric_tools_serial_type = NumericToolsSerial<value_type, device_type>;
  using numeric_tools_levelset_var0_type = NumericToolsLevelSet<value_type, device_type, 0>;
  using numeric_tools_levelset_var1_type = NumericToolsLevelSet<value_type, device_type, 1>;
  using numeric_tools_levelset_var2_type = NumericToolsLevelSet<value_type, device_type, 2>;
  using numeric_tools_levelset_var3_type = NumericToolsLevelSet<value_type, device_type, 3>;

  TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING;
  TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER;
  #define TACHO_LEVELSET_ON_HOST
  #if defined TACHO_LEVELSET_ON_HOST
  TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER;
  #endif

  void setBaseMember(const ordinal_type method,
                     // input matrix A
                     const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                     // input permutation
                     const ordinal_type_array &perm, const ordinal_type_array &peri,
                     // supernodes
                     const ordinal_type nsupernodes, const ordinal_type_array &supernodes,
                     const size_type_array &gid_ptr, const ordinal_type_array &gid_colidx,
                     const size_type_array &sid_ptr, const ordinal_type_array &sid_colidx,
                     const ordinal_type_array &blk_colidx, const ordinal_type_array &stree_parent,
                     const size_type_array &stree_ptr, const ordinal_type_array &stree_children,
                     const ordinal_type_array_host &stree_level, const ordinal_type_array_host &stree_roots,
                     const ordinal_type verbose) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER;
  }

  void setLevelSetMember(const ordinal_type variant, const ordinal_type device_level_cut,
                         const ordinal_type device_factor_thres, const ordinal_type device_solve_thres,
                         const ordinal_type nstreams) {
    #if defined TACHO_LEVELSET_ON_HOST
    TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER;
    #endif
  }

  void createObject(numeric_tools_base_type *&object) {
#if defined TACHO_LEVELSET_ON_HOST
    KOKKOS_IF_ON_HOST((
    switch (_variant) {
    case 0: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var0_type);
      break;
    }
    case 1: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var1_type);
      break;
    }
    case 2: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var2_type);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid variant input");
    }
    }))
#else
    KOKKOS_IF_ON_HOST((TACHO_NUMERIC_TOOLS_FACTORY_SERIAL_BODY;))
#endif
  }
};
#endif

#if defined(KOKKOS_ENABLE_CUDA)
template <typename ValueType> class NumericToolsFactory<ValueType, typename UseThisDevice<Kokkos::Cuda>::type> {
public:
  using value_type = ValueType;
  using device_type = typename UseThisDevice<Kokkos::Cuda>::type;
  using numeric_tools_base_type = NumericToolsBase<value_type, device_type>;
  using numeric_tools_serial_type = NumericToolsSerial<value_type, device_type>;
  using numeric_tools_levelset_var0_type = NumericToolsLevelSet<value_type, device_type, 0>;
  using numeric_tools_levelset_var1_type = NumericToolsLevelSet<value_type, device_type, 1>;
  using numeric_tools_levelset_var2_type = NumericToolsLevelSet<value_type, device_type, 2>;
  using numeric_tools_levelset_var3_type = NumericToolsLevelSet<value_type, device_type, 3>;

  TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING;
  TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER;
  TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER;

  void setBaseMember(const ordinal_type method,
                     // input matrix A
                     const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                     // input permutation
                     const ordinal_type_array &perm, const ordinal_type_array &peri,
                     // supernodes
                     const ordinal_type nsupernodes, const ordinal_type_array &supernodes,
                     const size_type_array &gid_ptr, const ordinal_type_array &gid_colidx,
                     const size_type_array &sid_ptr, const ordinal_type_array &sid_colidx,
                     const ordinal_type_array &blk_colidx, const ordinal_type_array &stree_parent,
                     const size_type_array &stree_ptr, const ordinal_type_array &stree_children,
                     const ordinal_type_array_host &stree_level, const ordinal_type_array_host &stree_roots,
                     const ordinal_type verbose) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER;
  }

  void setLevelSetMember(const ordinal_type variant, const ordinal_type device_level_cut,
                         const ordinal_type device_factor_thres, const ordinal_type device_solve_thres,
                         const ordinal_type nstreams) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER;
  }

  void createObject(numeric_tools_base_type *&object) {
    switch (_variant) {
    case 0: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var0_type);
      break;
    }
    case 1: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var1_type);
      break;
    }
    case 2: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var2_type);
      break;
    }
    case 3: {
      if (_method == 1 || _method == 3) {
        TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var3_type);
      } else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid variant input");
      }
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid variant input");
    }
    }
  }
};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <typename ValueType>
class NumericToolsFactory<ValueType, typename UseThisDevice<Kokkos::HIP>::type> {
public:
  using value_type = ValueType;
  using device_type = typename UseThisDevice<Kokkos::HIP>::type;
  using numeric_tools_base_type = NumericToolsBase<value_type, device_type>;
  using numeric_tools_serial_type = NumericToolsSerial<value_type, device_type>;
  using numeric_tools_levelset_var0_type = NumericToolsLevelSet<value_type, device_type, 0>;
  using numeric_tools_levelset_var1_type = NumericToolsLevelSet<value_type, device_type, 1>;
  using numeric_tools_levelset_var2_type = NumericToolsLevelSet<value_type, device_type, 2>;
  using numeric_tools_levelset_var3_type = NumericToolsLevelSet<value_type, device_type, 3>;

  TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING;
  TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER;
  TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER;

  void setBaseMember(const ordinal_type method,
                     // input matrix A
                     const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                     // input permutation
                     const ordinal_type_array &perm, const ordinal_type_array &peri,
                     // supernodes
                     const ordinal_type nsupernodes, const ordinal_type_array &supernodes,
                     const size_type_array &gid_ptr, const ordinal_type_array &gid_colidx,
                     const size_type_array &sid_ptr, const ordinal_type_array &sid_colidx,
                     const ordinal_type_array &blk_colidx, const ordinal_type_array &stree_parent,
                     const size_type_array &stree_ptr, const ordinal_type_array &stree_children,
                     const ordinal_type_array_host &stree_level, const ordinal_type_array_host &stree_roots,
                     const ordinal_type verbose) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER;
  }

  void setLevelSetMember(const ordinal_type variant, const ordinal_type device_level_cut,
                         const ordinal_type device_factor_thres, const ordinal_type device_solve_thres,
                         const ordinal_type nstreams) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER;
  }

  void createObject(numeric_tools_base_type *&object) {
    switch (_variant) {
    case 0: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var0_type);
      break;
    }
    case 1: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var1_type);
      break;
    }
    case 2: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var2_type);
      break;
    }
    case 3: {
      if (_method == 1 || _method == 3) {
        TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var3_type);
      } else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid variant input");
      }
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid variant input");
      break;
    }
    }
  }
};
#endif

#if defined(KOKKOS_ENABLE_SYCL)
template <typename ValueType>
class NumericToolsFactory<ValueType, typename UseThisDevice<Kokkos::Experimental::SYCL>::type> {
public:
  using value_type = ValueType;
  using device_type = typename UseThisDevice<Kokkos::Experimental::SYCL>::type;
  using numeric_tools_base_type = NumericToolsBase<value_type, device_type>;
  using numeric_tools_serial_type = NumericToolsSerial<value_type, device_type>;
  using numeric_tools_levelset_var0_type = NumericToolsLevelSet<value_type, device_type, 0>;
  using numeric_tools_levelset_var1_type = NumericToolsLevelSet<value_type, device_type, 1>;
  using numeric_tools_levelset_var2_type = NumericToolsLevelSet<value_type, device_type, 2>;
  using numeric_tools_levelset_var3_type = NumericToolsLevelSet<value_type, device_type, 3>;

  TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING;
  TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER;
  TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER;

  void setBaseMember(const ordinal_type method,
                     // input matrix A
                     const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                     // input permutation
                     const ordinal_type_array &perm, const ordinal_type_array &peri,
                     // supernodes
                     const ordinal_type nsupernodes, const ordinal_type_array &supernodes,
                     const size_type_array &gid_ptr, const ordinal_type_array &gid_colidx,
                     const size_type_array &sid_ptr, const ordinal_type_array &sid_colidx,
                     const ordinal_type_array &blk_colidx, const ordinal_type_array &stree_parent,
                     const size_type_array &stree_ptr, const ordinal_type_array &stree_children,
                     const ordinal_type_array_host &stree_level, const ordinal_type_array_host &stree_roots,
                     const ordinal_type verbose) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER;
  }

  void setLevelSetMember(const ordinal_type variant, const ordinal_type device_level_cut,
                         const ordinal_type device_factor_thres, const ordinal_type device_solve_thres,
                         const ordinal_type nstreams) {
    TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER;
  }

  void createObject(numeric_tools_base_type *&object) {
    switch (_variant) {
    case 0: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var0_type);
      break;
    }
    case 1: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var1_type);
      break;
    }
    case 2: {
      TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_BODY(numeric_tools_levelset_var2_type);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid variant input");
      break;
    }
    }
  }
};
#endif

#undef TACHO_NUMERIC_TOOLS_FACTORY_BASE_USING
#undef TACHO_NUMERIC_TOOLS_FACTORY_BASE_MEMBER
#undef TACHO_NUMERIC_TOOLS_FACTORY_SET_BASE_MEMBER
#undef TACHO_NUMERIC_TOOLS_FACTORY_LEVELSET_MEMBER
#undef TACHO_NUMERIC_TOOLS_FACTORY_SET_LEVELSET_MEMBER
#undef TACHO_NUMERIC_TOOLS_SERIAL_BODY

} // namespace Tacho
#endif
