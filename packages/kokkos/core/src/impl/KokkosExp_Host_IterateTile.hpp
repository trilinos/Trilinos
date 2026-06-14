// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HOST_EXP_ITERATE_TILE_HPP
#define KOKKOS_HOST_EXP_ITERATE_TILE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP) && !defined(__CUDA_ARCH__)
#define KOKKOS_MDRANGE_IVDEP
#endif

#ifdef KOKKOS_MDRANGE_IVDEP
#define KOKKOS_ENABLE_IVDEP_MDRANGE _Pragma("ivdep")
#else
#define KOKKOS_ENABLE_IVDEP_MDRANGE
#endif

#include <algorithm>

namespace Kokkos {
namespace Impl {

// FIXME: Implement generation of nested for loops for parallel_reduce
// using recursive templates (like it is done for parallel_for)

// Loop Macros for parallel_reduce

// non-tagged

#define KOKKOS_IMPL_APPLY_REDUX(val, func, ...) func(__VA_ARGS__, val);

// LayoutRight
// d = 0 to start
#define KOKKOS_IMPL_LOOP_R_1_REDUX(val, func, type, m_offset, extent, d, ...) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                 \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {          \
    KOKKOS_IMPL_APPLY_REDUX(val, func, __VA_ARGS__, i0 + m_offset[d])         \
  }

#define KOKKOS_IMPL_LOOP_R_2_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_LOOP_R_1_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i1 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_3_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_LOOP_R_2_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i2 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_4_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_LOOP_R_3_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i3 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_5_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_LOOP_R_4_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i4 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_6_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_LOOP_R_5_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i5 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_7_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_LOOP_R_6_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i6 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_8_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_LOOP_R_7_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i7 + m_offset[d])                 \
  }

// LayoutLeft
// d = rank-1 to start
#define KOKKOS_IMPL_LOOP_L_1_REDUX(val, func, type, m_offset, extent, d, ...) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                 \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {          \
    KOKKOS_IMPL_APPLY_REDUX(val, func, i0 + m_offset[d], __VA_ARGS__)         \
  }

#define KOKKOS_IMPL_LOOP_L_2_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_LOOP_L_1_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i1 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_3_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_LOOP_L_2_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i2 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_4_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_LOOP_L_3_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i3 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_5_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_LOOP_L_4_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i4 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_6_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_LOOP_L_5_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i5 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_7_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_LOOP_L_6_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i6 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_8_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_LOOP_L_7_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i7 + m_offset[d], __VA_ARGS__)                 \
  }

// Left vs Right
#define KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX(val, func, type, is_left, m_offset, \
                                        extent, rank)                       \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                               \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[0]); ++i0) {        \
    KOKKOS_IMPL_APPLY_REDUX(val, func, i0 + m_offset[0])                    \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[rank - 1]); ++i1) { \
      KOKKOS_IMPL_LOOP_L_1_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i1 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[0]); ++i1) {        \
      KOKKOS_IMPL_LOOP_R_1_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i1 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[rank - 1]); ++i2) { \
      KOKKOS_IMPL_LOOP_L_2_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i2 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[0]); ++i2) {        \
      KOKKOS_IMPL_LOOP_R_2_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i2 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[rank - 1]); ++i3) { \
      KOKKOS_IMPL_LOOP_L_3_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i3 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[0]); ++i3) {        \
      KOKKOS_IMPL_LOOP_R_3_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i3 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[rank - 1]); ++i4) { \
      KOKKOS_IMPL_LOOP_L_4_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i4 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[0]); ++i4) {        \
      KOKKOS_IMPL_LOOP_R_4_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i4 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[rank - 1]); ++i5) { \
      KOKKOS_IMPL_LOOP_L_5_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i5 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[0]); ++i5) {        \
      KOKKOS_IMPL_LOOP_R_5_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i5 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[rank - 1]); ++i6) { \
      KOKKOS_IMPL_LOOP_L_6_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i6 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[0]); ++i6) {        \
      KOKKOS_IMPL_LOOP_R_6_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i6 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[rank - 1]); ++i7) { \
      KOKKOS_IMPL_LOOP_L_7_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i7 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[0]); ++i7) {        \
      KOKKOS_IMPL_LOOP_R_7_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i7 + m_offset[0])                            \
    }                                                                         \
  }

// Partial vs Full Tile
#define KOKKOS_IMPL_TILE_LOOP_1_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_2_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_3_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_4_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_5_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_6_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_7_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_8_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }
// end untagged macros

// tagged

#define KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, ...) \
  func(tag, __VA_ARGS__, val);

// LayoutRight
// d = 0 to start
#define KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX(val, tag, func, type, m_offset, \
                                          extent, d, ...)                 \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                             \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {      \
    KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, __VA_ARGS__,           \
                                   i0 + m_offset[d])                      \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i1 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i2 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i3 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i4 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i5 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i6 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_8_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i7 + m_offset[d])   \
  }

// LayoutLeft
// d = rank-1 to start
#define KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX(val, tag, func, type, m_offset, \
                                          extent, d, ...)                 \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                             \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {      \
    KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, i0 + m_offset[d],      \
                                   __VA_ARGS__)                           \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i1 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i2 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i3 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i4 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i5 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i6 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_8_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i7 + m_offset[d], __VA_ARGS__)   \
  }

// Left vs Right
#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                 \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[0]); ++i0) {          \
    KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, i0 + m_offset[0])          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[rank - 1]); ++i1) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i1 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[0]); ++i1) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i1 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[rank - 1]); ++i2) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i2 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[0]); ++i2) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i2 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[rank - 1]); ++i3) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i3 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[0]); ++i3) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i3 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[rank - 1]); ++i4) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i4 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[0]); ++i4) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i4 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[rank - 1]); ++i5) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i5 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[0]); ++i5) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i5 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[rank - 1]); ++i6) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i6 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[0]); ++i6) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i6 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[rank - 1]); ++i7) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i7 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[0]); ++i7) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i7 + m_offset[0])          \
    }                                                                         \
  }

// Partial vs Full Tile
#define KOKKOS_IMPL_TAGGED_TILE_LOOP_1_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_2_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_3_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_4_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_5_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_6_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_7_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_8_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

// end tagged macros

// Structs for calling loops to compute the reductions

template <int Rank, bool IsLeft, typename IType, typename Tagged,
          typename Enable = void>
struct Tile_Loop_Type;

// non-tagged versions

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<1, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_1_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 1);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<2, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_2_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 2);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<3, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_3_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 3);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<4, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_4_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 4);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<5, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_5_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 5);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<6, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_6_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 6);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<7, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_7_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 7);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<8, IsLeft, IType, void, void> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_8_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 8);
  }
};

// tagged versions

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<1, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_1_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 1);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<2, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_2_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 2);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<3, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_3_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 3);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<4, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_4_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 4);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<5, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_5_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 5);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<6, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_6_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 6);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<7, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_7_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 7);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<8, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void_v<Tagged>>> {
  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_8_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 8);
  }
};

// end structs for calling loops for reduction

template <typename RP, typename Functor, typename Tag = void,
          typename ValueType = void, typename Enable = void>
struct HostIterateTile;

// For ParallelFor
template <typename RP, typename Functor, typename Tag, typename ValueType>
struct HostIterateTile<RP, Functor, Tag, ValueType,
                       std::enable_if_t<std::is_void_v<ValueType>>> {
  using index_type = typename RP::index_type;
  using point_type = typename RP::point_type;

  using value_type = ValueType;

  inline HostIterateTile(RP const& rp, Functor const& func)
      : m_rp(rp), m_func(func) {}

  inline void check_iteration_bounds(point_type& actual_tile,
                                     const point_type& offset) const {
    for (int i = 0; i < RP::rank; ++i) {
      if ((offset[i] + m_rp.m_tile[i]) <= m_rp.m_upper[i]) {
        actual_tile[i] = m_rp.m_tile[i];
      } else {
        actual_tile[i] =
            m_rp.m_upper[i] - offset[i];  // remaining elements in dimension i
      }
    }
  }  // end check_iteration_bounds

  template <int Rank>
  struct RankTag {
    using type = RankTag<Rank>;
    enum { value = (int)Rank };
  };

  // functor encapsulating the inner-most loop
  template <typename TileOffset, typename TileDims, typename... Idxs>
  void func_innermost_loop(TileOffset const& offset, TileDims const& tiledims,
                           Idxs&&... idxs) const {
    if constexpr (RP::inner_direction == Iterate::Left) {
      KOKKOS_ENABLE_IVDEP_MDRANGE
      for (index_type i = offset[0]; i < offset[0] + tiledims[0]; ++i) {
        if constexpr (std::is_void_v<Tag>) {
          m_func(i, (Idxs&&)idxs...);
        } else {
          m_func(Tag{}, i, (Idxs&&)idxs...);
        }
      }
    } else {
      KOKKOS_ENABLE_IVDEP_MDRANGE
      for (index_type i = offset[RP::rank - 1];
           i < offset[RP::rank - 1] + tiledims[RP::rank - 1]; ++i) {
        if constexpr (std::is_void_v<Tag>) {
          m_func((Idxs&&)idxs..., i);
        } else {
          m_func(Tag{}, (Idxs&&)idxs..., i);
        }
      }
    }
  }

  // ----------------------------------------------------------------------- //
  // \brief Nested loops with recursive template instantiation
  //
  // Nested for loop order depends on the iteration order:
  //  Iterate::Left:
  //    Outermost loop corresponds to right-most index.
  //  Iterate::Right:
  //    Outermost loop corresponds to left-most index.
  // The fastest changing index is always in innermost loop.
  //
  // Indices accumulated in parameter pack Idxs...
  //
  // Loop indices passed according to iteration order:
  //  Iterate::Left:
  //    functor(i_{inner_most}, ..., i_{outer_most})
  //  Iterate::Right:
  //    functor(i_{outer_most}, ..., i_{inner_most})
  //
  //  \tparam IterLevel iteration level of the nested loops
  //  \tparam Idxs... index pack
  template <unsigned IterLevel, typename TileOffset, typename TileDims,
            typename... Idxs>
  inline void iterate(std::integral_constant<unsigned, IterLevel>,
                      TileOffset const& offset, TileDims const& tiledims,
                      Idxs... idxs) const {
    const index_type start = (RP::inner_direction == Iterate::Left)
                                 ? offset[RP::rank - 1 - IterLevel]
                                 : offset[IterLevel];
    const index_type end   = (RP::inner_direction == Iterate::Left)
                                 ? offset[RP::rank - 1 - IterLevel] +
                                     tiledims[RP::rank - 1 - IterLevel]
                                 : offset[IterLevel] + tiledims[IterLevel];

    for (index_type idx = start; idx < end; ++idx) {
      if constexpr (RP::inner_direction == Iterate::Left) {
        iterate(std::integral_constant<unsigned, IterLevel + 1>(), offset,
                tiledims, idx, idxs...);
      } else {
        iterate(std::integral_constant<unsigned, IterLevel + 1>(), offset,
                tiledims, idxs..., idx);
      }
    }
  }

  template <typename TileOffset, typename TileDims, typename... Idxs>
  inline void iterate(std::integral_constant<unsigned, RP::rank - 1>,
                      TileOffset const& offset, TileDims const& tiledims,
                      Idxs... idxs) const {
    func_innermost_loop(offset, tiledims, idxs...);
  }

  template <typename IType>
  inline void operator()(IType tile_idx) const {
    point_type offset;
    point_type tiledims;

    if constexpr (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - return actual tile dims in tiledims
    check_iteration_bounds(tiledims, offset);

    iterate(std::integral_constant<unsigned, 0u>(), offset, tiledims);
  }

  RP const m_rp;
  Functor const m_func;
  std::conditional_t<std::is_void_v<Tag>, int, Tag> m_tag{};
};

// For ParallelReduce
// ValueType - scalar: For reductions
template <typename RP, typename Functor, typename Tag, typename ValueType>
struct HostIterateTile<RP, Functor, Tag, ValueType,
                       std::enable_if_t<!std::is_void_v<ValueType> &&
                                        !std::is_array_v<ValueType>>> {
  using index_type = typename RP::index_type;
  using point_type = typename RP::point_type;

  using value_type = ValueType;

  inline HostIterateTile(RP const& rp, Functor const& func)
      : m_rp(rp)  // Cuda 7.0 does not like braces...
        ,
        m_func(func) {
    // Errors due to braces rather than parenthesis for init (with cuda 7.0)
    //      /home/ndellin/kokkos/core/src/impl/KokkosExp_Host_IterateTile.hpp:1216:98:
    //      error: too many braces around initializer for ‘int’ [-fpermissive]
    //      /home/ndellin/kokkos/core/src/impl/KokkosExp_Host_IterateTile.hpp:1216:98:
    //      error: aggregate value used where an integer was expected
  }

  inline bool check_iteration_bounds(point_type& partial_tile,
                                     const point_type& offset) const {
    bool is_full_tile = true;

    for (int i = 0; i < RP::rank; ++i) {
      if ((offset[i] + m_rp.m_tile[i]) <= m_rp.m_upper[i]) {
        partial_tile[i] = m_rp.m_tile[i];
      } else {
        is_full_tile = false;
        partial_tile[i] =
            m_rp.m_upper[i] - offset[i];  // remaining elements in dimension i
      }
    }

    return is_full_tile;
  }  // end check bounds

  template <int Rank>
  struct RankTag {
    using type = RankTag<Rank>;
    enum { value = (int)Rank };
  };

  template <typename IType>
  inline void operator()(IType tile_idx, value_type& val) const {
    point_type m_offset;
    point_type m_tiledims;

    if constexpr (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    Tile_Loop_Type<RP::rank, (RP::inner_direction == Iterate::Left), index_type,
                   Tag>::apply(val, m_func.get_functor(), full_tile, m_offset,
                               m_rp.m_tile, m_tiledims);
  }

  RP const m_rp;
  Functor const m_func;
  std::conditional_t<std::is_void_v<Tag>, int, Tag> m_tag{};
};

// For ParallelReduce
// Extra specialization for array reductions
// ValueType[]: For array reductions
template <typename RP, typename Functor, typename Tag, typename ValueType>
struct HostIterateTile<RP, Functor, Tag, ValueType,
                       std::enable_if_t<!std::is_void_v<ValueType> &&
                                        std::is_array_v<ValueType>>> {
  using index_type = typename RP::index_type;
  using point_type = typename RP::point_type;

  using value_type =
      std::remove_extent_t<ValueType>;  // strip away the
                                        // 'array-ness' [], only
                                        // underlying type remains

  inline HostIterateTile(RP const& rp, Functor const& func)
      : m_rp(rp)  // Cuda 7.0 does not like braces...
        ,
        m_func(func) {}

  inline bool check_iteration_bounds(point_type& partial_tile,
                                     const point_type& offset) const {
    bool is_full_tile = true;

    for (int i = 0; i < RP::rank; ++i) {
      if ((offset[i] + m_rp.m_tile[i]) <= m_rp.m_upper[i]) {
        partial_tile[i] = m_rp.m_tile[i];
      } else {
        is_full_tile = false;
        partial_tile[i] =
            m_rp.m_upper[i] - offset[i];  // remaining elements in dimension i
      }
    }

    return is_full_tile;
  }  // end check bounds

  template <int Rank>
  struct RankTag {
    using type = RankTag<Rank>;
    enum { value = (int)Rank };
  };

  template <typename IType>
  inline void operator()(IType tile_idx, value_type* val) const {
    point_type m_offset;
    point_type m_tiledims;

    if constexpr (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    Tile_Loop_Type<RP::rank, (RP::inner_direction == Iterate::Left), index_type,
                   Tag>::apply(val, m_func, full_tile, m_offset, m_rp.m_tile,
                               m_tiledims);
  }

  RP const m_rp;
  Functor const m_func;
  std::conditional_t<std::is_void_v<Tag>, int, Tag> m_tag{};
};

// ------------------------------------------------------------------ //

#undef KOKKOS_IMPL_APPLY_REDUX
#undef KOKKOS_IMPL_LOOP_R_1_REDUX
#undef KOKKOS_IMPL_LOOP_R_2_REDUX
#undef KOKKOS_IMPL_LOOP_R_3_REDUX
#undef KOKKOS_IMPL_LOOP_R_4_REDUX
#undef KOKKOS_IMPL_LOOP_R_5_REDUX
#undef KOKKOS_IMPL_LOOP_R_6_REDUX
#undef KOKKOS_IMPL_LOOP_R_7_REDUX
#undef KOKKOS_IMPL_LOOP_R_8_REDUX
#undef KOKKOS_IMPL_LOOP_L_1_REDUX
#undef KOKKOS_IMPL_LOOP_L_2_REDUX
#undef KOKKOS_IMPL_LOOP_L_3_REDUX
#undef KOKKOS_IMPL_LOOP_L_4_REDUX
#undef KOKKOS_IMPL_LOOP_L_5_REDUX
#undef KOKKOS_IMPL_LOOP_L_6_REDUX
#undef KOKKOS_IMPL_LOOP_L_7_REDUX
#undef KOKKOS_IMPL_LOOP_L_8_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_1_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_2_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_3_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_4_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_5_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_6_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_7_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_8_REDUX
#undef KOKKOS_IMPL_TAGGED_APPLY_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_8_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_8_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_1_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_2_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_3_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_4_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_5_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_6_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_7_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_8_REDUX

}  // namespace Impl
}  // namespace Kokkos

#endif
