// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_ENTITYVALUES_HPP
#define STK_ENTITYVALUES_HPP

#include "stk_util/stk_config.h"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_mesh/base/FieldIndexTypes.hpp"
#include "Kokkos_Macros.hpp"

namespace stk::mesh {

//==============================================================================
// Device EntityValues
//==============================================================================

template <typename T, typename Space = stk::ngp::HostSpace,
          Layout DataLayout = DefaultLayoutSelector<Space>::layout>
class EntityValues
{
public:
  using value_type = T;
  using space = Space;
  using exec_space = typename Space::exec_space;
  using mem_space = typename Space::mem_space;
  static constexpr Layout layout = DataLayout;

  KOKKOS_INLINE_FUNCTION EntityValues(T* dataPtr, int numComponents, int numCopies, int scalarStride,
                                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_scalarStride(scalarStride)
  {
    static_assert(DataLayout == Layout::Left, "Only Layout::Left is supported for device data");
  }

  KOKKOS_DEFAULTED_FUNCTION ~EntityValues() = default;

  // The functions below may be used to iterate through your Entity data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
  //     for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::CopyIdx copy : entityValues.copies()) {
  //     for (stk::mesh::ComponentIdx component : entityValues.components()) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may use hard-coded numerical
  // indexing through user-defined literals, such as: 0_copy, 1_copy, 2_copy, etc. and
  // 0_comp, 1_comp, 2_comp, etc.

  KOKKOS_INLINE_FUNCTION int num_copies() const { return m_numCopies; }
  KOKKOS_INLINE_FUNCTION CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  KOKKOS_INLINE_FUNCTION int num_components() const { return m_numComponents; }
  KOKKOS_INLINE_FUNCTION ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  KOKKOS_INLINE_FUNCTION int num_scalars() const { return m_numCopies*m_numComponents; }
  KOKKOS_INLINE_FUNCTION ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  KOKKOS_INLINE_FUNCTION bool is_field_defined() const { return m_numComponents != 0; }

  KOKKOS_INLINE_FUNCTION T& operator()(
      const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);

    return *m_dataPtr;
  }

  KOKKOS_INLINE_FUNCTION T& operator()(ComponentIdx component,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[component*m_scalarStride];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(CopyIdx copy,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[copy*m_scalarStride];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(CopyIdx copy, ComponentIdx component,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[(copy()*m_numComponents + component())*m_scalarStride];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(ScalarIdx scalar,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[scalar*m_scalarStride];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Entity data on either the host
  // or device if you follow this pattern:
  //
  //   double* entityPtr = entityValues.pointer();
  //   const int copyStride = entityValues.copy_stride();
  //   const int componentStride = entityValues.component_stride();
  //
  //   for (int copy = 0; copy < entityValues.num_copies(); ++copy) {
  //     for (int component = 0; component < entityValues.num_components(); ++component) {
  //       entityPtr[copy*copyStride + component*componentStride] = 0.0;
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  KOKKOS_INLINE_FUNCTION T* pointer() const { return m_dataPtr; }

  KOKKOS_INLINE_FUNCTION int copy_stride() const { return m_numComponents*m_scalarStride; }

  KOKKOS_INLINE_FUNCTION int component_stride() const { return m_scalarStride; }

  KOKKOS_INLINE_FUNCTION int scalar_stride() const { return m_scalarStride; }

private:
#ifdef STK_FIELD_BOUNDS_CHECK
  KOKKOS_INLINE_FUNCTION void check_defined_field(const char* file, int line) const {
    if (not is_field_defined()) {
      if (line == -1) {
        printf("Error: Accessing EntityValues for Field '%s' that is not defined on this Entity.\n", m_fieldName);
      }
      else {
        printf("Error: %s:%i: Accessing EntityValues for Field '%s' that is not defined on this Entity.\n",
               file, line, m_fieldName);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_single_scalar_access(const char* file, int line) const {
    if (m_numComponents*m_numCopies != 1) {
      if (line == -1) {
        printf("Error: Accessing EntityValues for Field '%s' as a scalar when it has %i components and %i"
               " copies.  Please use an Entity operator() that takes appropriate index arguments.\n",
               m_fieldName, m_numComponents, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Accessing EntityValues for Field '%s' as a scalar when it has %i components and %i"
               " copies.  Please use an Entity operator() that takes appropriate index arguments.\n",
               file, line, m_fieldName, m_numComponents, m_numCopies);

      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_single_component_access(const char* file, int line) const {
    if (m_numComponents != 1) {
      if (line == -1) {
        printf("Error: Accessing EntityValues for Field '%s' as if it only has one component when it actually"
               " has %i components.  Please use an Entity operator() that also has a component argument.\n",
               m_fieldName, m_numComponents);
      }
      else {
        printf("Error: %s:%i: Accessing EntityValues for Field '%s' as if it only has one component when it actually"
               " has %i components.  Please use an Entity operator() that also has a component argument.\n",
               file, line, m_fieldName, m_numComponents);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_single_copy_access(const char* file, int line) const {
    if (m_numCopies != 1) {
      if (line == -1) {
        printf("Error: Accessing EntityValues for Field '%s' as if it only has one copy when it actually has"
               " %i copies.  Please use an Entity operator() that also has a copy argument.\n",
               m_fieldName, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Accessing EntityValues for Field '%s' as if it only has one copy when it actually has"
               " %i copies.  Please use an Entity operator() that also has a copy argument.\n",
               file, line, m_fieldName, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_component_bounds(int component, const char* file, int line) const {
    if ((component < 0) || (component >= m_numComponents)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to EntityValues for Field '%s' with component index %i for an"
               " Entity with %i components.\n", m_fieldName, component, m_numComponents);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to EntityValues for Field '%s' with component index %i for an"
               " Entity with %i components.\n", file, line, m_fieldName, component, m_numComponents);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_copy_bounds(int copy, const char* file, int line) const {
    if ((copy < 0) || (copy >= m_numCopies)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to EntityValues for Field '%s' with copy index %i for an Entity"
               " with %i copies.\n", m_fieldName, copy, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to EntityValues for Field '%s' with copy index %i for an Entity"
               " with %i copies.\n", file, line, m_fieldName, copy, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_copy_and_component_bounds(int copy, int component,
                                                              const char* file, int line) const {
    if (((copy < 0) || (copy >= m_numCopies)) || ((component < 0) || (component >= m_numComponents))) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to EntityValues for Field '%s' with component index %i and copy"
               " index %i for an Entity with %i components and %i copies.\n", m_fieldName, component, copy,
               m_numComponents, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to EntityValues for Field '%s' with component index %i and copy"
               " index %i for an Entity with %i components and %i copies.\n", file, line, m_fieldName, component, copy,
               m_numComponents, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_scalar_bounds(int scalar, const char* file, int line) const {
    if ((scalar < 0) || (scalar >= m_numCopies*m_numComponents)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to EntityValues for Field '%s' with scalar index %i for an"
               " Entity with %i scalars.\n", m_fieldName, scalar, m_numCopies*m_numComponents);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to EntityValues for Field '%s' with scalar index %i for an"
               " Entity with %i scalars.\n", file, line, m_fieldName, scalar, m_numCopies*m_numComponents);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
#else
  KOKKOS_INLINE_FUNCTION void check_defined_field(const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_single_scalar_access(const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_single_component_access(const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_single_copy_access(const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_component_bounds(int, const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_copy_bounds(int, const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_copy_and_component_bounds(int, int, const char*, int) const {}
  KOKKOS_INLINE_FUNCTION void check_scalar_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.
  KOKKOS_INLINE_FUNCTION T& operator()(int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_scalarStride;
};


//==============================================================================
// Host EntityValues: Layout::Right
//==============================================================================

template <typename T>
class EntityValues<T, stk::ngp::HostSpace, Layout::Right>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Right;

  inline EntityValues(T* dataPtr, int numComponents, int numCopies, [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies)
  {}

  ~EntityValues() = default;

  // The functions below may be used to iterate through your Entity data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
  //     for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::CopyIdx copy : entityValues.copies()) {
  //     for (stk::mesh::ComponentIdx component : entityValues.components()) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may use hard-coded numerical
  // indexing through user-defined literals, such as: 0_copy, 1_copy, 2_copy, etc. and
  // 0_comp, 1_comp, 2_comp, etc.

  inline int num_copies() const { return m_numCopies; }
  inline CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  inline int num_components() const { return m_numComponents; }
  inline ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  inline int num_scalars() const { return m_numCopies*m_numComponents; }
  inline ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  inline bool is_field_defined() const { return m_numComponents != 0; }

  inline T& operator()(
      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);

    return *m_dataPtr;
  }

  inline T& operator()(ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[component];
  }

  inline T& operator()(CopyIdx copy,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[copy];
  }

  inline T& operator()(CopyIdx copy, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[copy()*m_numComponents + component()];
  }

  inline T& operator()(ScalarIdx scalar,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[scalar];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Entity data on either the host
  // or device if you follow this pattern:
  //
  //   double* entityPtr = entityValues.pointer();
  //   const int copyStride = entityValues.copy_stride();
  //   const int componentStride = entityValues.component_stride();
  //
  //   for (int copy = 0; copy < entityValues.num_copies(); ++copy) {
  //     for (int component = 0; component < entityValues.num_components(); ++component) {
  //       entityPtr[copy*copyStride + component*componentStride] = 0.0;
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  inline T* pointer() const { return m_dataPtr; }

  inline int copy_stride() const { return m_numComponents; }

  inline int component_stride() const { return 1; }

  inline int scalar_stride() const { return 1; }

private:
#ifdef STK_FIELD_BOUNDS_CHECK
  inline std::string location_string(const char* file, int line) const {
    if (line != -1) {
      std::string fileName(file);
      std::size_t pathDelimeter = fileName.find_last_of("/");
      if (pathDelimeter < fileName.size()) {
        fileName = fileName.substr(pathDelimeter+1);
      }
      return fileName + ":" + std::to_string(line) + ": ";
    }
    else {
      return "";
    }
  }
  inline void check_defined_field(const char* file, int line) const {
    STK_ThrowRequireMsg(is_field_defined(),
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " that is not defined on this Entity.");
  }
  inline void check_single_scalar_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents*m_numCopies == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as a scalar when it has " << m_numComponents << " components and " << m_numCopies <<
                        " copies.  Please use an Entity operator() that takes appropriate index arguments.");
  }
  inline void check_single_component_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as if it only has one component when it actually has " << m_numComponents <<
                        " components.  Please use an Entity operator() that also has a component argument.");
  }
  inline void check_single_copy_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numCopies == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as if it only has one copy when it actually has " << m_numCopies <<
                        " copies.  Please use an Entity operator() that also has a copy argument.");
  }
  inline void check_component_bounds(int component, const char* file, int line) const {
    STK_ThrowRequireMsg((component >= 0) && (component < m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with component index " << component << " for an Entity with " <<
                        m_numComponents << " components.");
  }
  inline void check_copy_bounds(int copy, const char* file, int line) const {
    STK_ThrowRequireMsg((copy >= 0) && (copy < m_numCopies),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with copy index " << copy << " for an Entity with " << m_numCopies <<
                        " copies.");
  }
  inline void check_copy_and_component_bounds(int copy, int component, const char* file, int line) const {
    STK_ThrowRequireMsg(((copy >= 0) && (copy < m_numCopies)) && ((component >= 0) && (component < m_numComponents)),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with component index " << component << " and copy index " << copy <<
                        " for an Entity with " << m_numComponents << " components and " << m_numCopies << " copies.");
  }
  inline void check_scalar_bounds(int scalar, const char* file, int line) const {
    STK_ThrowRequireMsg((scalar >= 0) && (scalar < m_numCopies*m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with scalar index " << scalar << " for an Entity with " <<
                        m_numCopies*m_numComponents << " scalars.");
  }
#else
  inline void check_defined_field(const char*, int) const {}
  inline void check_single_scalar_access(const char*, int) const {}
  inline void check_single_component_access(const char*, int) const {}
  inline void check_single_copy_access(const char*, int) const {}
  inline void check_component_bounds(int, const char*, int) const {}
  inline void check_copy_bounds(int, const char*, int) const {}
  inline void check_copy_and_component_bounds(int, int, const char*, int) const {}
  inline void check_scalar_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.
  inline T& operator()(int) const = delete;
  inline T& operator()(int, int) const = delete;
  inline T& operator()(ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
};


//==============================================================================
// Host EntityValues: Layout::Left
//==============================================================================

template <typename T>
class EntityValues<T, stk::ngp::HostSpace, Layout::Left>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Left;

  inline EntityValues(T* dataPtr, int numComponents, int numCopies, int scalarStride,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_scalarStride(scalarStride)
  {}

  ~EntityValues() = default;

  // The functions below may be used to iterate through your Entity data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
  //     for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::CopyIdx copy : entityValues.copies()) {
  //     for (stk::mesh::ComponentIdx component : entityValues.components()) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may use hard-coded numerical
  // indexing through user-defined literals, such as: 0_copy, 1_copy, 2_copy, etc. and
  // 0_comp, 1_comp, 2_comp, etc.

  inline int num_copies() const { return m_numCopies; }
  inline CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  inline int num_components() const { return m_numComponents; }
  inline ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  inline int num_scalars() const { return m_numCopies*m_numComponents; }
  inline ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  inline bool is_field_defined() const { return m_numComponents != 0; }

  inline T& operator()(
      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);

    return *m_dataPtr;
  }

  inline T& operator()(ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[component*m_scalarStride];
  }

  inline T& operator()(CopyIdx copy,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[copy*m_scalarStride];
  }

  inline T& operator()(CopyIdx copy, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[(copy()*m_numComponents + component())*m_scalarStride];
  }

  inline T& operator()(ScalarIdx scalar,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[scalar*m_scalarStride];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Entity data on either the host
  // or device if you follow this pattern:
  //
  //   double* entityPtr = entityValues.pointer();
  //   const int copyStride = entityValues.copy_stride();
  //   const int componentStride = entityValues.component_stride();
  //
  //   for (int copy = 0; copy < entityValues.num_copies(); ++copy) {
  //     for (int component = 0; component < entityValues.num_components(); ++component) {
  //       entityPtr[copy*copyStride + component*componentStride] = 0.0;
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  inline T* pointer() const { return m_dataPtr; }

  inline int copy_stride() const { return m_numComponents*m_scalarStride; }

  inline int component_stride() const { return m_scalarStride; }

  inline int scalar_stride() const { return m_scalarStride; }

private:
#ifdef STK_FIELD_BOUNDS_CHECK
  inline std::string location_string(const char* file, int line) const {
    if (line != -1) {
      std::string fileName(file);
      std::size_t pathDelimeter = fileName.find_last_of("/");
      if (pathDelimeter < fileName.size()) {
        fileName = fileName.substr(pathDelimeter+1);
      }
      return fileName + ":" + std::to_string(line) + ": ";
    }
    else {
      return "";
    }
  }
  inline void check_defined_field(const char* file, int line) const {
    STK_ThrowRequireMsg(is_field_defined(),
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " that is not defined on this Entity.");
  }
  inline void check_single_scalar_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents*m_numCopies == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as a scalar when it has " << m_numComponents << " components and " << m_numCopies <<
                        " copies.  Please use an Entity operator() that takes appropriate index arguments.");
  }
  inline void check_single_component_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as if it only has one component when it actually has " << m_numComponents <<
                        " components.  Please use an Entity operator() that also has a component argument.");
  }
  inline void check_single_copy_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numCopies == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as if it only has one copy when it actually has " << m_numCopies <<
                        " copies.  Please use an Entity operator() that also has a copy argument.");
  }
  inline void check_component_bounds(int component, const char* file, int line) const {
    STK_ThrowRequireMsg((component >= 0) && (component < m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with component index " << component << " for an Entity with " <<
                        m_numComponents << " components.");
  }
  inline void check_copy_bounds(int copy, const char* file, int line) const {
    STK_ThrowRequireMsg((copy >= 0) && (copy < m_numCopies),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with copy index " << copy << " for an Entity with " << m_numCopies <<
                        " copies.");
  }
  inline void check_copy_and_component_bounds(int copy, int component, const char* file, int line) const {
    STK_ThrowRequireMsg(((copy >= 0) && (copy < m_numCopies)) && ((component >= 0) && (component < m_numComponents)),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with component index " << component << " and copy index " << copy <<
                        " for an Entity with " << m_numComponents << " components and " << m_numCopies << " copies.");
  }
  inline void check_scalar_bounds(int scalar, const char* file, int line) const {
    STK_ThrowRequireMsg((scalar >= 0) && (scalar < m_numCopies*m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with scalar index " << scalar << " for an Entity with " <<
                        m_numCopies*m_numComponents << " scalars.");
  }
#else
  inline void check_defined_field(const char*, int) const {}
  inline void check_single_scalar_access(const char*, int) const {}
  inline void check_single_component_access(const char*, int) const {}
  inline void check_single_copy_access(const char*, int) const {}
  inline void check_component_bounds(int, const char*, int) const {}
  inline void check_copy_bounds(int, const char*, int) const {}
  inline void check_copy_and_component_bounds(int, int, const char*, int) const {}
  inline void check_scalar_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.
  inline T& operator()(int) const = delete;
  inline T& operator()(int, int) const = delete;
  inline T& operator()(ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_scalarStride;
};

//==============================================================================
// Host EntityValues: Layout::Auto
//==============================================================================

template <typename T>
class EntityValues<T, stk::ngp::HostSpace, Layout::Auto>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Auto;

  inline EntityValues(T* dataPtr, int numComponents, int numCopies,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_scalarStride(1),
      m_isLayoutRight(true)
  {}

  inline EntityValues(T* dataPtr, int numComponents, int numCopies, int scalarStride,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_scalarStride(scalarStride),
      m_isLayoutRight(false)
  {}

  ~EntityValues() = default;

  // The functions below may be used to iterate through your Entity data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
  //     for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::CopyIdx copy : entityValues.copies()) {
  //     for (stk::mesh::ComponentIdx component : entityValues.components()) {
  //       entityValues(copy, component) = 0.0;
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may use hard-coded numerical
  // indexing through user-defined literals, such as: 0_copy, 1_copy, 2_copy, etc. and
  // 0_comp, 1_comp, 2_comp, etc.

  inline int num_copies() const { return m_numCopies; }
  inline CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  inline int num_components() const { return m_numComponents; }
  inline ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  inline int num_scalars() const { return m_numCopies*m_numComponents; }
  inline ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  inline bool is_field_defined() const { return m_numComponents != 0; }

  inline T& operator()(
      const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);

    return *m_dataPtr;
  }

  inline T& operator()(ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_component_bounds(component, file, line);

    // Branching improves performance significantly for the Layout::Right case over always doing
    // the more-expensive striding math.  (10% performance hit instead of a 25% performance hit.
    // Layout::Left data has an extra performance hit of 15% instead of 10%, so this roughly
    // normalized the auto-layout penalty between the two layouts.)

    if (m_isLayoutRight) {
      return m_dataPtr[component];
    }
    else {
      return m_dataPtr[component*m_scalarStride];
    }
  }

  inline T& operator()(CopyIdx copy,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_copy_bounds(copy, file, line);

    // Branching improves performance significantly for the Layout::Right case over always doing
    // the more-expensive striding math.  (10% performance hit instead of a 25% performance hit.
    // Layout::Left data has an extra performance hit of 15% instead of 10%, so this roughly
    // normalized the auto-layout penalty between the two layouts.)

    if (m_isLayoutRight) {
      return m_dataPtr[copy];
    }
    else {
      return m_dataPtr[copy*m_scalarStride];
    }
  }

  inline T& operator()(CopyIdx copy, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    // The striding math here is complex enough for both layouts that branching would likely not help.
    // This is very similar to BucketValues access, where branching never helps.

    return m_dataPtr[(copy()*m_numComponents + component())*m_scalarStride];
  }

  inline T& operator()(ScalarIdx scalar,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_scalar_bounds(scalar, file, line);

    // Branching improves performance significantly for the Layout::Right case over always doing
    // the more-expensive striding math.  (10% performance hit instead of a 25% performance hit.
    // Layout::Left data has an extra performance hit of 15% instead of 10%, so this roughly
    // normalized the auto-layout penalty between the two layouts.)

    if (m_isLayoutRight) {
      return m_dataPtr[scalar];
    }
    else {
      return m_dataPtr[scalar*m_scalarStride];
    }
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Entity data on either the host
  // or device if you follow this pattern:
  //
  //   double* entityPtr = entityValues.pointer();
  //   const int copyStride = entityValues.copy_stride();
  //   const int componentStride = entityValues.component_stride();
  //
  //   for (int copy = 0; copy < entityValues.num_copies(); ++copy) {
  //     for (int component = 0; component < entityValues.num_components(); ++component) {
  //       entityPtr[copy*copyStride + component*componentStride] = 0.0;
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  inline T* pointer() const { return m_dataPtr; }

  inline int copy_stride() const { return m_numComponents*m_scalarStride; }

  inline int component_stride() const { return m_scalarStride; }

  inline int scalar_stride() const { return m_scalarStride; }

private:
#ifdef STK_FIELD_BOUNDS_CHECK
  inline std::string location_string(const char* file, int line) const {
    if (line != -1) {
      std::string fileName(file);
      std::size_t pathDelimeter = fileName.find_last_of("/");
      if (pathDelimeter < fileName.size()) {
        fileName = fileName.substr(pathDelimeter+1);
      }
      return fileName + ":" + std::to_string(line) + ": ";
    }
    else {
      return "";
    }
  }
  inline void check_defined_field(const char* file, int line) const {
    STK_ThrowRequireMsg(is_field_defined(),
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " that is not defined on this Entity.");
  }
  inline void check_single_scalar_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents*m_numCopies == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as a scalar when it has " << m_numComponents << " components and " << m_numCopies <<
                        " copies.  Please use an Entity operator() that takes appropriate index arguments.");
  }
  inline void check_single_component_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as if it only has one component when it actually has " << m_numComponents <<
                        " components.  Please use an Entity operator() that also has a component argument.");
  }
  inline void check_single_copy_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numCopies == 1,
                        location_string(file, line) << "Accessing EntityValues for Field '" << m_fieldName << "'"
                        " as if it only has one copy when it actually has " << m_numCopies <<
                        " copies.  Please use an Entity operator() that also has a copy argument.");
  }
  inline void check_component_bounds(int component, const char* file, int line) const {
    STK_ThrowRequireMsg((component >= 0) && (component < m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with component index " << component << " for an Entity with " <<
                        m_numComponents << " components.");
  }
  inline void check_copy_bounds(int copy, const char* file, int line) const {
    STK_ThrowRequireMsg((copy >= 0) && (copy < m_numCopies),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with copy index " << copy << " for an Entity with " << m_numCopies <<
                        " copies.");
  }
  inline void check_copy_and_component_bounds(int copy, int component, const char* file, int line) const {
    STK_ThrowRequireMsg(((copy >= 0) && (copy < m_numCopies)) && ((component >= 0) && (component < m_numComponents)),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with component index " << component << " and copy index " << copy <<
                        " for an Entity with " << m_numComponents << " components and " << m_numCopies << " copies.");
  }
  inline void check_scalar_bounds(int scalar, const char* file, int line) const {
    STK_ThrowRequireMsg((scalar >= 0) && (scalar < m_numCopies*m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to EntityValues for Field '" <<
                        m_fieldName << "' with scalar index " << scalar << " for an Entity with " <<
                        m_numCopies*m_numComponents << " scalars.");
  }
#else
  inline void check_defined_field(const char*, int) const {}
  inline void check_single_scalar_access(const char*, int) const {}
  inline void check_single_component_access(const char*, int) const {}
  inline void check_single_copy_access(const char*, int) const {}
  inline void check_component_bounds(int, const char*, int) const {}
  inline void check_copy_bounds(int, const char*, int) const {}
  inline void check_copy_and_component_bounds(int, int, const char*, int) const {}
  inline void check_scalar_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.
  inline T& operator()(int) const = delete;
  inline T& operator()(int, int) const = delete;
  inline T& operator()(ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_scalarStride;
  bool m_isLayoutRight;
};

}

#endif // STK_ENTITYVALUES_HPP
