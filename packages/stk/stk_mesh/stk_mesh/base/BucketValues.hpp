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

#ifndef STK_BUCKETVALUES_HPP
#define STK_BUCKETVALUES_HPP

#include "stk_util/stk_config.h"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_mesh/base/FieldIndexTypes.hpp"
#include "Kokkos_Macros.hpp"

namespace stk::mesh {

//==============================================================================
// Device BucketValues
//==============================================================================

template <typename T, typename Space = stk::ngp::HostSpace,
          Layout DataLayout = DefaultLayoutSelector<Space>::layout>
class BucketValues
{
public:
  using value_type = T;
  using space = Space;
  using exec_space = typename Space::exec_space;
  using mem_space = typename Space::mem_space;
  static constexpr Layout layout = DataLayout;

  KOKKOS_INLINE_FUNCTION BucketValues(T* dataPtr, int numComponents, int numCopies, int numEntities,
                                      int scalarStride, [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_numEntities(numEntities),
      m_scalarStride(scalarStride)
  {}

  KOKKOS_DEFAULTED_FUNCTION ~BucketValues() = default;

  // The functions below may be used to iterate through your Bucket data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::EntityIdx entity(0); entity < bucketValues.num_entities(); ++entity) {
  //     for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
  //       for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::EntityIdx entity : bucketValues.entities()) {
  //     for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
  //       for (stk::mesh::ComponentIdx component : bucketValues.components()) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may also reuse the Entity
  // index between multiple Fields if they are from the same Bucket.  You may use hard-coded
  // numerical indexing through user-defined literals, such as: 0_entity, 1_entity, 2_entity, etc.,
  // 0_copy, 1_copy, 2_copy, etc., and 0_comp, 1_comp, 2_comp, etc.

  KOKKOS_INLINE_FUNCTION int num_entities() const { return m_numEntities; }
  KOKKOS_INLINE_FUNCTION EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  KOKKOS_INLINE_FUNCTION int num_copies() const { return m_numCopies; }
  KOKKOS_INLINE_FUNCTION CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  KOKKOS_INLINE_FUNCTION int num_components() const { return m_numComponents; }
  KOKKOS_INLINE_FUNCTION ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  KOKKOS_INLINE_FUNCTION int num_scalars() const { return m_numCopies*m_numComponents; }
  KOKKOS_INLINE_FUNCTION ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  KOKKOS_INLINE_FUNCTION bool is_field_defined() const { return m_numComponents != 0; }

  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx entity,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);
    check_entity_bounds(entity, file, line);

    return m_dataPtr[entity];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx entity, ComponentIdx component,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_entity_bounds(entity, file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[component()*m_scalarStride + entity()];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx entity, CopyIdx copy,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[copy()*m_scalarStride + entity()];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx entity, CopyIdx copy, ComponentIdx component,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[(copy()*m_numComponents + component())*m_scalarStride + entity()];
  }

  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx entity, ScalarIdx scalar,
                                       const char* file = STK_DEVICE_FILE, int line = STK_DEVICE_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[scalar()*m_scalarStride + entity()];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Bucket data on either the host
  // or device if you follow this pattern:
  //
  //   double* bucketPtr = bucketValues.pointer();
  //   const int entityStride = bucketValues.entity_stride();
  //   const int copyStride = bucketValues.copy_stride();
  //   const int componentStride = bucketValues.component_stride();
  //
  //   for (int entity = 0; entity < bucketValues.num_entities(); ++entity) {
  //     for (int copy = 0; copy < bucketValues.num_copies(); ++copy) {
  //       for (int component = 0; component < bucketValues.num_components(); ++component) {
  //         entityPtr[entity*entityStride + copy*copyStride + component*componentStride] = 0.0;
  //       }
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  KOKKOS_INLINE_FUNCTION T* pointer() const { return m_dataPtr; }

  KOKKOS_INLINE_FUNCTION int entity_stride() const { return 1; }

  KOKKOS_INLINE_FUNCTION int copy_stride() const { return m_numComponents*m_scalarStride; }

  KOKKOS_INLINE_FUNCTION int component_stride() const { return m_scalarStride; }

  KOKKOS_INLINE_FUNCTION int scalar_stride() const { return m_scalarStride; }

private:
#ifdef STK_FIELD_BOUNDS_CHECK
  KOKKOS_INLINE_FUNCTION void check_defined_field(const char* file, int line) const {
    if (not is_field_defined()) {
      if (line == -1) {
        printf("Error: Accessing BucketValues for Field '%s' that is not defined on this Bucket.\n", m_fieldName);
      }
      else {
        printf("Error: %s:%i: Accessing BucketValues for Field '%s' that is not defined on this Bucket.\n",
               file, line, m_fieldName);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_single_scalar_access(const char* file, int line) const {
    if (m_numCopies*m_numComponents != 1) {
      if (line == -1) {
        printf("Error: Accessing BucketValues for Field '%s' as a scalar when it has %i components and %i copies. "
               " Please use a Bucket operator() that takes appropriate index arguments.\n", m_fieldName,
               m_numComponents, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Accessing BucketValues Field '%s' as a scalar when it has %i components and %i copies. "
               " Please use a Bucket operator() that takes appropriate index arguments.\n", file, line, m_fieldName,
               m_numComponents, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_single_component_access(const char* file, int line) const {
    if (m_numComponents != 1) {
      if (line == -1) {
        printf("Error: Accessing BucketValues for Field '%s' as if it only has one component when it actually"
               " has %i components.  Please use a Bucket operator() that also has a component argument.\n",
               m_fieldName, m_numComponents);
      }
      else {
        printf("Error: %s:%i: Accessing BucketValues for Field '%s' as if it only has one component when it actually"
               " has %i components.  Please use a Bucket operator() that also has a component argument.\n",
               file, line, m_fieldName, m_numComponents);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_single_copy_access(const char* file, int line) const {
    if (m_numCopies != 1) {
      if (line == -1) {
        printf("Error: Accessing BucketValues for Field '%s' as if it only has one copy when it actually"
               " has %i copies.  Please use a Bucket operator() that also has a copy argument.\n",
               m_fieldName, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Accessing BucketValues for Field '%s' as if it only has one copy when it actually"
               " has %i copies.  Please use a Bucket operator() that also has a copy argument.\n",
               file, line, m_fieldName, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field access error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_component_bounds(int component, const char* file, int line) const {
    if ((component < 0) || (component >= m_numComponents)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to BucketValues for Field '%s' with component index %i for a"
               " Bucket with %i components.\n", m_fieldName, component, m_numComponents);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to BucketValues for Field '%s' with component index %i for a"
               " Bucket with %i components.\n", file, line, m_fieldName, component, m_numComponents);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_copy_bounds(int copy, const char* file, int line) const {
    if ((copy < 0) || (copy >= m_numCopies)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to BucketValues for Field '%s' with copy index %i for a Bucket"
               " with %i copies.\n", m_fieldName, copy, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to BucketValues for Field '%s' with copy index %i for a Bucket"
               " with %i copies.\n", file, line, m_fieldName, copy, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_copy_and_component_bounds(int copy, int component, const char* file,
                                                              int line) const {
    if (((copy < 0) || (copy >= m_numCopies)) || ((component < 0) || (component >= m_numComponents))) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to BucketValues for Field '%s' with component index %i and copy"
               " index %i for a Bucket with %i components and %i copies.\n", m_fieldName, component, copy,
               m_numComponents, m_numCopies);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to BucketValues for Field '%s' with component index %i and copy"
               " index %i for a Bucket with %i components and %i copies.\n", file, line, m_fieldName, component, copy,
               m_numComponents, m_numCopies);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_scalar_bounds(int scalar, const char* file, int line) const {
    if ((scalar < 0) || (scalar >= m_numCopies*m_numComponents)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to BucketValues for Field '%s' with scalar index %i for a"
               " Bucket with %i scalars.\n", m_fieldName, scalar, m_numCopies*m_numComponents);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to BucketValues for Field '%s' with scalar index %i for a"
               " Bucket with %i scalars.\n", file, line, m_fieldName, scalar, m_numCopies*m_numComponents);
      }
      STK_NGP_ThrowErrorMsg("Field bounds error.");
    }
  }
  KOKKOS_INLINE_FUNCTION void check_entity_bounds(int entity, const char* file, int line) const {
    if ((entity < 0) || (entity >= m_numEntities)) {
      if (line == -1) {
        printf("Error: Out-of-bounds access to BucketValues for Field '%s' with Entity index %i for a Bucket"
               " with %i Entities.\n", m_fieldName, entity, m_numEntities);
      }
      else {
        printf("Error: %s:%i: Out-of-bounds access to BucketValues for Field '%s' with Entity index %i for a Bucket"
               " with %i Entities.\n", file, line, m_fieldName, entity, m_numEntities);
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
  KOKKOS_INLINE_FUNCTION void check_entity_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.  Also prevent direct passing of an int as the
  // EntityIdx argument because it allows casting from an int to support Kokkos indexing.
  KOKKOS_INLINE_FUNCTION T& operator()(int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, int, int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, ComponentIdx) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, CopyIdx) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, CopyIdx, ComponentIdx) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(int, ScalarIdx) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx, int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx, int, int) const = delete;
  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx, ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_numEntities;
  int m_scalarStride;
};


//==============================================================================
// Host BucketValues: Layout::Right
//==============================================================================

template<typename T>
class BucketValues<T, stk::ngp::HostSpace, Layout::Right>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Right;

  inline BucketValues(T* dataPtr, int numComponents, int numCopies, int numEntities,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_numEntities(numEntities)
  {}

  ~BucketValues() = default;

  // The functions below may be used to iterate through your Bucket data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::EntityIdx entity(0); entity < bucketValues.num_entities(); ++entity) {
  //     for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
  //       for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::EntityIdx entity : bucketValues.entities()) {
  //     for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
  //       for (stk::mesh::ComponentIdx component : bucketValues.components()) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may also reuse the Entity
  // index between multiple Fields if they are from the same Bucket.  You may use hard-coded
  // numerical indexing through user-defined literals, such as: 0_entity, 1_entity, 2_entity, etc.,
  // 0_copy, 1_copy, 2_copy, etc., and 0_comp, 1_comp, 2_comp, etc.

  inline int num_entities() const { return m_numEntities; }
  inline EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  inline int num_copies() const { return m_numCopies; }
  inline CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  inline int num_components() const { return m_numComponents; }
  inline ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  inline int num_scalars() const { return m_numCopies*m_numComponents; }
  inline ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  inline bool is_field_defined() const { return m_numComponents != 0; }

  inline T& operator()(EntityIdx entity,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);
    check_entity_bounds(entity, file, line);

    return m_dataPtr[entity];
  }

  inline T& operator()(EntityIdx entity, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_entity_bounds(entity, file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[entity()*m_numComponents + component()];
  }

  inline T& operator()(EntityIdx entity, CopyIdx copy,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[entity()*m_numCopies + copy()];
  }

  inline T& operator()(EntityIdx entity, CopyIdx copy, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[(entity()*m_numCopies + copy())*m_numComponents + component()];
  }

  inline T& operator()(EntityIdx entity, ScalarIdx scalar,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[entity()*m_numCopies*m_numComponents + scalar()];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Bucket data on either the host
  // or device if you follow this pattern:
  //
  //   double* bucketPtr = bucketValues.pointer();
  //   const int entityStride = bucketValues.entity_stride();
  //   const int copyStride = bucketValues.copy_stride();
  //   const int componentStride = bucketValues.component_stride();
  //
  //   for (int entity = 0; entity < bucketValues.num_entities(); ++entity) {
  //     for (int copy = 0; copy < bucketValues.num_copies(); ++copy) {
  //       for (int component = 0; component < bucketValues.num_components(); ++component) {
  //         entityPtr[entity*entityStride + copy*copyStride + component*componentStride] = 0.0;
  //       }
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  inline T* pointer() const { return m_dataPtr; }

  inline int entity_stride() const { return m_numCopies*m_numComponents; }

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
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " that is not defined on this Bucket.");
  }
  inline void check_single_scalar_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents*m_numCopies == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as a scalar when it has " << m_numComponents << " components and " << m_numCopies <<
                        " copies.  Please use a Bucket operator() that takes appropriate index arguments.");
  }
  inline void check_single_component_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as if it only has one component when it actually has " << m_numComponents <<
                        " components.  Please use a Bucket operator() that also has a component argument.");
  }
  inline void check_single_copy_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numCopies == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as if it only has one copy when it actually has " << m_numCopies <<
                        " copies.  Please use a Bucket operator() that also has a copy argument.");
  }
  inline void check_component_bounds(int component, const char* file, int line) const {
    STK_ThrowRequireMsg((component >= 0) && (component < m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with component index " << component << " for a Bucket with " <<
                        m_numComponents << " components.");
  }
  inline void check_copy_bounds(int copy, const char* file, int line) const {
    STK_ThrowRequireMsg((copy >= 0) && (copy < m_numCopies),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with copy index " << copy << " for a Bucket with " << m_numCopies <<
                        " copies.");
  }
  inline void check_copy_and_component_bounds(int copy, int component, const char* file, int line) const {
    STK_ThrowRequireMsg(((copy >= 0) && (copy < m_numCopies)) && ((component >= 0) && (component < m_numComponents)),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with component index " << component << " and copy index " << copy <<
                        " for a Bucket with " << m_numComponents << " components and " << m_numCopies << " copies.");
  }
  inline void check_scalar_bounds(int scalar, const char* file, int line) const {
    STK_ThrowRequireMsg((scalar >= 0) && (scalar < m_numCopies*m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with scalar index " << scalar << " for a Bucket with " <<
                        m_numCopies*m_numComponents << " scalars.");
  }
  inline void check_entity_bounds(int entity, const char* file, int line) const {
    STK_ThrowRequireMsg((entity >= 0) && (entity < m_numEntities),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with Entity index " << entity << " for a Bucket with " << m_numEntities <<
                        " Entities.");
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
  inline void check_entity_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.  Also prevent direct passing of an int as the
  // EntityIdx argument because it allows casting from an int to support Kokkos indexing.
  inline T& operator()(int) const = delete;
  inline T& operator()(int, int) const = delete;
  inline T& operator()(int, int, int) const = delete;
  inline T& operator()(int, ComponentIdx) const = delete;
  inline T& operator()(int, CopyIdx) const = delete;
  inline T& operator()(int, CopyIdx, ComponentIdx) const = delete;
  inline T& operator()(int, ScalarIdx) const = delete;
  inline T& operator()(EntityIdx, int) const = delete;
  inline T& operator()(EntityIdx, int, int) const = delete;
  inline T& operator()(EntityIdx, ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_numEntities;
};


//==============================================================================
// Host BucketValues: Layout::Left
//==============================================================================

template<typename T>
class BucketValues<T, stk::ngp::HostSpace, Layout::Left>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Left;

  inline BucketValues(T* dataPtr, int numComponents, int numCopies, int numEntities, int scalarStride,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_numEntities(numEntities),
      m_scalarStride(scalarStride)
  {}

  ~BucketValues() = default;

  // The functions below may be used to iterate through your Bucket data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::EntityIdx entity(0); entity < bucketValues.num_entities(); ++entity) {
  //     for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
  //       for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::EntityIdx entity : bucketValues.entities()) {
  //     for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
  //       for (stk::mesh::ComponentIdx component : bucketValues.components()) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may also reuse the Entity
  // index between multiple Fields if they are from the same Bucket.  You may use hard-coded
  // numerical indexing through user-defined literals, such as: 0_entity, 1_entity, 2_entity, etc.,
  // 0_copy, 1_copy, 2_copy, etc., and 0_comp, 1_comp, 2_comp, etc.

  inline int num_entities() const { return m_numEntities; }
  inline EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  inline int num_copies() const { return m_numCopies; }
  inline CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  inline int num_components() const { return m_numComponents; }
  inline ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  inline int num_scalars() const { return m_numCopies*m_numComponents; }
  inline ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  inline bool is_field_defined() const { return m_numComponents != 0; }

  inline T& operator()(EntityIdx entity,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);
    check_entity_bounds(entity, file, line);

    return m_dataPtr[entity];
  }

  inline T& operator()(EntityIdx entity, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_entity_bounds(entity, file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[component()*m_scalarStride + entity()];
  }

  inline T& operator()(EntityIdx entity, CopyIdx copy,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[copy()*m_scalarStride + entity()];
  }

  inline T& operator()(EntityIdx entity, CopyIdx copy, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[(copy()*m_numComponents + component())*m_scalarStride + entity()];
  }

  inline T& operator()(EntityIdx entity, ScalarIdx scalar,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[scalar()*m_scalarStride + entity()];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Bucket data on either the host
  // or device if you follow this pattern:
  //
  //   double* bucketPtr = bucketValues.pointer();
  //   const int entityStride = bucketValues.entity_stride();
  //   const int copyStride = bucketValues.copy_stride();
  //   const int componentStride = bucketValues.component_stride();
  //
  //   for (int entity = 0; entity < bucketValues.num_entities(); ++entity) {
  //     for (int copy = 0; copy < bucketValues.num_copies(); ++copy) {
  //       for (int component = 0; component < bucketValues.num_components(); ++component) {
  //         entityPtr[entity*entityStride + copy*copyStride + component*componentStride] = 0.0;
  //       }
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  inline T* pointer() const { return m_dataPtr; }

  inline int entity_stride() const { return 1; }

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
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " that is not defined on this Bucket.");
  }
  inline void check_single_scalar_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents*m_numCopies == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as a scalar when it has " << m_numComponents << " components and " << m_numCopies <<
                        " copies.  Please use a Bucket operator() that takes appropriate index arguments.");
  }
  inline void check_single_component_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as if it only has one component when it actually has " << m_numComponents <<
                        " components.  Please use a Bucket operator() that also has a component argument.");
  }
  inline void check_single_copy_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numCopies == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as if it only has one copy when it actually has " << m_numCopies <<
                        " copies.  Please use a Bucket operator() that also has a copy argument.");
  }
  inline void check_component_bounds(int component, const char* file, int line) const {
    STK_ThrowRequireMsg((component >= 0) && (component < m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with component index " << component << " for a Bucket with " <<
                        m_numComponents << " components.");
  }
  inline void check_copy_bounds(int copy, const char* file, int line) const {
    STK_ThrowRequireMsg((copy >= 0) && (copy < m_numCopies),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with copy index " << copy << " for a Bucket with " << m_numCopies <<
                        " copies.");
  }
  inline void check_copy_and_component_bounds(int copy, int component, const char* file, int line) const {
    STK_ThrowRequireMsg(((copy >= 0) && (copy < m_numCopies)) && ((component >= 0) && (component < m_numComponents)),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with component index " << component << " and copy index " << copy <<
                        " for a Bucket with " << m_numComponents << " components and " << m_numCopies << " copies.");
  }
  inline void check_scalar_bounds(int scalar, const char* file, int line) const {
    STK_ThrowRequireMsg((scalar >= 0) && (scalar < m_numCopies*m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with scalar index " << scalar << " for a Bucket with " <<
                        m_numCopies*m_numComponents << " scalars.");
  }
  inline void check_entity_bounds(int entity, const char* file, int line) const {
    STK_ThrowRequireMsg((entity >= 0) && (entity < m_numEntities),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with Entity index " << entity << " for a Bucket with " << m_numEntities <<
                        " Entities.");
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
  inline void check_entity_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.  Also prevent direct passing of an int as the
  // EntityIdx argument because it allows casting from an int to support Kokkos indexing.
  inline T& operator()(int) const = delete;
  inline T& operator()(int, int) const = delete;
  inline T& operator()(int, int, int) const = delete;
  inline T& operator()(int, ComponentIdx) const = delete;
  inline T& operator()(int, CopyIdx) const = delete;
  inline T& operator()(int, CopyIdx, ComponentIdx) const = delete;
  inline T& operator()(int, ScalarIdx) const = delete;
  inline T& operator()(EntityIdx, int) const = delete;
  inline T& operator()(EntityIdx, int, int) const = delete;
  inline T& operator()(EntityIdx, ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_numEntities;
  int m_scalarStride;
};


//==============================================================================
// Host BucketValues: Layout::Auto
//==============================================================================

template<typename T>
class BucketValues<T, stk::ngp::HostSpace, Layout::Auto>
{
public:
  using value_type = T;
  using space = stk::ngp::HostSpace;
  using mem_space = stk::ngp::HostSpace::mem_space;
  using exec_space = stk::ngp::HostSpace::exec_space;
  static constexpr Layout layout = Layout::Auto;

  inline BucketValues(T* dataPtr, int numComponents, int numCopies, int numEntities,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_numEntities(numEntities),
      m_scalarStride(1),
      m_entityStride(numCopies*numComponents)
  {}

  inline BucketValues(T* dataPtr, int numComponents, int numCopies, int numEntities, int scalarStride,
                      [[maybe_unused]] const char* fieldName)
    : m_dataPtr(dataPtr),
    #ifdef STK_FIELD_BOUNDS_CHECK
      m_fieldName(fieldName),
    #endif
      m_numComponents(numComponents),
      m_numCopies(numCopies),
      m_numEntities(numEntities),
      m_scalarStride(scalarStride),
      m_entityStride(1)
  {}

  ~BucketValues() = default;

  // The functions below may be used to iterate through your Bucket data on either the host or
  // device with strongly-typed indexing.  You may use either a traditional for-loop:
  //
  //   for (stk::mesh::EntityIdx entity(0); entity < bucketValues.num_entities(); ++entity) {
  //     for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
  //       for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // or a range-based for-loop:
  //
  //   for (stk::mesh::EntityIdx entity : bucketValues.entities()) {
  //     for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
  //       for (stk::mesh::ComponentIdx component : bucketValues.components()) {
  //         entityValues(entity, copy, component) = 0.0;
  //       }
  //     }
  //   }
  //
  // You may safely reuse copy and component indices between multiple Fields in the same loop
  // as long as they all have the same extents in each dimension.  You may also reuse the Entity
  // index between multiple Fields if they are from the same Bucket.  You may use hard-coded
  // numerical indexing through user-defined literals, such as: 0_entity, 1_entity, 2_entity, etc.,
  // 0_copy, 1_copy, 2_copy, etc., and 0_comp, 1_comp, 2_comp, etc.

  inline int num_entities() const { return m_numEntities; }
  inline EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  inline int num_copies() const { return m_numCopies; }
  inline CopyIdxProxy copies() const { return CopyIdxProxy(m_numCopies); }

  inline int num_components() const { return m_numComponents; }
  inline ComponentIdxProxy components() const { return ComponentIdxProxy(m_numComponents); }

  inline int num_scalars() const { return m_numCopies*m_numComponents; }
  inline ScalarIdxProxy scalars() const { return ScalarIdxProxy(m_numCopies*m_numComponents); }

  inline bool is_field_defined() const { return m_numComponents != 0; }

  inline T& operator()(EntityIdx entity,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_scalar_access(file, line);
    check_entity_bounds(entity, file, line);

    return m_dataPtr[entity];
  }

  inline T& operator()(EntityIdx entity, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_copy_access(file, line);
    check_entity_bounds(entity, file, line);
    check_component_bounds(component, file, line);

    return m_dataPtr[entity()*m_entityStride + component()*m_scalarStride];
  }

  inline T& operator()(EntityIdx entity, CopyIdx copy,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_single_component_access(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_bounds(copy, file, line);

    return m_dataPtr[entity()*m_entityStride + copy()*m_scalarStride];
  }

  inline T& operator()(EntityIdx entity, CopyIdx copy, ComponentIdx component,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_copy_and_component_bounds(copy, component, file, line);

    return m_dataPtr[entity()*m_entityStride + (copy()*m_numComponents + component())*m_scalarStride];
  }

  inline T& operator()(EntityIdx entity, ScalarIdx scalar,
                       const char* file = STK_HOST_FILE, int line = STK_HOST_LINE) const
  {
    check_defined_field(file, line);
    check_entity_bounds(entity, file, line);
    check_scalar_bounds(scalar, file, line);

    return m_dataPtr[entity()*m_entityStride + scalar()*m_scalarStride];
  }


  // The functions below are intended for expert users only, as they bypass all correctness and
  // consistency safeguards and perform no bounds checking.  Use these if you must have raw pointer
  // access to the data.  You can use them to iterate through the Bucket data on either the host
  // or device if you follow this pattern:
  //
  //   double* bucketPtr = bucketValues.pointer();
  //   const int entityStride = bucketValues.entity_stride();
  //   const int copyStride = bucketValues.copy_stride();
  //   const int componentStride = bucketValues.component_stride();
  //
  //   for (int entity = 0; entity < bucketValues.num_entities(); ++entity) {
  //     for (int copy = 0; copy < bucketValues.num_copies(); ++copy) {
  //       for (int component = 0; component < bucketValues.num_components(); ++component) {
  //         entityPtr[entity*entityStride + copy*copyStride + component*componentStride] = 0.0;
  //       }
  //     }
  //   }
  //
  // You may leave off the copy indexing if you are certain you have a single-copy Field and you
  // may leave off the component indexing if you are certain that you have a single-component Field.

  inline T* pointer() const { return m_dataPtr; }

  inline int entity_stride() const { return m_entityStride; }

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
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " that is not defined on this Bucket.");
  }
  inline void check_single_scalar_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents*m_numCopies == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as a scalar when it has " << m_numComponents << " components and " << m_numCopies <<
                        " copies.  Please use a Bucket operator() that takes appropriate index arguments.");
  }
  inline void check_single_component_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numComponents == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as if it only has one component when it actually has " << m_numComponents <<
                        " components.  Please use a Bucket operator() that also has a component argument.");
  }
  inline void check_single_copy_access(const char* file, int line) const {
    STK_ThrowRequireMsg(m_numCopies == 1,
                        location_string(file, line) << "Accessing BucketValues for Field '" << m_fieldName << "'"
                        " as if it only has one copy when it actually has " << m_numCopies <<
                        " copies.  Please use a Bucket operator() that also has a copy argument.");
  }
  inline void check_component_bounds(int component, const char* file, int line) const {
    STK_ThrowRequireMsg((component >= 0) && (component < m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with component index " << component << " for a Bucket with " <<
                        m_numComponents << " components.");
  }
  inline void check_copy_bounds(int copy, const char* file, int line) const {
    STK_ThrowRequireMsg((copy >= 0) && (copy < m_numCopies),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with copy index " << copy << " for a Bucket with " << m_numCopies <<
                        " copies.");
  }
  inline void check_copy_and_component_bounds(int copy, int component, const char* file, int line) const {
    STK_ThrowRequireMsg(((copy >= 0) && (copy < m_numCopies)) && ((component >= 0) && (component < m_numComponents)),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with component index " << component << " and copy index " << copy <<
                        " for a Bucket with " << m_numComponents << " components and " << m_numCopies << " copies.");
  }
  inline void check_scalar_bounds(int scalar, const char* file, int line) const {
    STK_ThrowRequireMsg((scalar >= 0) && (scalar < m_numCopies*m_numComponents),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with scalar index " << scalar << " for a Bucket with " <<
                        m_numCopies*m_numComponents << " scalars.");
  }
  inline void check_entity_bounds(int entity, const char* file, int line) const {
    STK_ThrowRequireMsg((entity >= 0) && (entity < m_numEntities),
                        location_string(file, line) << "Out-of-bounds access to BucketValues for Field '" <<
                        m_fieldName << "' with Entity index " << entity << " for a Bucket with " << m_numEntities <<
                        " Entities.");
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
  inline void check_entity_bounds(int, const char*, int) const {}
#endif

  // Disallow some improper overload resolutions where a 0 mistakenly passed as a strongly-typed
  // index argument is implicitly cast to a "const char*" for the internal-use debugging output
  // argument of a lower-dimension overload.  Also prevent direct passing of an int as the
  // EntityIdx argument because it allows casting from an int to support Kokkos indexing.
  inline T& operator()(int) const = delete;
  inline T& operator()(int, int) const = delete;
  inline T& operator()(int, int, int) const = delete;
  inline T& operator()(int, ComponentIdx) const = delete;
  inline T& operator()(int, CopyIdx) const = delete;
  inline T& operator()(int, CopyIdx, ComponentIdx) const = delete;
  inline T& operator()(int, ScalarIdx) const = delete;
  inline T& operator()(EntityIdx, int) const = delete;
  inline T& operator()(EntityIdx, int, int) const = delete;
  inline T& operator()(EntityIdx, ComponentIdx, int) const = delete;

  T* m_dataPtr;
#ifdef STK_FIELD_BOUNDS_CHECK
  const char* m_fieldName;
#endif
  int m_numComponents;
  int m_numCopies;
  int m_numEntities;
  int m_scalarStride;
  int m_entityStride;
};

}

#endif // STK_BUCKETVALUES_HPP
