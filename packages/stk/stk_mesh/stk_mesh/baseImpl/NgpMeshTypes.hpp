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

#ifndef STK_MESH_NGPMESHTYPES_HPP
#define STK_MESH_NGPMESHTYPES_HPP

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "Kokkos_Core.hpp"
#include <limits>
#include <ostream>
#include <vector>

namespace stk::mesh::impl {

constexpr unsigned INVALID_INDEX = std::numeric_limits<unsigned>::max();

struct DevicePartOrdinalLess
{
  KOKKOS_DEFAULTED_FUNCTION
  DevicePartOrdinalLess() = default;

  template <typename PartOrdinal>
  KOKKOS_INLINE_FUNCTION
  bool operator()(PartOrdinal const& lhs, PartOrdinal const& rhs) {
    if (lhs.extent(0) != rhs.extent(0)) {
      return lhs.extent(0) < rhs.extent(0);
    }

    for (unsigned i = 0; i < lhs.extent(0); ++i) {
      if (lhs(i) != rhs(i)) {
        return lhs(i) < rhs(i);
      }
    }
    return false;
  }

  bool operator()(const std::vector<stk::mesh::PartOrdinal>& lhs, const std::vector<stk::mesh::PartOrdinal>& rhs)
  {
    using VectorWrapper = Kokkos::View<const stk::mesh::PartOrdinal*, stk::ngp::HostSpace::mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    VectorWrapper lhsWrapper(lhs.data(), lhs.size());
    VectorWrapper rhsWrapper(rhs.data(), rhs.size());

    return (*this)(lhsWrapper, rhsWrapper);
  }
};

template <typename T = Entity>
struct EntityWrapper
{
  KOKKOS_FUNCTION
  EntityWrapper()
    : value(T{}),
      isUserInputEntity(false),
      isForPartInduction(false)
  {}

  KOKKOS_FUNCTION
  EntityWrapper(T e, bool u = false, bool p = false)
    : value(e),
      isUserInputEntity(u),
      isForPartInduction(p)
  {}

  KOKKOS_INLINE_FUNCTION
  operator T() const {
    return value;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(EntityWrapper other) const {
    return (value == other.value);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator<(EntityWrapper other) const {
    if (value == other.value) {
      return isForPartInduction == false;
    } else {
      return value < other.value;
    }
  }

  T value;
  bool isUserInputEntity;
  bool isForPartInduction;
};

struct DeniedPartRemoval
{
  Entity entity;
  PartOrdinal part;

  KOKKOS_INLINE_FUNCTION
  bool operator<(DeniedPartRemoval const& rhs) const {
    if (entity.local_offset() != rhs.entity.local_offset()) {
      return entity.local_offset() < rhs.entity.local_offset();
    }
    return part < rhs.part;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(DeniedPartRemoval const& rhs) const {
    return (entity == rhs.entity) && (part == rhs.part);
  }
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const EntityWrapper<T>& entityWrapper)
{
  os << static_cast<T>(entityWrapper);
  return os;
}

template <typename WrappedEntityViewType, typename PartOrdinalsViewType, typename PartOrdinalsProxyViewType>
struct NewPartListsPerEntity
{
  WrappedEntityViewType wrappedEntities;
  PartOrdinalsViewType newPartOrdinalsPerEntity;
  PartOrdinalsProxyViewType partOrdinalsProxy;
  unsigned maxNewNumPartsPerEntity = 0;
};

struct PartOrdinalsProxyIndices {
  KOKKOS_DEFAULTED_FUNCTION
  PartOrdinalsProxyIndices() = default;

  KOKKOS_INLINE_FUNCTION
  PartOrdinalsProxyIndices(EntityRank entityRank, Ordinal* ptr, unsigned len)
    : rank(entityRank),
      startPtr(ptr),
      length(len)
  {}

  KOKKOS_INLINE_FUNCTION
  bool operator==(PartOrdinalsProxyIndices const& other) const
  {
    if (rank != other.rank) { return false; }
    if (length != other.length) { return false; }

    for (unsigned i = 0; i < length; ++i) {
      auto thisOrd = startPtr + i;
      auto otherOrd = other.startPtr + i;
      if (*thisOrd != *otherOrd) { return false; }
    }
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator<(PartOrdinalsProxyIndices const& other) const
  {
    if (rank != other.rank) { return rank < other.rank; }
    if (length != other.length) { return length < other.length; }

    for (unsigned i = 0; i < length; ++i) {
      auto thisOrd = startPtr + i;
      auto otherOrd = other.startPtr + i;
      if (*thisOrd != *otherOrd) {
        return *thisOrd < *otherOrd;
      }
    }
    return false;
  }

  EntityRank rank;
  Ordinal* startPtr;
  unsigned length;
};

} // namespace stk::mesh::impl

#endif
