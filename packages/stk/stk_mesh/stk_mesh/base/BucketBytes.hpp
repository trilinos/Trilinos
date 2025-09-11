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

#ifndef STK_BUCKETBYTES_HPP
#define STK_BUCKETBYTES_HPP

#include "stk_util/stk_config.h"
#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_mesh/base/FieldIndexTypes.hpp"
#include "Kokkos_Macros.hpp"

namespace stk::mesh {

//==============================================================================
// Device BucketBytes
//==============================================================================

template<typename T = std::byte, typename MemSpace = stk::ngp::HostMemSpace>
class BucketBytes {
public:
  KOKKOS_INLINE_FUNCTION BucketBytes(T* bytePtr, int numBytesPerEntity, int numBytesPerScalar, int numEntities,
                                     int scalarStride)
    : m_bytePtr(bytePtr),
      m_numBytesPerEntity(numBytesPerEntity),
      m_numBytesPerScalar(numBytesPerScalar),
      m_numEntities(numEntities),
      m_scalarByteStride(scalarStride*numBytesPerScalar)
  {}

  KOKKOS_DEFAULTED_FUNCTION ~BucketBytes() = default;

  KOKKOS_INLINE_FUNCTION int num_bytes() const { return m_numBytesPerEntity; }
  KOKKOS_INLINE_FUNCTION ByteIdxProxy bytes() const { return ByteIdxProxy(m_numBytesPerEntity); }

  KOKKOS_INLINE_FUNCTION int num_entities() const { return m_numEntities; }
  KOKKOS_INLINE_FUNCTION EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  KOKKOS_INLINE_FUNCTION bool is_field_defined() const { return m_numBytesPerEntity != 0; }

  KOKKOS_INLINE_FUNCTION T& operator()(EntityIdx entity, ByteIdx byte) const {
    const int scalar = static_cast<int>(byte) / m_numBytesPerScalar;
    const int byteInScalar = static_cast<int>(byte) % m_numBytesPerScalar;
    return m_bytePtr[entity*m_numBytesPerScalar + scalar*m_scalarByteStride + byteInScalar];
  }

private:
  T* m_bytePtr;
  int m_numBytesPerEntity;
  int m_numBytesPerScalar;
  int m_numEntities;
  int m_scalarByteStride;
};


//==============================================================================
// Host BucketBytes
//==============================================================================

template<typename T>
class BucketBytes<T, stk::ngp::HostMemSpace> {
public:

  inline BucketBytes(T* bytePtr, int numBytesPerEntity, int numBytesPerScalar, int numEntities, int scalarStride)
    : m_bytePtr(bytePtr),
      m_numBytesPerEntity(numBytesPerEntity),
      m_numEntities(numEntities),
      m_numBytesPerScalar(numBytesPerScalar),
      m_scalarByteStride(scalarStride*numBytesPerScalar),
      m_isLayoutRight(false)
  {}

  inline BucketBytes(T* bytePtr, int numBytesPerEntity, int numEntities)
    : m_bytePtr(bytePtr),
      m_numBytesPerEntity(numBytesPerEntity),
      m_numEntities(numEntities),
      m_isLayoutRight(true)
  {}

  ~BucketBytes() = default;

  inline int num_bytes() const { return m_numBytesPerEntity; }
  inline ByteIdxProxy bytes() const { return ByteIdxProxy(m_numBytesPerEntity); }

  inline int num_entities() const { return m_numEntities; }
  inline EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  inline bool is_field_defined() const { return m_numBytesPerEntity != 0; }

  inline T& operator()(EntityIdx entity, ByteIdx byte) const {
    if (m_isLayoutRight) {
      return m_bytePtr[static_cast<int>(entity)*m_numBytesPerEntity + static_cast<int>(byte)];
    }
    else {
      const int scalar = static_cast<int>(byte) / m_numBytesPerScalar;
      const int byteInScalar = static_cast<int>(byte) % m_numBytesPerScalar;
      return m_bytePtr[entity*m_numBytesPerScalar + scalar*m_scalarByteStride + byteInScalar];
    }
  }

private:
  T* m_bytePtr;
  int m_numBytesPerEntity;
  int m_numEntities;
  int m_numBytesPerScalar;
  int m_scalarByteStride;
  bool m_isLayoutRight;
};


//==============================================================================
// Host BucketBytesLeft
//==============================================================================

template<typename T>
class BucketBytesLeft {
public:

  inline BucketBytesLeft(T* bytePtr, int numBytesPerEntity, int numBytesPerScalar, int numEntities, int scalarStride)
    : m_bytePtr(bytePtr),
      m_numBytesPerEntity(numBytesPerEntity),
      m_numEntities(numEntities),
      m_numBytesPerScalar(numBytesPerScalar),
      m_scalarByteStride(scalarStride*numBytesPerScalar)
  {}

  ~BucketBytesLeft() = default;

  inline int num_bytes() const { return m_numBytesPerEntity; }
  inline ByteIdxProxy bytes() const { return ByteIdxProxy(m_numBytesPerEntity); }

  inline int num_entities() const { return m_numEntities; }
  inline EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  inline bool is_field_defined() const { return m_numBytesPerEntity != 0; }

  inline T& operator()(EntityIdx entity, ByteIdx byte) const {
    const int scalar = static_cast<int>(byte) / m_numBytesPerScalar;
    const int byteInScalar = static_cast<int>(byte) % m_numBytesPerScalar;
    return m_bytePtr[entity*m_numBytesPerScalar + scalar*m_scalarByteStride + byteInScalar];
  }

private:
  T* m_bytePtr;
  int m_numBytesPerEntity;
  int m_numEntities;
  int m_numBytesPerScalar;
  int m_scalarByteStride;
};


//==============================================================================
// Host BucketBytesRight
//==============================================================================

template<typename T>
class BucketBytesRight {
public:

  inline BucketBytesRight(T* bytePtr, int numBytesPerEntity, int numEntities)
    : m_bytePtr(bytePtr),
      m_numBytesPerEntity(numBytesPerEntity),
      m_numEntities(numEntities)
  {}

  ~BucketBytesRight() = default;

  inline int num_bytes() const { return m_numBytesPerEntity; }
  inline ByteIdxProxy bytes() const { return ByteIdxProxy(m_numBytesPerEntity); }

  inline int num_entities() const { return m_numEntities; }
  inline EntityIdxProxy entities() const { return EntityIdxProxy(m_numEntities); }

  inline bool is_field_defined() const { return m_numBytesPerEntity != 0; }

  inline T& operator()(EntityIdx entity, ByteIdx byte) const {
    return m_bytePtr[static_cast<int>(entity)*m_numBytesPerEntity + static_cast<int>(byte)];
  }

private:
  T* m_bytePtr;
  int m_numBytesPerEntity;
  int m_numEntities;
};

}

#endif // STK_BUCKETBYTES_HPP
