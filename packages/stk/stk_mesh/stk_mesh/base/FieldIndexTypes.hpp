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
//

#ifndef FIELDINDEXTYPES_HPP
#define FIELDINDEXTYPES_HPP

#include "Kokkos_Macros.hpp"
#include <type_traits>
#include <iostream>  // TODO: DEBUG

namespace stk::mesh {

// Kokkos likes to increment/decrement index variables with both itegral types and unnamed enum values
template <typename T> using EnableIfIntegral = std::enable_if_t<std::is_integral_v<T> || std::is_enum_v<T>>;

//==============================================================================
// Base class for all of the strongly-typed Field index types.  This is meant
// to interoperate cleanly with other of the same index types as well as with
// plain integral types.
//
template <typename Derived, typename IndexType>
class FieldIndex {
public:
  using index_type = IndexType;

  // Implicit conversion constructor from integral types
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION FieldIndex(T value) noexcept : m_value(static_cast<IndexType>(value)) {}
  KOKKOS_DEFAULTED_FUNCTION ~FieldIndex() = default;

  constexpr KOKKOS_DEFAULTED_FUNCTION FieldIndex(const FieldIndex&) noexcept = default;
  constexpr KOKKOS_DEFAULTED_FUNCTION FieldIndex(FieldIndex&&) noexcept = default;
  constexpr KOKKOS_DEFAULTED_FUNCTION FieldIndex& operator=(const FieldIndex&) noexcept = default;
  constexpr KOKKOS_DEFAULTED_FUNCTION FieldIndex& operator=(FieldIndex&&) noexcept = default;

  // Direct accessor for the internal integral value
  constexpr KOKKOS_INLINE_FUNCTION index_type operator()() const noexcept { return m_value; }

  // Implicit conversion operators to integral types
  constexpr KOKKOS_INLINE_FUNCTION operator int() const noexcept { return static_cast<int>(m_value); }
  constexpr KOKKOS_INLINE_FUNCTION operator unsigned int() const noexcept { return static_cast<unsigned int>(m_value); }
  constexpr KOKKOS_INLINE_FUNCTION operator long int() const noexcept { return static_cast<long int>(m_value); }
  constexpr KOKKOS_INLINE_FUNCTION operator unsigned long int() const noexcept { return static_cast<unsigned long int>(m_value); }

  // Mutating operators
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator++() { ++m_value; return static_cast<Derived&>(*this); }
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator--() { --m_value; return static_cast<Derived&>(*this); }

  constexpr KOKKOS_INLINE_FUNCTION Derived operator++(int) {
    Derived temp = static_cast<Derived&>(*this); ++m_value; return temp;
  }
  constexpr KOKKOS_INLINE_FUNCTION Derived operator--(int) {
    Derived temp = static_cast<Derived&>(*this); --m_value; return temp;
  }

  constexpr KOKKOS_INLINE_FUNCTION Derived& operator+=(const Derived& rhs) {
    m_value += rhs.m_value; return static_cast<Derived&>(*this);
  }
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator-=(const Derived& rhs) {
    m_value -= rhs.m_value; return static_cast<Derived&>(*this);
  }

  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator+=(T rhs) { m_value += rhs; return static_cast<Derived&>(*this); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator-=(T rhs) { m_value -= rhs; return static_cast<Derived&>(*this); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator*=(T rhs) { m_value *= rhs; return static_cast<Derived&>(*this); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator/=(T rhs) { m_value /= rhs; return static_cast<Derived&>(*this); }

  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived& operator=(T rhs) { m_value = rhs; return static_cast<Derived&>(*this); }

  // Arithmetic operators
  constexpr KOKKOS_INLINE_FUNCTION Derived operator+(const Derived& other) const { return Derived(m_value + other.m_value); }
  constexpr KOKKOS_INLINE_FUNCTION Derived operator-(const Derived& other) const { return Derived(m_value - other.m_value); }
  constexpr KOKKOS_INLINE_FUNCTION Derived operator*(const Derived& other) const { return Derived(m_value * other.m_value); }
  constexpr KOKKOS_INLINE_FUNCTION Derived operator/(const Derived& other) const { return Derived(m_value / other.m_value); }

  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived operator+(T rhs) const { return Derived(m_value + rhs); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived operator-(T rhs) const { return Derived(m_value - rhs); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived operator*(T rhs) const { return Derived(m_value * rhs); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION Derived operator/(T rhs) const { return Derived(m_value / rhs); }

  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend Derived operator+(T lhs, const Derived& rhs) { return Derived(lhs + rhs.m_value); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend Derived operator-(T lhs, const Derived& rhs) { return Derived(lhs - rhs.m_value); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend Derived operator*(T lhs, const Derived& rhs) { return Derived(lhs * rhs.m_value); }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend Derived operator/(T lhs, const Derived& rhs) { return Derived(lhs / rhs.m_value); }

  // Comparison operators
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator==(const Derived& lhs, const Derived& rhs) { return lhs.m_value == rhs.m_value; }
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator!=(const Derived& lhs, const Derived& rhs) { return lhs.m_value != rhs.m_value; }
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator< (const Derived& lhs, const Derived& rhs) { return lhs.m_value <  rhs.m_value; }
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator<=(const Derived& lhs, const Derived& rhs) { return lhs.m_value <= rhs.m_value; }
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator> (const Derived& lhs, const Derived& rhs) { return lhs.m_value >  rhs.m_value; }
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator>=(const Derived& lhs, const Derived& rhs) { return lhs.m_value >= rhs.m_value; }

  // Comparison operators with integral types
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator==(const Derived& lhs, T rhs) { return lhs.m_value == rhs; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator!=(const Derived& lhs, T rhs) { return lhs.m_value != rhs; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator< (const Derived& lhs, T rhs) { return lhs.m_value <  rhs; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator<=(const Derived& lhs, T rhs) { return lhs.m_value <= rhs; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator> (const Derived& lhs, T rhs) { return lhs.m_value >  rhs; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator>=(const Derived& lhs, T rhs) { return lhs.m_value >= rhs; }

  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator==(T lhs, const Derived& rhs) { return lhs == rhs.m_value; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator!=(T lhs, const Derived& rhs) { return lhs != rhs.m_value; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator< (T lhs, const Derived& rhs) { return lhs <  rhs.m_value; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator<=(T lhs, const Derived& rhs) { return lhs <= rhs.m_value; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator> (T lhs, const Derived& rhs) { return lhs >  rhs.m_value; }
  template <typename T, typename = EnableIfIntegral<T>>
  constexpr KOKKOS_INLINE_FUNCTION friend bool operator>=(T lhs, const Derived& rhs) { return lhs >= rhs.m_value; }

private:
  int m_value;
};

template <typename Derived, typename IndexType>
std::ostream& operator<<(std::ostream& os, const FieldIndex<Derived, IndexType>& rhs)
{
  os << static_cast<int>(rhs);
  return os;
}

template <typename T, typename Derived, typename IndexType, typename = EnableIfIntegral<T>>
constexpr KOKKOS_INLINE_FUNCTION T& operator+=(T& lhs, const FieldIndex<Derived, IndexType>& rhs) { return lhs += static_cast<T>(rhs); }

template <typename T, typename Derived, typename IndexType, typename = EnableIfIntegral<T>>
constexpr KOKKOS_INLINE_FUNCTION T& operator-=(T& lhs, const FieldIndex<Derived, IndexType>& rhs) { return lhs -= static_cast<T>(rhs); }

template <typename T, typename Derived, typename IndexType, typename = EnableIfIntegral<T>>
constexpr KOKKOS_INLINE_FUNCTION T& operator*=(T& lhs, const FieldIndex<Derived, IndexType>& rhs) { return lhs *= static_cast<T>(rhs); }

template <typename T, typename Derived, typename IndexType, typename = EnableIfIntegral<T>>
constexpr KOKKOS_INLINE_FUNCTION T& operator/=(T& lhs, const FieldIndex<Derived, IndexType>& rhs) { return lhs /= static_cast<T>(rhs); }

template <typename Derived, typename IndexType>
class FieldIndexIterator {
public:
  using index_type = typename IndexType::index_type;

  constexpr KOKKOS_INLINE_FUNCTION FieldIndexIterator(index_type value) : m_value(value) {}
  constexpr KOKKOS_INLINE_FUNCTION FieldIndexIterator& operator++() { ++m_value; return *this; }
  constexpr KOKKOS_INLINE_FUNCTION FieldIndexIterator& operator--() { --m_value; return *this; }
  constexpr KOKKOS_INLINE_FUNCTION IndexType operator*() const { return IndexType(m_value); }
  constexpr KOKKOS_INLINE_FUNCTION bool operator!=(const FieldIndexIterator& other) const { return m_value != other.m_value; }
private:
  index_type m_value;
};

template <typename Derived, typename IteratorType,
          typename = std::enable_if<std::is_signed_v<typename IteratorType::index_type>>>
class FieldIndexProxy {
public:
  using index_type = typename IteratorType::index_type;

  constexpr KOKKOS_INLINE_FUNCTION FieldIndexProxy(index_type size) : m_size(size) {}
  constexpr KOKKOS_INLINE_FUNCTION IteratorType begin() { return IteratorType(0); }
  constexpr KOKKOS_INLINE_FUNCTION IteratorType end() { return IteratorType(m_size); }
  constexpr KOKKOS_INLINE_FUNCTION IteratorType rbegin() { return IteratorType(m_size-1); }
  constexpr KOKKOS_INLINE_FUNCTION IteratorType rend() { return IteratorType(-1); }
private:
  index_type m_size;
};


//==============================================================================
class ComponentIdx : public FieldIndex<ComponentIdx, int> {
public:
  constexpr KOKKOS_INLINE_FUNCTION explicit ComponentIdx(int value) noexcept : FieldIndex(value) {}
  constexpr KOKKOS_INLINE_FUNCTION ComponentIdx& operator=(int rhs) noexcept { FieldIndex::operator=(rhs); return *this; }
};

class ComponentIdxIterator : public FieldIndexIterator<ComponentIdxIterator, ComponentIdx> {
public:
  constexpr KOKKOS_INLINE_FUNCTION ComponentIdxIterator(ComponentIdx::index_type value) : FieldIndexIterator(value) {}
};

class ComponentIdxProxy : public FieldIndexProxy<ComponentIdxProxy, ComponentIdxIterator> {
public:
  constexpr KOKKOS_INLINE_FUNCTION ComponentIdxProxy(ComponentIdxIterator::index_type size) : FieldIndexProxy(size) {}
};


//==============================================================================
class CopyIdx : public FieldIndex<CopyIdx, int> {
public:
  constexpr KOKKOS_INLINE_FUNCTION explicit CopyIdx(int value) noexcept : FieldIndex(value) {}
  constexpr KOKKOS_INLINE_FUNCTION CopyIdx& operator=(int rhs) noexcept { FieldIndex::operator=(rhs); return *this; }
};

class CopyIdxIterator : public FieldIndexIterator<CopyIdxIterator, CopyIdx> {
public:
  constexpr KOKKOS_INLINE_FUNCTION CopyIdxIterator(CopyIdx::index_type value) : FieldIndexIterator(value) {}
};

class CopyIdxProxy : public FieldIndexProxy<CopyIdxProxy, CopyIdxIterator> {
public:
  constexpr KOKKOS_INLINE_FUNCTION CopyIdxProxy(CopyIdxIterator::index_type size) : FieldIndexProxy(size) {}
};


//==============================================================================
class ScalarIdx : public FieldIndex<ScalarIdx, int> {
public:
  constexpr KOKKOS_INLINE_FUNCTION explicit ScalarIdx(int value) noexcept : FieldIndex(value) {}
  constexpr KOKKOS_INLINE_FUNCTION ScalarIdx& operator=(int rhs) noexcept { FieldIndex::operator=(rhs); return *this; }
};

class ScalarIdxIterator : public FieldIndexIterator<ScalarIdxIterator, ScalarIdx> {
public:
  constexpr KOKKOS_INLINE_FUNCTION ScalarIdxIterator(ScalarIdx::index_type value) : FieldIndexIterator(value) {}
};

class ScalarIdxProxy : public FieldIndexProxy<ScalarIdxProxy, ScalarIdxIterator> {
public:
  constexpr KOKKOS_INLINE_FUNCTION ScalarIdxProxy(ScalarIdxIterator::index_type size) : FieldIndexProxy(size) {}
};


//==============================================================================
// This must be implicitly-convertible from int to interoperate properly with
// Kokkos thread team iteration
class EntityIdx : public FieldIndex<EntityIdx, int> {
public:
  constexpr KOKKOS_INLINE_FUNCTION EntityIdx(int value) noexcept : FieldIndex(value) {}
  constexpr KOKKOS_INLINE_FUNCTION EntityIdx& operator=(int rhs) noexcept { FieldIndex::operator=(rhs); return *this; }
};

class EntityIdxIterator : public FieldIndexIterator<EntityIdxIterator, EntityIdx> {
public:
  constexpr KOKKOS_INLINE_FUNCTION EntityIdxIterator(EntityIdx::index_type value) : FieldIndexIterator(value) {}
};

class EntityIdxProxy : public FieldIndexProxy<EntityIdxProxy, EntityIdxIterator> {
public:
  constexpr KOKKOS_INLINE_FUNCTION EntityIdxProxy(EntityIdxIterator::index_type size) : FieldIndexProxy(size) {}
};

//==============================================================================
class ByteIdx : public FieldIndex<ByteIdx, int> {
public:
  constexpr KOKKOS_INLINE_FUNCTION explicit ByteIdx(int value) noexcept : FieldIndex(value) {}
  constexpr KOKKOS_INLINE_FUNCTION ByteIdx& operator=(int rhs) noexcept { FieldIndex::operator=(rhs); return *this; }
};

class ByteIdxIterator : public FieldIndexIterator<ByteIdxIterator, ByteIdx> {
public:
  constexpr KOKKOS_INLINE_FUNCTION ByteIdxIterator(ByteIdx::index_type value) : FieldIndexIterator(value) {}
};

class ByteIdxProxy : public FieldIndexProxy<ByteIdxProxy, ByteIdxIterator> {
public:
  constexpr KOKKOS_INLINE_FUNCTION ByteIdxProxy(ByteIdxIterator::index_type size) : FieldIndexProxy(size) {}
};

}

//==============================================================================
constexpr KOKKOS_INLINE_FUNCTION stk::mesh::ComponentIdx operator ""_comp(unsigned long long int comp) { return stk::mesh::ComponentIdx(comp); }
constexpr KOKKOS_INLINE_FUNCTION stk::mesh::CopyIdx operator ""_copy(unsigned long long int copy) { return stk::mesh::CopyIdx(copy); }
constexpr KOKKOS_INLINE_FUNCTION stk::mesh::ScalarIdx operator ""_scalar(unsigned long long int scalar) { return stk::mesh::ScalarIdx(scalar); }
constexpr KOKKOS_INLINE_FUNCTION stk::mesh::EntityIdx operator ""_entity(unsigned long long int entity) { return stk::mesh::EntityIdx(entity); }
constexpr KOKKOS_INLINE_FUNCTION stk::mesh::ByteIdx operator ""_byte(unsigned long long int byte) { return stk::mesh::ByteIdx(byte); }

constexpr stk::mesh::ComponentIdx x_comp(0);
constexpr stk::mesh::ComponentIdx y_comp(1);
constexpr stk::mesh::ComponentIdx z_comp(2);

//==============================================================================
#endif // FIELDINDEXTYPES_HPP
