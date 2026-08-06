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


#ifndef STK_MESH_VIEWVECTOR_HPP
#define STK_MESH_VIEWVECTOR_HPP

#include "stk_util/ngp/NgpSpaces.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include <string>
#include <utility>
#include <cstddef>

namespace stk::mesh::impl {

template <typename T, typename MemSpace = stk::ngp::HostMemSpace, typename SizeT = unsigned>
class ViewVector {
public:
  using value_type = T;
  using index_type = SizeT;
  using mem_space = MemSpace;
  using buffer_type = Kokkos::View<T*, MemSpace>;
  using unmanaged_buffer_type = Kokkos::View<T*, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using size_view_type = Kokkos::View<SizeT, stk::ngp::HostPinnedSpace>;
  static constexpr SizeT growth_scale = 2;

  explicit ViewVector(const std::string& name);
  explicit ViewVector(const std::string& name, SizeT size);
  KOKKOS_INLINE_FUNCTION ~ViewVector();

  KOKKOS_DEFAULTED_FUNCTION ViewVector() = default;
  KOKKOS_DEFAULTED_FUNCTION ViewVector(const ViewVector& other) = default;
  KOKKOS_DEFAULTED_FUNCTION ViewVector(ViewVector&& other) = default;
  KOKKOS_DEFAULTED_FUNCTION ViewVector& operator=(const ViewVector& rhs) = default;
  KOKKOS_DEFAULTED_FUNCTION ViewVector& operator=(ViewVector&& rhs) = default;

  KOKKOS_INLINE_FUNCTION SizeT size() const;
  KOKKOS_INLINE_FUNCTION SizeT capacity() const;
  KOKKOS_INLINE_FUNCTION bool empty() const;

  KOKKOS_INLINE_FUNCTION T& operator[](SizeT idx) const;

  KOKKOS_INLINE_FUNCTION T* data() const;

  std::string name() const;
  void reserve(SizeT newCapacity);   // Change capacity to requested value; Increase only
  void resize(SizeT newSize);        // Change size to requested value; Grow capacity to precisely fit if needed
  void resize_scale(SizeT newSize);  // Change size to requested value; Grow capacity by scale factor if needed
  void clear();
  void swap(ViewVector& other);

  void push_back(const T& value);
  void push_back(T&& value);

  template <typename... Args>
  void emplace_back(Args&&... args);

  KOKKOS_INLINE_FUNCTION auto& get_view() const;
  KOKKOS_INLINE_FUNCTION bool is_allocated() const;

private:
  void change_capacity(SizeT newCapacity);

  buffer_type m_buffer;
  size_view_type m_size;
};

//------------------------------------------------------------------------------
template <typename BufferType>
void construct_entries(const BufferType& buffer,
                       typename BufferType::index_type beginIdx,  typename BufferType::index_type endIdx)
{
  using T = typename BufferType::element_type;
  using SizeT = typename BufferType::index_type;
  using ExecSpace = typename BufferType::memory_space::execution_space;

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(beginIdx, endIdx),
    KOKKOS_LAMBDA(const SizeT i) {
      new (&buffer[i]) T{};
    }
  );
}

//------------------------------------------------------------------------------
template <typename BufferType>
requires stk::ngp::is_host_exec_space<typename BufferType::memory_space::execution_space>
void destroy_entries(const BufferType& buffer,
                     typename BufferType::index_type beginIdx, typename BufferType::index_type endIdx)
{
  using T = typename BufferType::element_type;
  using SizeT = typename BufferType::index_type;
  using ExecSpace = typename BufferType::memory_space::execution_space;

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(beginIdx, endIdx),
    [&](const SizeT i) {
      (&buffer[i])->~T();
    }
  );
}

template <typename BufferType>
requires stk::ngp::is_device_exec_space<typename BufferType::memory_space::execution_space>
void destroy_entries(const BufferType& buffer,
                     typename BufferType::index_type beginIdx, typename BufferType::index_type endIdx)
{
  using T = typename BufferType::element_type;
  using SizeT = typename BufferType::index_type;
  using ExecSpace = typename BufferType::memory_space::execution_space;

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(beginIdx, endIdx),
    KOKKOS_LAMBDA(const SizeT i) {
      (&buffer[i])->~T();
    }
  );
}

//------------------------------------------------------------------------------
template <typename BufferType>
void move_entries(const BufferType& srcBuffer, const BufferType& dstBuffer,
                  typename BufferType::index_type beginIdx, typename BufferType::index_type endIdx)
{
  using T = typename BufferType::element_type;
  using SizeT = typename BufferType::index_type;
  using ExecSpace = typename BufferType::memory_space::execution_space;

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(beginIdx, endIdx),
    KOKKOS_LAMBDA(const SizeT i) {
      new (&dstBuffer[i]) T(std::move(srcBuffer[i]));
    }
  );
}

template <typename BufferType>
KOKKOS_INLINE_FUNCTION
bool is_last_reference(const BufferType& buffer)
{
  return (buffer.use_count() == 1);
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
ViewVector<T, MemSpace, SizeT>::ViewVector(const std::string& name)
  : m_buffer(Kokkos::view_alloc(name, Kokkos::WithoutInitializing), 0),
    m_size(name + " size")
{
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
ViewVector<T, MemSpace, SizeT>::ViewVector(const std::string& name, SizeT size)
  : m_buffer(Kokkos::view_alloc(name, Kokkos::WithoutInitializing), size),
    m_size(name + " size")
{
  construct_entries(m_buffer, 0, size);
  m_size() = size;
  Kokkos::fence();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION ViewVector<T, MemSpace, SizeT>::~ViewVector()
{
  KOKKOS_IF_ON_HOST(
    bool viewIsTracked = const_cast<Kokkos::Impl::SharedAllocationTracker&>(m_buffer.impl_track()).tracking_enabled();
    if (viewIsTracked && is_last_reference(m_buffer)) {
       destroy_entries(m_buffer, 0, size());
       Kokkos::fence();
    }
  )
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
SizeT ViewVector<T, MemSpace, SizeT>::size() const
{
  if (!m_size.is_allocated()) {
    return 0;
  }
  return m_size();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
SizeT ViewVector<T, MemSpace, SizeT>::capacity() const
{
  return m_buffer.extent(0);
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
bool ViewVector<T, MemSpace, SizeT>::empty() const
{
  return size() == 0u;
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
T& ViewVector<T, MemSpace, SizeT>::operator[](SizeT idx) const
{
  return m_buffer[idx];
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
T* ViewVector<T, MemSpace, SizeT>::data() const
{
  return m_buffer.data();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
std::string ViewVector<T, MemSpace, SizeT>::name() const
{
  return m_buffer.label();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::reserve(SizeT newCapacity)
{
  if (newCapacity > capacity()) {
    change_capacity(newCapacity);
  }
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::resize(SizeT newSize)
{
  if (newSize > capacity()) {
    change_capacity(newSize);
  }

  if (newSize > size()) {
    construct_entries(m_buffer, size(), newSize);
  }
  else if (newSize < size()) {
    destroy_entries(m_buffer, newSize, size());
  }
  Kokkos::fence();

  if (!m_size.is_allocated()) {
    m_size = size_view_type("size");
  }
  m_size() = newSize;
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::resize_scale(SizeT newSize)
{
  const SizeT oldCapacity = capacity();
  if (newSize > oldCapacity) {
    SizeT newCapacity = std::max<SizeT>(oldCapacity*growth_scale, 1u);
    while (newCapacity < newSize) {
      newCapacity *= growth_scale;
    }
    change_capacity(newCapacity);
  }

  if (newSize > size()) {
    construct_entries(m_buffer, size(), newSize);
  }
  else if (newSize < size()) {
    destroy_entries(m_buffer, newSize, size());
  }
  Kokkos::fence();

  if (!m_size.is_allocated()) {
    m_size = size_view_type("size");
  }
  m_size() = newSize;
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::clear()
{
  resize(0);
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::swap(ViewVector& other)
{
  std::swap(m_buffer, other.m_buffer);
  std::swap(m_size, other.m_size);
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::push_back(const T& value)
{
  static_assert(Kokkos::SpaceAccessibility<ngp::HostExecSpace, MemSpace>::accessible);
  KOKKOS_IF_ON_DEVICE((
    Kokkos::abort("Cannot call ViewVector::push_back() on device because it may need to allocate memory.");
  ))

  const SizeT oldCapacity = capacity();
  if (size() == oldCapacity) {
    SizeT newCapacity = std::max<SizeT>(oldCapacity*growth_scale, 1u);
    change_capacity(newCapacity);
  }

  new (&m_buffer[size()]) T(value);
  if (!m_size.is_allocated()) {
    m_size = size_view_type("size");
  }
  ++m_size();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::push_back(T&& value)
{
  static_assert(Kokkos::SpaceAccessibility<ngp::HostExecSpace, MemSpace>::accessible);
  KOKKOS_IF_ON_DEVICE((
    Kokkos::abort("Cannot call ViewVector::push_back() on device because it may need to allocate memory.");
  ))

  const SizeT oldCapacity = capacity();
  if (size() == oldCapacity) {
    SizeT newCapacity = std::max<SizeT>(oldCapacity*growth_scale, 1u);
    change_capacity(newCapacity);
  }

  new (&m_buffer[size()]) T(std::move(value));
  if (!m_size.is_allocated()) {
    m_size = size_view_type("size");
  }
  ++m_size();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
template <typename... Args>
void ViewVector<T, MemSpace, SizeT>::emplace_back(Args&&... args)
{
  static_assert(Kokkos::SpaceAccessibility<ngp::HostExecSpace, MemSpace>::accessible);
  KOKKOS_IF_ON_DEVICE((
    Kokkos::abort("Cannot call ViewVector::emplace_back() on device because it may need to allocate memory.");
  ))

  const SizeT oldCapacity = capacity();
  if (size() == oldCapacity) {
    SizeT newCapacity = std::max<SizeT>(oldCapacity*growth_scale, 1u);
    change_capacity(newCapacity);
  }

  new (&m_buffer[size()]) T(std::forward<Args>(args)...);
  if (!m_size.is_allocated()) {
    m_size = size_view_type("size");
  }
  ++m_size();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
auto& ViewVector<T, MemSpace, SizeT>::get_view() const
{
  return m_buffer;
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
KOKKOS_INLINE_FUNCTION
bool ViewVector<T, MemSpace, SizeT>::is_allocated() const
{
  return m_buffer.is_allocated() && m_size.is_allocated();
}

//------------------------------------------------------------------------------
template <typename T, typename MemSpace, typename SizeT>
void ViewVector<T, MemSpace, SizeT>::change_capacity(SizeT newCapacity)
{
  if (newCapacity == capacity()) return;
  STK_ThrowRequireMsg(m_buffer.use_count() <= 1, "When resizing a ViewVector with multiple copies, the capcity must be large enough to include the new size.  Make sure to preallocate a large enough buffer using `.reserve()` function before making any copies.");

  buffer_type newBuffer(Kokkos::view_alloc(m_buffer.label(), Kokkos::WithoutInitializing), newCapacity);

  const SizeT numToMove = std::min(size(), newCapacity);
  move_entries(m_buffer, newBuffer, 0, numToMove);
  destroy_entries(m_buffer, 0, size());
  Kokkos::fence();

  if (!m_size.is_allocated()) {
    m_size = size_view_type("size");
  }
  m_size() = numToMove;
  m_buffer = newBuffer;
}

template <typename T, typename DestMemSpace, typename SrcMemSpace, typename DestSizeT, typename SrcSizeT>
void deep_copy(ViewVector<T, DestMemSpace, DestSizeT>& destView, ViewVector<T, SrcMemSpace, SrcSizeT>& srcView)
{
  destView.resize(srcView.size());
  auto destViewSize = typename ViewVector<T, DestMemSpace, DestSizeT>::unmanaged_buffer_type(destView.data(), destView.size());
  auto srcViewSize = typename ViewVector<T, SrcMemSpace, SrcSizeT>::unmanaged_buffer_type(srcView.data(), srcView.size());
  Kokkos::deep_copy(destViewSize, srcViewSize);
}

//------------------------------------------------------------------------------
}  // namespace stk::mesh::impl

#endif // STK_MESH_VIEWVECTOR_HPP
