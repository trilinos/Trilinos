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

#ifndef STK_UTIL_STK_UTIL_UTIL_ALLOCATOR_MEMORY_USAGE_HPP
#define STK_UTIL_STK_UTIL_UTIL_ALLOCATOR_MEMORY_USAGE_HPP


namespace stk {

// Based off of code example in "The C++ Standard Library - A Tutorial and Reference"

template <typename Tag = void>
struct allocator_memory_usage
{
  typedef allocator_memory_usage<void> default_usage;

  static void allocate(size_t num_bytes);
  static void deallocate(size_t num_bytes);

  static size_t peak_memory() { return m_peak_memory; }
  static size_t current_memory() { return m_current_memory; }
  static size_t num_allocations() { return m_num_allocations; }
  static size_t num_deallocations() { return m_num_deallocations; }

  static void reset()
  {
    m_peak_memory = 0u;
    m_current_memory = 0u;
    m_num_allocations = 0u;
    m_num_allocations = 0u;
  }

private:
  static size_t m_peak_memory;
  static size_t m_current_memory;
  static size_t m_num_allocations;
  static size_t m_num_deallocations;
};

// initialize static members
template <typename Tag>
size_t allocator_memory_usage<Tag>::m_peak_memory = 0;

template <typename Tag>
size_t allocator_memory_usage<Tag>::m_current_memory = 0;

template <typename Tag>
size_t allocator_memory_usage<Tag>::m_num_allocations = 0;

template <typename Tag>
size_t allocator_memory_usage<Tag>::m_num_deallocations = 0;


template<>
inline void allocator_memory_usage<void>::allocate(size_t num_bytes)
{
  ++m_num_allocations;
  m_current_memory += num_bytes;
  m_peak_memory =  (m_current_memory > m_peak_memory) ? m_current_memory : m_peak_memory;
}

template<>
inline void allocator_memory_usage<void>::deallocate(size_t num_bytes)
{
  m_current_memory -= num_bytes;
  ++m_num_deallocations;
}

template<typename Tag>
inline void allocator_memory_usage<Tag>::allocate(size_t num_bytes)
{
  allocator_memory_usage<void>::allocate(num_bytes);
  ++m_num_allocations;
  m_current_memory += num_bytes;
  m_peak_memory =  (m_current_memory > m_peak_memory) ? m_current_memory : m_peak_memory;
}

template<typename Tag>
inline void allocator_memory_usage<Tag>::deallocate(size_t num_bytes)
{
  allocator_memory_usage<void>::deallocate(num_bytes);
  m_current_memory -= num_bytes;
  ++m_num_deallocations;
}

} // namespace stk



#endif //STK_UTIL_STK_UTIL_UTIL_ALLOCATOR_MEMORY_USAGE_HPP
