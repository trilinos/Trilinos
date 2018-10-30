// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STKTOPOLOGY_DETAIL_FILL_CONTAINER_HPP
#define STKTOPOLOGY_DETAIL_FILL_CONTAINER_HPP

#include <vector>

namespace stk { namespace topology_detail {

template <typename OrdinalOutputIterator>
struct fill_ordinal_container {

  template <typename Ordinal>
  void operator()(Ordinal i)
  { *m_itr = i; ++m_itr; }

  fill_ordinal_container( OrdinalOutputIterator itr)
    : m_itr(itr)
  {}

  OrdinalOutputIterator m_itr;
};

template <typename T, typename A>
struct fill_ordinal_container< std::vector<T,A> >
{
  template <typename Ordinal>
  void operator()(Ordinal i)
  { *m_itr = i; ++m_itr; }

  fill_ordinal_container( std::vector<T,A> & vec)
    : m_itr(vec.begin())
  {}

  typename std::vector<T,A>::iterator m_itr;

};

template <typename NodeArray, typename NodeOutputIterator>
struct fill_node_container {

  template <typename Ordinal>
  void operator()(Ordinal i)
  { *m_itr = m_nodes[i]; ++m_itr; }

  fill_node_container( const NodeArray & nodes, NodeOutputIterator itr)
    : m_nodes(nodes)
    , m_itr(itr)
  {}

  const NodeArray    & m_nodes;
  NodeOutputIterator   m_itr;
};

template <typename NodeArray, typename T, typename A>
struct fill_node_container<NodeArray, std::vector<T,A> > {

  template <typename Ordinal>
  void operator()(Ordinal i)
  { *m_itr = m_nodes[i]; ++m_itr; }

  STK_FUNCTION
  fill_node_container( const NodeArray & nodes, std::vector<T,A> & vec)
    : m_nodes(nodes)
    , m_itr(vec.begin())
  {}

  const NodeArray    & m_nodes;
  typename std::vector<T,A>::iterator   m_itr;
};

}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_FILL_UNSIGNED_CONTAINER_HPP

