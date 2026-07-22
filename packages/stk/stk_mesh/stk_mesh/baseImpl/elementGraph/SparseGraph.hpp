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

#ifndef STK_ELEM_SPARSE_GRAPH_HPP
#define STK_ELEM_SPARSE_GRAPH_HPP

#include <vector>
#include <map>
#include "ElemElemGraphImpl.hpp"

namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }
namespace stk { namespace mesh { namespace impl { class IdMapper; } } }

namespace stk
{
namespace mesh
{
namespace impl
{

class SparseGraph
{
public:
    typedef std::map<LocalId, std::vector<GraphEdge>> ElemToCoincidentEdgesMap;
    typedef ElemToCoincidentEdgesMap::value_type value_type;

    size_t get_num_elements_in_graph() const { return m_graphEdges.size(); }
    const std::vector<GraphEdge> & get_edges_for_element(impl::LocalId elem) const;
    void reserve_edges(impl::LocalId /*localElemId*/, size_t /*numEdges*/) {}
    void add_edge(const GraphEdge &graphEdge);
    void delete_edge(const GraphEdge &graphEdge);
    void delete_all_edges(impl::LocalId elem);
    void clear()
    {
        m_graphEdges.clear();
    }

    ElemToCoincidentEdgesMap::const_iterator begin() const { return m_graphEdges.begin(); }
    ElemToCoincidentEdgesMap::const_iterator end() const { return m_graphEdges.end(); }

private:
    ElemToCoincidentEdgesMap m_graphEdges;
    static const std::vector<GraphEdge> emptyVector;

    void remove_coincident_edge(const GraphEdge& graphEdge, std::vector<GraphEdge>& graphEdgesForElem, size_t i);
};


}}} // end namespaces stk mesh

#endif
