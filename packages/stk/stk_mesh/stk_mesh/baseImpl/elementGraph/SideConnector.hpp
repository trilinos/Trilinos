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

#ifndef STK_SIDE_CONNECTOR_HPP
#define STK_SIDE_CONNECTOR_HPP

#include <vector>
#include <stk_mesh/base/Types.hpp>
#include "BulkDataIdMapper.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "GraphEdgeData.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace impl { class ElementLocalIdMapper; } } }

namespace stk
{
namespace mesh
{

class SideCreationElementChooser
{
public:
    SideCreationElementChooser(const stk::mesh::BulkData &b,
                               const stk::mesh::impl::ElementLocalIdMapper &im,
                               const stk::mesh::Graph &g,
                               const stk::mesh::impl::SparseGraph &cg)
    : bulk(b), localIdMapper(im), graph(g), coincidentGraph(cg) { }

    GraphEdge get_chosen_graph_edge(stk::mesh::Entity elem, int elemSide) const;
private:
    const stk::mesh::BulkData &bulk;
    const stk::mesh::impl::ElementLocalIdMapper &localIdMapper;
    const stk::mesh::Graph &graph;
    const stk::mesh::impl::SparseGraph &coincidentGraph;
};

class SideIdChooser
{
public:
    SideIdChooser(const stk::mesh::BulkData &b,
                  const stk::mesh::impl::ElementLocalIdMapper &im,
                  const stk::mesh::Graph &g,
                  const stk::mesh::impl::SparseGraph &cg)
    : elemChooser(b, im, g, cg), idMapper(b, im) { }
    stk::mesh::EntityId get_chosen_side_id(stk::mesh::Entity elem, int elemSide) const;
private:
    const stk::mesh::SideCreationElementChooser elemChooser;
    stk::mesh::impl::BulkDataIdMapper idMapper;
};

class SideNodeConnector
{
public:
    SideNodeConnector(stk::mesh::BulkData &b,
                      const stk::mesh::Graph &g,
                      const stk::mesh::impl::SparseGraph &cg,
                      const stk::mesh::ParallelInfoForGraphEdges &p,
                      const stk::mesh::impl::ElementLocalIdMapper & lm)
    : bulk(b), graph(g), coincidentGraph(cg), parallelGraph(p), localMapper(lm),
      m_scratchOrdinals1(), m_scratchOrdinals2(), m_scratchOrdinals3()
    { }
    void connect_side_to_nodes(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide);
private:
    void connect_side_to_elements_nodes(stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide);
    void connect_side_to_other_elements_nodes(const GraphEdge &edgeWithMinId, stk::mesh::Entity sideEntity, stk::mesh::Entity elemEntity, int elemSide);
    stk::mesh::EntityVector get_permuted_side_nodes(stk::mesh::Entity elemEntity, int elemSide, const stk::mesh::EntityVector &sideNodes, int permutation);
    void declare_relations_to_nodes(stk::mesh::Entity sideEntity, const stk::mesh::EntityVector &sideNodes);
private:
    stk::mesh::BulkData &bulk;
    const stk::mesh::Graph &graph;
    const stk::mesh::impl::SparseGraph &coincidentGraph;
    const stk::mesh::ParallelInfoForGraphEdges &parallelGraph;
    const stk::mesh::impl::ElementLocalIdMapper &localMapper;
    OrdinalVector m_scratchOrdinals1, m_scratchOrdinals2, m_scratchOrdinals3;
};

class SideConnector
{
public:
    SideConnector(stk::mesh::BulkData &b,
                  const stk::mesh::Graph &g,
                  const stk::mesh::impl::SparseGraph &cg,
                  const stk::mesh::impl::ElementLocalIdMapper & localMapper) :
            m_bulk_data(b),
            m_graph(g),
            m_coincidentGraph(cg),
            m_localMapper(localMapper),
            m_scratchOrdinals1(), m_scratchOrdinals2(), m_scratchOrdinals3()
    {
    }

    void connect_side_to_all_elements(stk::mesh::Entity sideEntity,
                                      stk::mesh::Entity elemEntity,
                                      int elemSide);
private:
    SideConnector();

    void connect_side_entity_to_other_element(stk::mesh::Entity sideEntity,
                                              const stk::mesh::GraphEdge &graphEdge);
    void connect_side_to_coincident_elements(stk::mesh::Entity sideEntity, stk::mesh::impl::LocalId elemLocalId, int elemSide);
    stk::mesh::Permutation get_permutation_for_side(stk::mesh::Entity sideEntity,
                                                    stk::mesh::Entity element,
                                                    int sideOrd);
    stk::mesh::EntityVector get_nodes_of_elem_side(stk::mesh::Entity element, int sideOrd);
    void connect_side_to_elem(stk::mesh::Entity sideEntity,
                              stk::mesh::Entity otherElement,
                              int sideOrd);
    void connect_side_to_adjacent_elements(stk::mesh::Entity sideEntity, stk::mesh::impl::LocalId elemLocalId, int elemSide);
    stk::mesh::EntityId get_min_id(const GraphEdge& graphEdge,
                                   int elemSide,
                                   const stk::mesh::impl::BulkDataIdMapper& idMapper,
                                   stk::mesh::EntityId minId);

    stk::mesh::BulkData &m_bulk_data;
    const stk::mesh::Graph &m_graph;
    const stk::mesh::impl::SparseGraph &m_coincidentGraph;
    const stk::mesh::impl::ElementLocalIdMapper & m_localMapper;
    OrdinalVector m_scratchOrdinals1, m_scratchOrdinals2, m_scratchOrdinals3;
};

}
}

#endif
