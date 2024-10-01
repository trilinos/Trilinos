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

#ifndef STK_ELEM_GRAPH_ELEM_FILTER_HPP
#define STK_ELEM_GRAPH_ELEM_FILTER_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include "BulkDataIdMapper.hpp"
#include <vector>

namespace stk { namespace mesh {
class BulkData;

namespace impl {

inline bool is_solid_shell_connection(stk::topology t1, stk::topology t2)
{
    return (!t1.is_shell() && t2.is_shell());
}

inline bool is_shell_solid_connection(stk::topology t1, stk::topology t2)
{
    return (t1.is_shell() && !t2.is_shell());
}

class ElementFilter
{
public:
  ElementFilter(const stk::mesh::BulkData& bulk,
                const impl::ElementLocalIdMapper& localMapper,
                bool populateSideNodesForConnectedElems)
  : m_bulk(bulk), m_localMapper(localMapper), m_populateSideNodes(populateSideNodesForConnectedElems)
  {}

  template<typename SideData>
  void add_local_elements_to_connected_list(stk::topology elementTopology,
                                             const stk::mesh::EntityVector & elementsAttachedToSideNodes,
                                             const stk::mesh::EntityVector & sideNodes,
                                             std::vector<SideData> & connectedElements)
  {
      SideData connectedElemData;
      stk::mesh::EntityVector connectedElemSideNodes;
      for (const stk::mesh::Entity& otherElement : elementsAttachedToSideNodes)
      {
          impl::LocalId localId = m_localMapper.entity_to_local(otherElement);
          if (localId != impl::INVALID_LOCAL_ID)
          {
              stk::mesh::OrdinalAndPermutation connectedOrdAndPerm = stk::mesh::get_ordinal_and_permutation(m_bulk, otherElement, m_bulk.mesh_meta_data().side_rank(), sideNodes);
              if (INVALID_CONNECTIVITY_ORDINAL != connectedOrdAndPerm.first)
              {
                  const stk::mesh::Bucket & connectedBucket = m_bulk.bucket(otherElement);
      
                  STK_ThrowAssertMsg(connectedBucket.topology().side_topology(connectedOrdAndPerm.first).num_nodes() == sideNodes.size(),
                                "Error, number of nodes on sides of adjacent elements do not agree:  " <<
                                 sideNodes.size() << " != " << connectedBucket.topology().side_topology(connectedOrdAndPerm.first).num_nodes());
  
                  stk::topology connectedTopology = connectedBucket.topology();
                  stk::topology connectedSideTopology = connectedTopology.side_topology(connectedOrdAndPerm.first);
                  const bool isConnectingSolidToWrongSideOfShell = is_shell_solid_connection(elementTopology, connectedTopology)
                               && static_cast<unsigned>(connectedOrdAndPerm.second) < connectedSideTopology.num_positive_permutations();
                  if(!isConnectingSolidToWrongSideOfShell)
                  {
                      if(is_solid_shell_connection(elementTopology, connectedTopology))
                          connectedOrdAndPerm = flip_shell_to_get_opposing_normal(connectedOrdAndPerm, connectedSideTopology);
  
                      int sideIndex = connectedOrdAndPerm.first;
                      if (m_populateSideNodes) {
                        connectedElemSideNodes.resize(sideNodes.size());
                        const stk::mesh::Entity* connectedElemNodes = m_bulk.begin_nodes(otherElement);
                        connectedTopology.side_nodes(connectedElemNodes, sideIndex, connectedElemSideNodes.data());
                        connectedElemData.set_element_side_nodes(connectedElemSideNodes);
                      }
                      connectedElemData.set_element_local_id(localId);
                      connectedElemData.set_element_identifier(m_bulk.identifier(otherElement));
                      connectedElemData.set_element_topology(connectedTopology);
                      connectedElemData.set_element_side_index(sideIndex);
                      connectedElemData.set_permutation(connectedOrdAndPerm.second);
                      connectedElements.push_back(connectedElemData);
                  }
              }
          }
      }
  }

private:
  const stk::mesh::BulkData& m_bulk;
  const impl::ElementLocalIdMapper& m_localMapper;
  bool m_populateSideNodes;
};

template<typename SideData>
void get_elements_connected_via_sidenodes(const stk::mesh::BulkData& bulkData, stk::mesh::EntityId elementId, stk::topology elementTopology, const impl::ElementLocalIdMapper & localMapper,
                                                           const stk::mesh::EntityVector &sideNodesOfReceivedElement,
                                                           stk::mesh::EntityVector& scratchEntityVector,
                                                           std::vector<SideData>& connectedElementDataVector)
{
    impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, sideNodesOfReceivedElement.size(), sideNodesOfReceivedElement.data(), scratchEntityVector);
    stk::util::sort_and_unique(scratchEntityVector);
    bool populateSideNodesForConnectedElems = true;
    ElementFilter elemFilter(bulkData, localMapper, populateSideNodesForConnectedElems);
    elemFilter.add_local_elements_to_connected_list(elementTopology, scratchEntityVector, sideNodesOfReceivedElement, connectedElementDataVector);
}

template<typename SideData>
void get_elements_with_larger_ids_connected_via_sidenodes(const stk::mesh::BulkData& bulkData, stk::mesh::EntityId elementId,
                                                                           stk::topology elementTopology, const impl::ElementLocalIdMapper & localMapper,
                                                                           const stk::mesh::EntityVector &sideNodesOfReceivedElement,
                                                                           stk::mesh::EntityVector& localElementsConnectedToReceivedSideNodes,
                                                                           std::vector<SideData>& connectedElementDataVector)
{
    impl::find_entities_with_larger_ids_these_nodes_have_in_common_and_locally_owned(elementId, bulkData, stk::topology::ELEMENT_RANK, sideNodesOfReceivedElement.size(), sideNodesOfReceivedElement.data(), localElementsConnectedToReceivedSideNodes);
    bool populateSideNodesForConnectedElems = false;
    ElementFilter elemFilter(bulkData, localMapper, populateSideNodesForConnectedElems);
    elemFilter.add_local_elements_to_connected_list(elementTopology, localElementsConnectedToReceivedSideNodes, sideNodesOfReceivedElement, connectedElementDataVector);
}

} // end impl

}} // end stk mesh namespaces

#endif
