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

#ifndef elementgraphtester_hpp
#define elementgraphtester_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequire
#include <vector>                       // for allocator, vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>  // for ElemElemGraph
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>
namespace stk { namespace mesh { class Part; } }

struct GraphEdgeMock
{
    stk::mesh::EntityId element1;
    stk::mesh::EntityId element2;
    int sideOrdinalConnectingElement1ToElement2;
};
typedef std::vector<GraphEdgeMock> GraphEdges;


class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
    ElemElemGraphTester(stk::mesh::BulkData& bulkData) : ElemElemGraph(bulkData) {}

    virtual ~ElemElemGraphTester() {}

    void fill_graph() { ElemElemGraph::fill_graph(); }

    void write_graph(std::ostream& out) const { ElemElemGraph::write_graph(out); }

    stk::mesh::Graph & get_graph() { return m_graph; }
    stk::mesh::impl::ParallelGraphInfo & get_parallel_graph_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
    stk::mesh::ParallelInfoForGraphEdges& get_parallel_graph() { return m_parallelInfoForGraphEdges; }

    stk::mesh::BulkData & get_bulk_data() { return m_bulk_data; }
    size_t get_graph_size() { return m_graph.get_num_elements_in_graph(); }

    stk::mesh::EntityId get_entity_id(stk::mesh::impl::LocalId localId)
    {
        stk::mesh::EntityId id = -localId;
        if(localId >= 0)
            id = get_bulk_data().identifier(m_idMapper.local_to_entity(localId));
        return id;
    }

    int get_first_encountered_side_between_elems(stk::mesh::impl::LocalId elemId, stk::mesh::EntityId elem2Id)
    {
        for(size_t i=0; i<m_graph.get_num_edges_for_element(elemId); i++)
        {
            if(get_entity_id(m_graph.get_edge_for_element(elemId, i).elem2()) == elem2Id)
            {
                return m_graph.get_edge_for_element(elemId, i).side1();
            }
        }
        return -1;
    }

    int check_local_connectivity(stk::mesh::Entity elem1, stk::mesh::Entity elem2)
    {
        int side=-1;
        if (is_valid_graph_element(elem1) && is_valid_graph_element(elem2)) {
            side = get_first_encountered_side_between_elems(get_local_element_id(elem1), get_bulk_data().identifier(elem2));
        }
        return side;
    }

    int check_remote_connectivity(stk::mesh::Entity elem, stk::mesh::EntityId other_elem_id)
    {
        int side=-1;
        if (is_valid_graph_element(elem)) {
            side = get_first_encountered_side_between_elems(get_local_element_id(elem), other_elem_id);
        }
        return side;
    }
    int get_side_from_element1_to_element2(stk::mesh::EntityId elem1_id, stk::mesh::EntityId elem2_id)
    {
        int side = -1;
        stk::mesh::Entity elem1 = m_bulk_data.get_entity(stk::topology::ELEM_RANK, elem1_id);
        stk::mesh::Entity elem2 = m_bulk_data.get_entity(stk::topology::ELEM_RANK, elem2_id);
        bool isElem1Local = m_bulk_data.is_valid(elem1) && m_bulk_data.bucket(elem1).owned();
        bool isElem2Local = m_bulk_data.is_valid(elem2) && m_bulk_data.bucket(elem2).owned();

        ThrowRequire(isElem1Local);

        if(isElem2Local)
        {
            side = check_local_connectivity(elem1, elem2);
        }
        else
        {
            side = check_remote_connectivity(elem1, elem2_id);
        }

        return side;
    }
};

#endif
