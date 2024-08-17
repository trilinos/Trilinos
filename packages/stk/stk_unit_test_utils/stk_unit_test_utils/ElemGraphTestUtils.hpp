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

#ifndef Stk_UnitTest_ElemDeath_Utils
#define Stk_UnitTest_ElemDeath_Utils

#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

namespace ElemGraphTestUtils
{

class ElementDeathBulkDataTester : public stk::mesh::BulkData
{
public:

    ElementDeathBulkDataTester(stk::mesh::MetaData &mesh_meta_data,
                               MPI_Comm comm,
                               enum stk::mesh::BulkData::AutomaticAuraOption auraOption) :
      stk::mesh::BulkData(std::shared_ptr<stk::mesh::MetaData>(&mesh_meta_data, [](auto pointerWeWontDelete){}), comm, auraOption)
    {
    }

    void my_de_induce_unranked_part_from_nodes(const stk::mesh::EntityVector & deactivatedElements,
                                               stk::mesh::Part & activePart)
    {
        this->de_induce_parts_from_nodes(deactivatedElements, activePart);
    }
    void my_remove_boundary_faces_from_part(stk::mesh::ElemElemGraph &graph,
                                            const stk::mesh::EntityVector & deactivatedElements,
                                            stk::mesh::Part & activePart)
    {
        this->remove_boundary_faces_from_part(graph, deactivatedElements, activePart);
    }
};

inline void deactivate_elements(const stk::mesh::EntityVector &deactivated_elems, stk::mesh::BulkData &bulkData, stk::mesh::Part& active)
{
    bulkData.modification_begin();

    for(size_t i = 0; i < deactivated_elems.size(); ++i)
    {
        bulkData.change_entity_parts(deactivated_elems[i], stk::mesh::PartVector(), stk::mesh::PartVector(1, &active));
    }

    bulkData.modification_end();
}

inline int get_side_between_elements(const stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& graph, stk::mesh::Entity elem1, stk::mesh::EntityId elem2Id)
{
    size_t numConnected = graph.get_num_connected_elems(elem1);
    for(size_t i=0; i<numConnected; ++i)
    {
        stk::mesh::EntityId id = 0;
        int side = -1;
        if (graph.is_connected_elem_locally_owned(elem1, i))
        {
            stk::mesh::impl::ElementViaSidePair elemViaSidePair = graph.get_connected_element_and_via_side(elem1, i);
            id = bulkData.identifier(elemViaSidePair.element);
            side = elemViaSidePair.side;
        }
        else
        {
            stk::mesh::impl::IdViaSidePair idViaSidePair = graph.get_connected_remote_id_and_via_side(elem1, i);
            id = idViaSidePair.id;
            side = idViaSidePair.side;
        }
        if (id == elem2Id) {
            return side;
        }
    }
    return -1;
}

inline stk::mesh::Entity get_face_between_element_ids(stk::mesh::ElemElemGraph& graph, stk::mesh::BulkData& bulkData, stk::mesh::EntityId elem1Id, stk::mesh::EntityId elem2Id)
{
    stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, elem1Id);
    stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem2Id);

    bool isElem1LocallyOwnedAndValid = bulkData.is_valid(elem1) && bulkData.bucket(elem1).owned();
    bool isElem2LocallyOwnedAndValid = bulkData.is_valid(elem2) && bulkData.bucket(elem2).owned();

    stk::mesh::Entity face_between_elem1_and_elem2;

    if(isElem1LocallyOwnedAndValid && isElem2LocallyOwnedAndValid)
    {
        int side = get_side_between_elements(bulkData, graph, elem1, elem2Id);
        EXPECT_TRUE(side != -1);
        face_between_elem1_and_elem2 = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, elem1, side);
    }
    else if(isElem1LocallyOwnedAndValid)
    {
        int side = get_side_between_elements(bulkData, graph, elem1, elem2Id);
        EXPECT_TRUE(side != -1);
        face_between_elem1_and_elem2 = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, elem1, side);
    }
    else if(isElem2LocallyOwnedAndValid)
    {
        int side = get_side_between_elements(bulkData, graph, elem2, elem1Id);
        EXPECT_TRUE(side != -1);
        face_between_elem1_and_elem2 = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, elem2, side);
    }
    return face_between_elem1_and_elem2;
}

inline void skin_boundary(stk::mesh::BulkData& bulkData, stk::mesh::Part &partToSkin, const stk::mesh::PartVector& putSkinInTheseParts)
{
    stk::mesh::Selector sel = partToSkin;
    stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, putSkinInTheseParts);
}

inline void skin_part(stk::mesh::BulkData& bulkData, const stk::mesh::Part &partToSkin, const stk::mesh::PartVector& putSkinInTheseParts)
{
    stk::mesh::Selector sel = partToSkin;
    stk::mesh::Selector air = !partToSkin;
    stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, putSkinInTheseParts, air);
}

inline void test_num_faces_on_this_element(const stk::mesh::BulkData& bulkData, stk::mesh::EntityId id, size_t gold_num_faces_this_elem)
{
    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, id);
    if(bulkData.is_valid(element))
    {
        unsigned num_faces_this_elem = bulkData.num_faces(element);
        EXPECT_EQ(gold_num_faces_this_elem, num_faces_this_elem);
    }
}

inline void test_num_faces_per_element(const stk::mesh::BulkData& bulkData, const std::vector<size_t>& gold_num_faces_per_elem)
{
    for(size_t i=0;i<gold_num_faces_per_elem.size();++i)
    {
        stk::mesh::EntityId element_id = i+1;
        test_num_faces_on_this_element(bulkData, element_id, gold_num_faces_per_elem[i]);
    }
}

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
ElementDeathBulkDataTester : public ElemGraphTestUtils::ElementDeathBulkDataTester {};

} // namespace simple_fields

} // end namespace


#endif
