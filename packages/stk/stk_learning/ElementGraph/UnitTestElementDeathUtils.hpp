#ifndef Stk_UnitTest_ElemDeath_Utils
#define Stk_UnitTest_ElemDeath_Utils

#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

namespace ElementDeathUtils
{

class ElementDeathBulkDataTester : public stk::mesh::BulkData
{
public:

    ElementDeathBulkDataTester(stk::mesh::MetaData &mesh_meta_data,
                               MPI_Comm comm,
                               enum stk::mesh::BulkData::AutomaticAuraOption auraOption) :
            stk::mesh::BulkData(mesh_meta_data, comm, auraOption)
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
        face_between_elem1_and_elem2 = stk::mesh::impl::get_side_for_element(bulkData, elem1, side);
    }
    else if(isElem1LocallyOwnedAndValid)
    {
        int side = get_side_between_elements(bulkData, graph, elem1, elem2Id);
        EXPECT_TRUE(side != -1);
        face_between_elem1_and_elem2 = stk::mesh::impl::get_side_for_element(bulkData, elem1, side);
    }
    else if(isElem2LocallyOwnedAndValid)
    {
        int side = get_side_between_elements(bulkData, graph, elem2, elem1Id);
        EXPECT_TRUE(side != -1);
        face_between_elem1_and_elem2 = stk::mesh::impl::get_side_for_element(bulkData, elem2, side);
    }
    return face_between_elem1_and_elem2;
}

inline void skin_boundary(stk::mesh::BulkData& bulkData, stk::mesh::Part &active, const stk::mesh::PartVector& skin_parts)
{
    stk::mesh::Selector sel = active;
    stk::mesh::ElemElemGraph elem_elem_graph(bulkData, sel);
    elem_elem_graph.skin_mesh(skin_parts);
}

inline void skin_part(stk::mesh::BulkData& bulkData, const stk::mesh::Part &active, const stk::mesh::PartVector& skin_parts)
{
    stk::mesh::Selector sel = active;
    stk::mesh::Selector air = !active;

    stk::mesh::ElemElemGraph elem_elem_graph(bulkData, sel, &air);
    elem_elem_graph.skin_mesh(skin_parts);
}

} // end namespace


#endif
