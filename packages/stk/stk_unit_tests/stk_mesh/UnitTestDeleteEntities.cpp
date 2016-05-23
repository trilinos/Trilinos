#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include "../../stk_mesh/stk_mesh/base/FEMHelpers.hpp"
#include "../../stk_mesh/stk_mesh/base/GetEntities.hpp"

namespace
{

void expect_num_elements_and_nodes(stk::mesh::BulkData &bulk, size_t goldNumElems, size_t goldNumNodes)
{
    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(bulk, entityCounts);
    EXPECT_EQ(goldNumElems, entityCounts[stk::topology::ELEM_RANK]);
    EXPECT_EQ(goldNumNodes, entityCounts[stk::topology::NODE_RANK]);
}

void expect_num_elements_and_faces_and_nodes(stk::mesh::BulkData &bulk, size_t goldNumElems, size_t goldNumFaces, size_t goldNumNodes)
{
    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(bulk, entityCounts);
    EXPECT_EQ(goldNumElems, entityCounts[stk::topology::ELEM_RANK]);
    EXPECT_EQ(goldNumFaces, entityCounts[stk::topology::FACE_RANK]);
    EXPECT_EQ(goldNumNodes, entityCounts[stk::topology::NODE_RANK]);
}

class HexShellHexMesh : public stk::unit_test_util::MeshFixture
{
protected:
    HexShellHexMesh()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        std::string meshDesc =
            "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
             1,2,HEX_8,5,6,7,8,9,10,11,12\n\
             0,3,SHELL_QUAD_4,5,6,7,8";
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }
};

TEST_F(HexShellHexMesh, DeleteShell_OnlyHexesRemain)
{
    expect_num_elements_and_nodes(get_bulk(), 3u, 12u);
    get_bulk().destroy_elements_of_topology(stk::topology::SHELL_QUAD_4);
    expect_num_elements_and_nodes(get_bulk(), 2u, 12u);
}

class HexHexShellMesh : public stk::unit_test_util::MeshFixture
{
protected:
    HexHexShellMesh()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        std::string meshDesc =
            "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
             1,2,HEX_8,5,6,7,8,9,10,11,12\n\
             0,3,SHELL_QUAD_4,9,10,11,12";
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }
};

TEST_F(HexHexShellMesh, DeleteAllHexes_OnlyShellRemains)
{
    expect_num_elements_and_nodes(get_bulk(), 3u, 12u);
    get_bulk().destroy_elements_of_topology(stk::topology::HEX_8);
    expect_num_elements_and_nodes(get_bulk(), 1u, 4u);

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(get_bulk(), stk::topology::NODE_RANK, nodes);
    for(stk::mesh::Entity node : nodes)
        EXPECT_TRUE(!get_bulk().bucket(node).shared());
}

class HexWedgeHexMesh : public stk::unit_test_util::MeshFixture
{
protected:
    HexWedgeHexMesh()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        std::string meshDesc =
            "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
             1,2,WEDGE_6,5,9,8,6,10,7\n\
             1,3,HEX_8,11,12,13,14,5,9,10,6";
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }
};

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteAllHexes_OnlyWedgeRemains)
{
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
    get_bulk().destroy_elements_of_topology(stk::topology::HEX_8);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 1u, 5u, 6u);
}

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteWedge_TwoHexesRemain)
{
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
    get_bulk().destroy_elements_of_topology(stk::topology::WEDGE_6);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 2u, 12u, 14u);
}

void add_hexes_back(stk::mesh::BulkData &bulk)
{
    bulk.modification_begin();
    stk::mesh::PartVector topologyParts = {&bulk.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8)};
    if(bulk.parallel_rank() == 0 || bulk.parallel_size() == 1)
        stk::mesh::declare_element(bulk, topologyParts, 1, {1,2,3,4,5,6,7,8});
    if(bulk.parallel_rank() == 1 || bulk.parallel_size() == 1)
        stk::mesh::declare_element(bulk, topologyParts, 3, {11,12,13,14,5,9,10,6});
    if(bulk.parallel_size() > 1 && bulk.parallel_rank() < 2)
    {
        int otherProcRank = bulk.parallel_rank() == 0 ? 1 : 0;
        bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 5), otherProcRank);
        bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 6), otherProcRank);
        bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 7), otherProcRank);
        bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 8), otherProcRank);
    }
    bulk.modification_end();
}

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteAllHexesThenCreateFaces_OnlyWedgeRemains)
{
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);

    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);

    get_bulk().destroy_elements_of_topology(stk::topology::HEX_8);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 1u, 5u, 6u);

    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 1u, 5u, 6u);

    add_hexes_back(get_bulk());
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 5u, 14u);

//    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
//    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
}

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteWedgeThenCreateFaces_TwoHexesRemain)
{
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
    get_bulk().destroy_elements_of_topology(stk::topology::WEDGE_6);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 2u, 12u, 14u);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
    expect_num_elements_and_faces_and_nodes(get_bulk(), 2u, 12u, 14u);
}

class SingleHexMesh : public stk::unit_test_util::MeshFixture
{
protected:
    void create_hex_with_one_face_on_proc_zero()
    {
        get_bulk().modification_begin();
        stk::mesh::PartVector topologyParts = {&get_meta().get_topology_root_part(stk::topology::HEX_8)};
        if(get_bulk().parallel_rank() == 0)
        {
            stk::mesh::Entity elem = stk::mesh::declare_element(get_bulk(), topologyParts, 1, {1,2,3,4,5,6,7,8});
            stk::mesh::declare_element_side(get_bulk(), elem, 5, {});
        }
        get_bulk().modification_end();
    }
    void create_adjacent_hex_on_proc_one()
    {
        stk::mesh::PartVector topologyParts = {&get_meta().get_topology_root_part(stk::topology::HEX_8)};
        get_bulk().modification_begin();
        if(get_bulk().parallel_rank() == 1)
            stk::mesh::declare_element(get_bulk(), topologyParts, 2, {5,6,7,8,9,10,11,12});

        int otherProcRank = get_bulk().parallel_rank() == 0 ? 1 : 0;
        get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, 5), otherProcRank);
        get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, 6), otherProcRank);
        get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, 7), otherProcRank);
        get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, 8), otherProcRank);
        get_bulk().modification_end();
    }
    void expect_face_connected_to_element_with_id(stk::mesh::EntityId id)
    {
        if(get_bulk().parallel_rank() == 1)
        {
            stk::mesh::EntityVector faces;
            stk::mesh::get_entities(get_bulk(), stk::topology::FACE_RANK, faces);
            ASSERT_EQ(1u, faces.size());
            unsigned numElems = get_bulk().num_elements(faces[0]);
            ASSERT_EQ(1u, numElems);
            const stk::mesh::Entity * elems = get_bulk().begin_elements(faces[0]);
            EXPECT_EQ(id, get_bulk().identifier(elems[0]));
        }
    }
};
// pre-existing face is not being attached to newly created element
TEST_F(SingleHexMesh, DISABLED_CreateFacesThenCreateAnotherElement_ConnectivityIsWrong)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        create_hex_with_one_face_on_proc_zero();
        create_adjacent_hex_on_proc_one();
        expect_face_connected_to_element_with_id(get_bulk().parallel_rank() + 1);
    }
}

}
