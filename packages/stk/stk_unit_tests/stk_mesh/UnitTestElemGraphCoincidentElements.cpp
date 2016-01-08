#include "gtest/gtest.h"
#include <mpi.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

void setup_node_sharing(stk::mesh::BulkData &mesh, const std::vector< std::vector<unsigned> > & shared_nodeIDs_and_procs )
{
    const unsigned p_rank = mesh.parallel_rank();

    for (size_t nodeIdx = 0, end = shared_nodeIDs_and_procs.size(); nodeIdx < end; ++nodeIdx) {
        if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
            stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
            int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
            mesh.add_node_sharing(node, sharingProc);
        }
    }
}

class HexShellShell : public stk::unit_test_util::MeshFixture
{
protected:
    HexShellShell()
    {
        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    }

    void setup_hex_shell_shell_on_procs(std::vector<int> owningProcs)
    {
        stk::mesh::Part* hexPart = &get_meta().declare_part_with_topology("hex_part", stk::topology::HEX_8);
        stk::mesh::Part* shellPart = &get_meta().declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
        stk::mesh::PartVector parts = {hexPart, shellPart, shellPart};
        declare_elements_on_procs_and_setup_node_sharing(owningProcs, parts);
    }

private:
    void declare_elements_on_procs_and_setup_node_sharing(const std::vector<int>& owningProcs, const stk::mesh::PartVector& parts)
    {
        get_bulk().modification_begin();
        declare_elements_on_procs(owningProcs, parts);
        setup_node_sharing(get_bulk(), shared_nodeIDs_and_procs);
        get_bulk().modification_end();
    }

    void declare_elements_on_procs(const std::vector<int>& owningProcs, const stk::mesh::PartVector& parts)
    {
        for(size_t i = 0; i < nodeIDs.size(); ++i)
            if(owningProcs[i] == stk::parallel_machine_rank(get_comm()))
                stk::mesh::declare_element(get_bulk(), *parts[i], elemIDs[i], nodeIDs[i]);
    }

    std::vector<stk::mesh::EntityIdVector> nodeIDs { {1, 2, 3, 4, 5, 6, 7, 8}, {5, 6, 7, 8}, {5, 6, 7, 8}};
    stk::mesh::EntityIdVector elemIDs = {1, 2, 3};
    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector<std::vector<unsigned> > shared_nodeIDs_and_procs = {
             {0, 5, 1}, {0, 6, 1}, {0, 7, 1}, {0, 8, 1},  // proc 0
             {1, 5, 0}, {1, 6, 0}, {1, 7, 0}, {1, 8, 0}}; // proc 1
};


TEST_F(HexShellShell, Hex0Shell1Shell1Parallel)
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.1|
    //       |   |          |3.1|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    if(stk::parallel_machine_size(get_comm()) == 2u)
    {
        setup_hex_shell_shell_on_procs({0, 1, 1});

        stk::mesh::ElemElemGraph elemElemGraph(get_bulk(), get_meta().universal_part());

        if(stk::parallel_machine_rank(get_comm()) == 0)
        {
            const stk::mesh::Entity hex1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
            ASSERT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
            EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
            EXPECT_EQ(5, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
            EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
            EXPECT_EQ(5, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);
            EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
            EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));
            EXPECT_EQ(2u, elemElemGraph.num_edges());
            EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
        }
        else
        {
            const stk::mesh::Entity shell2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
            ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
            EXPECT_EQ(1, elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).side);
            EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).id);
            EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

            const stk::mesh::Entity shell3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
            ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell3));
            EXPECT_EQ(1, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
            EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
            EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        }
    }
}

TEST_F(HexShellShell, Hex0Shell0Shell1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |3.1|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    if(stk::parallel_machine_size(get_comm()) == 2u)
    {
        setup_hex_shell_shell_on_procs({0, 0, 1});

        stk::mesh::ElemElemGraph elemElemGraph(get_bulk(), get_meta().universal_part());

        if(stk::parallel_machine_rank(get_comm()) == 0)
        {
            const stk::mesh::Entity hex1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
            const stk::mesh::Entity shell2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
            ASSERT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
            ASSERT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
            EXPECT_EQ(5, elemElemGraph.get_connected_element_and_via_side(hex1, 0).side);
            EXPECT_EQ(shell2, elemElemGraph.get_connected_element_and_via_side(hex1, 0).element);
            ASSERT_TRUE(!elemElemGraph.is_connected_elem_locally_owned(hex1, 1));
            EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
            EXPECT_EQ(5, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);

            ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
            EXPECT_EQ(1, elemElemGraph.get_connected_element_and_via_side(shell2, 0).side);
            EXPECT_EQ(hex1, elemElemGraph.get_connected_element_and_via_side(shell2, 0).element);
            EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
            EXPECT_EQ(3u, elemElemGraph.num_edges());
            EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
        }
        else
        {
            const stk::mesh::Entity shell3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
            ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell3));
            EXPECT_EQ(1, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
            EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
            EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
            EXPECT_EQ(1u, elemElemGraph.num_edges());
            EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
        }
    }
}

TEST_F(HexShellShell, Skin)
{
    if(stk::parallel_machine_size(get_comm()) == 2u)
    {
        setup_hex_shell_shell_on_procs({0, 1, 0});

        stk::mesh::ElemElemGraph elemElemGraph(get_bulk(), get_meta().universal_part());
        elemElemGraph.skin_mesh({});

        stk::mesh::Selector ownedOrShared = get_meta().locally_owned_part() | get_meta().globally_shared_part();

        if(stk::parallel_machine_rank(get_comm()) == 0)
            EXPECT_EQ(6u, stk::mesh::count_selected_entities(ownedOrShared, get_bulk().buckets(stk::topology::FACE_RANK)));
        else
            EXPECT_EQ(1u, stk::mesh::count_selected_entities(ownedOrShared, get_bulk().buckets(stk::topology::FACE_RANK)));
    }
}

void expect_correct_connected_element_via_side(stk::mesh::ElemElemGraph& elemElemGraph, stk::mesh::Entity elem, int k, stk::mesh::Entity otherElem, int viaSide)
{
    stk::mesh::impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(elem, k);
    EXPECT_EQ(viaSide,      elem_via_side.side);
    EXPECT_EQ(otherElem, elem_via_side.element);
}


TEST( ElementGraph, HexAddShellAddShellSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |3.0|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2, 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    mesh.modification_begin();
    stk::mesh::EntityVector added_shells;
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
    }
    mesh.modification_end();

    elemElemGraph.add_elements(added_shells);

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell2, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell3, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
    expect_correct_connected_element_via_side(elemElemGraph, shell2, 0, hex1, 1);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexAddShellAddShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Added two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    mesh.modification_begin();
    stk::mesh::EntityVector added_shells;
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
    }
    mesh.modification_end();

    elemElemGraph.add_elements(added_shells);

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell3, 4);
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 4);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 1, hex2, 0);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 1);
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 1, hex2, 0);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellShellSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |3.0|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2, 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell2, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell3, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
    expect_correct_connected_element_via_side(elemElemGraph, shell2, 0, hex1, 1);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

    EXPECT_EQ(4u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell3, 4);
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 4);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 1, hex2, 0);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 1);
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 1, hex2, 0);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

    EXPECT_EQ(8u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellReversedShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 8, 7, 6 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell3, 4);
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 4);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 1, hex2, 0);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 0);
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 1, hex2, 1);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

    EXPECT_EQ(8u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0Shell0Shell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
        expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
        EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
        EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

        // Connectivity for Shell Element 4
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
        expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 1);
        EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
        EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

        EXPECT_EQ(6u, elemElemGraph.num_edges());
        EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,      elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
        EXPECT_EQ(4,      elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
        EXPECT_EQ(3u,   elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
        EXPECT_EQ(4u,   elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

        EXPECT_EQ(2u, elemElemGraph.num_edges());
        EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
    }
}

TEST( ElementGraph, Hex0Shell0Shell1Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |4.1|          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
        EXPECT_EQ(5,      elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);
        EXPECT_EQ(4u,     elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
        EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
        EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 4
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
        expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex2, 0);
        EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
        EXPECT_EQ(1u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

        // Connectivity for Hex Element 2
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
        expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell4, 4);
        EXPECT_EQ(4,      elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
        EXPECT_EQ(3u,     elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }

    EXPECT_EQ(4u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0Shell0ReversedShell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells, opposite orientation

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 8, 7, 6 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
        expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
        EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
        EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

        // Connectivity for Shell Element 4
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
        expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 0);
        EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
        EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
        EXPECT_EQ(6u, elemElemGraph.num_edges());
        EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(3u,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
        EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
        EXPECT_EQ(4u,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
        EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
        EXPECT_EQ(2u, elemElemGraph.num_edges());
        EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
    }
}

TEST( ElementGraph, Hex1Shell0Shell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.1-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.1-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.1    |3.0|   2.1    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.1---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.1------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 1, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh, meta.universal_part());

    const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Shell Element 3
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
        EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
        EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
        EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

        // Connectivity for Shell Element 4
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell4));
        EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell4, 0).side);
        EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
        EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell4, 0).id);
        EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
        EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);
        EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
        EXPECT_EQ(4u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Hex Element 2
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
        EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
        EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
        EXPECT_EQ(4u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }

    EXPECT_EQ(4u, elemElemGraph.num_edges());
    EXPECT_EQ(4u, elemElemGraph.num_parallel_edges());
}

}
