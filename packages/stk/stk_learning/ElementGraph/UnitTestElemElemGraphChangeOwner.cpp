
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <ostream>                      // for basic_ostream::operator<<, etc
#include <stdexcept>                    // for logic_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/ElemElemGraph.hpp>  // for change_entity_owner, etc
#include <stk_mesh/base/ElemElemGraphImpl.hpp>  // for parallel_info
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include <stk_util/environment/WallTime.hpp>  // for wall_time
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <string>                       // for string
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "ElementGraphTester.hpp"       // for ElemElemGraphTester
#include "gtest/gtest-message.h"        // for Message
#include "mpi.h"                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace stk { namespace mesh { class Part; } }

namespace {

class ElemGraphChangeOwner : public stk::unit_test_util::MeshTestFixture
{
protected:
    typedef std::pair<stk::mesh::EntityId, int> EntityIdProc;
    typedef std::vector<EntityIdProc> EntityIdProcVector;

    void expect_initial_graph_correct()
    {
        if(get_bulk().parallel_rank() == 0)
            check_element2_connected_to_element1_locally_and_element3_remotely();
        else if(get_bulk().parallel_rank() == 1)
            check_element3_conencted_to_element4_locally_and_element2_remotely();
    }

    void check_element2_connected_to_element1_locally_and_element3_remotely()
    {
        stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
        ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem2));
        expect_connected_to_local_elem_id(elem2, 0, 1);
        expect_element2_connected_to_3_remotely_via_side_5();
    }

    void check_element3_conencted_to_element4_locally_and_element2_remotely()
    {
        stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
        ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem3));
        expect_connected_to_local_elem_id(elem3, 0, 4);
        expect_element3_connected_to_2_remotely_via_side_4();
    }

    void expect_element2_connected_to_3_remotely_via_side_5()
    {
        stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
        expect_connected_to_remote_elem_id(elem2, 1, 3);
        EXPECT_EQ(5, get_elem_graph().get_side_from_element1_to_remote_element2(elem2, 3));
        expect_parallel_info_from_elem2_to_3();
    }

    void expect_element3_connected_to_2_remotely_via_side_4()
    {
        stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
        expect_connected_to_remote_elem_id(elem3, 1, 2);
        EXPECT_EQ(4, get_elem_graph().get_side_from_element1_to_remote_element2(elem3, 2));
        expect_parallel_info_from_elem3_to_2();
    }

    void expect_parallel_info_from_elem2_to_3()
    {
        stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
        const stk::mesh::impl::parallel_info &parInfo = get_elem_graph().get_parallel_edge_info(elem2, 3);
        expect_otherProc_sideOrdinal_permutation_chosenId(parInfo, 1, 4, 4, 1);
        expect_parallel_info_part_memberships(parInfo);
    }

    void expect_parallel_info_from_elem3_to_2()
    {
        stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
        const stk::mesh::impl::parallel_info &parInfo = get_elem_graph().get_parallel_edge_info(elem3, 2);
        expect_otherProc_sideOrdinal_permutation_chosenId(parInfo, 0, 5, 4, 1);
        expect_parallel_info_part_memberships(parInfo);
    }

    void expect_connected_to_local_elem_id(stk::mesh::Entity elem,
                                           size_t connectedIndex,
                                           stk::mesh::EntityId connectedId)
    {
        stk::mesh::EntityId elemId = get_bulk().identifier(elem);
        ASSERT_TRUE(get_elem_graph().is_connected_elem_locally_owned(elem, connectedIndex))
                << "elem " << elemId << " expected local elem " << connectedId;
        EXPECT_EQ(connectedId, get_bulk().identifier(get_elem_graph().get_connected_element(elem, connectedIndex)))
                << "elem " << elemId;
    }

    void expect_connected_to_remote_elem_id(stk::mesh::Entity elem,
                                            size_t connectedIndex,
                                            stk::mesh::EntityId connectedId)
    {
        stk::mesh::EntityId elemId = get_bulk().identifier(elem);
        ASSERT_TRUE(!get_elem_graph().is_connected_elem_locally_owned(elem, connectedIndex))
                << "elem " << elemId << " expected remote elem " << connectedId;
        EXPECT_EQ(connectedId, get_elem_graph().get_entity_id_of_remote_element(elem, connectedIndex))
                << "elem " << elemId;
    }

    void expect_otherProc_sideOrdinal_permutation_chosenId(const stk::mesh::impl::parallel_info &parInfo,
                                                              int otherProc,
                                                              int sideOrdinal,
                                                              int perm,
                                                              stk::mesh::EntityId chosenId)
    {
        EXPECT_EQ(otherProc, parInfo.m_other_proc);
        EXPECT_EQ(sideOrdinal, parInfo.m_other_side_ord);
        EXPECT_EQ(perm, parInfo.m_permutation);
        EXPECT_EQ(chosenId, parInfo.m_chosen_side_id);
    }

    void expect_parallel_info_part_memberships(const stk::mesh::impl::parallel_info &parInfo)
    {
        EXPECT_TRUE(parInfo.m_in_body_to_be_skinned);
        EXPECT_FALSE(parInfo.m_is_air);
    }

    void move_elements(const EntityIdProcVector &elementIdProcsToMove)
    {
        stk::mesh::EntityProcVec elemProcPairsToMove;
        for(const EntityIdProc &entityIdProc : elementIdProcsToMove)
            append_element_if_owned(entityIdProc.first, entityIdProc.second, elemProcPairsToMove);
        change_entity_owner(get_bulk(), get_elem_graph(), elemProcPairsToMove);
    }

    void append_element_if_owned(stk::mesh::EntityId elementId, int destProc, stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elementId);
        if(is_owned_on_this_proc(element))
            elemProcPairsToMove.push_back(stk::mesh::EntityProc(element, destProc));
    }

    bool is_owned_on_this_proc(stk::mesh::Entity element)
    {
        return (get_bulk().is_valid(element) && get_bulk().bucket(element).owned());
    }

    void create_elem_graph()
    {
        elemElemGraph = new ElemElemGraphTester(get_bulk());
    }

    ElemElemGraphTester &get_elem_graph()
    {
        return *elemElemGraph;
    }

    ElemGraphChangeOwner() : elemElemGraph(nullptr)
    {
    }

    ~ElemGraphChangeOwner()
    {
        delete elemElemGraph;
    }

protected:
    ElemElemGraphTester *elemElemGraph;
};


class ElemGraphChangeOwnerMoveFrom1To0 : public ElemGraphChangeOwner
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh("generated:1x1x4", auraOption);
        test_graph_updated_with_elem_3_moving_to_proc_0();
    }

    void test_graph_updated_with_elem_3_moving_to_proc_0()
    {
        create_elem_graph();
        expect_initial_graph_correct();
        move_elements({EntityIdProc(3, 0)});
        expect_graph_updated_after_elem_3_moved_to_0();
    }

    void expect_graph_updated_after_elem_3_moved_to_0()
    {
        if(get_bulk().parallel_rank() == 0)
        {
            check_element2_connected_to_element1_and_element3_locally();
            check_element3_connected_to_element2_locally_and_element4_remotely();
        }
        else
        {
            check_element4_connected_to_element3_remotely();
        }
    }

    void check_element2_connected_to_element1_and_element3_locally()
    {
        stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
        ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem2));
        expect_connected_to_local_elem_id(elem2, 0, 1);
        expect_connected_to_local_elem_id(elem2, 1, 3);
    }

    void check_element4_connected_to_element3_remotely()
    {
        expect_element4_connected_to_3_remotely_via_side_3();
        expect_parallel_info_from_elem4_to_3();
    }

    void check_element3_connected_to_element2_locally_and_element4_remotely()
    {
        stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
        ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem3));
        expect_element3_connected_to_4_remotely_via_side_5();
        expect_connected_to_local_elem_id(elem3, 1, 2);
    }

    void expect_element3_connected_to_4_remotely_via_side_5()
    {
        stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
        expect_connected_to_remote_elem_id(elem3, 0, 4);
        EXPECT_EQ(5, get_elem_graph().get_side_from_element1_to_remote_element2(elem3, 4));
        expect_parallel_info_from_elem3_to_4();
    }

    void expect_parallel_info_from_elem3_to_4()
    {
        stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
        const stk::mesh::impl::parallel_info &parInfo = get_elem_graph().get_parallel_edge_info(elem3, 4);
        expect_otherProc_sideOrdinal_permutation_chosenId(parInfo, 1, 4, 4, 17);
        expect_parallel_info_part_memberships(parInfo);
    }

    void expect_element4_connected_to_3_remotely_via_side_3()
    {
        stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4);
        ASSERT_EQ(1u, get_elem_graph().get_num_connected_elems(elem4));
        expect_connected_to_remote_elem_id(elem4, 0, 3);
        EXPECT_EQ(4, get_elem_graph().get_side_from_element1_to_remote_element2(elem4, 3));
    }

    void expect_parallel_info_from_elem4_to_3()
    {
        stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4);
        const stk::mesh::impl::parallel_info &parInfo = get_elem_graph().get_parallel_edge_info(elem4, 3);
        expect_otherProc_sideOrdinal_permutation_chosenId(parInfo, 0, 5, 4, 17);
        expect_parallel_info_part_memberships(parInfo);
    }
};
TEST_F(ElemGraphChangeOwnerMoveFrom1To0, withAura)
{
    run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveFrom1To0, withoutAura)
{
    run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerMoveEverythingFromProc1 : public ElemGraphChangeOwner
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh("generated:1x1x4", auraOption);
        expect_graph_correct_after_moving_everything_to_proc0();
    }
    void expect_graph_correct_after_moving_everything_to_proc0()
    {
        create_elem_graph();
        expect_initial_graph_correct();
        move_elements({EntityIdProc(3, 0), EntityIdProc(4, 0)});

        ASSERT_TRUE(false) << "Need to add tests of graph data (expectations) if the code actually got this far.";
    }
};
TEST_F(ElemGraphChangeOwnerMoveEverythingFromProc1, DISABLED_withAura)
{
    run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveEverythingFromProc1, DISABLED_withoutAura)
{
    run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerLeapFrog : public ElemGraphChangeOwner
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh("generated:1x1x6", auraOption);
        expect_graph_correct_after_leaps();
    }
    void expect_graph_correct_after_leaps()
    {
        create_elem_graph();
        expect_initial_graph_correct();
        move_elements({EntityIdProc(2, 1), EntityIdProc(3, 2)});

        ASSERT_TRUE(false) << "Need to add tests of graph data (expectations) if the code actually got this far.";
    }
};
TEST_F(ElemGraphChangeOwnerLeapFrog, DISABLED_withAura)
{
    run_test_on_num_procs(3, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerLeapFrog, DISABLED_withoutAura)
{
    run_test_on_num_procs(3, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerMoveNeighborsToEnd : public ElemGraphChangeOwner
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh("generated:1x1x6", auraOption);
        expect_graph_correct_after_moving_neighbors_to_last_proc();
    }
    void expect_graph_correct_after_moving_neighbors_to_last_proc()
    {
        create_elem_graph();
        expect_initial_graph_correct();
        move_elements({EntityIdProc(2, 2), EntityIdProc(3, 2)});

        ASSERT_TRUE(false) << "Need to add tests of graph data (expectations) if the code actually got this far.";
    }
};
TEST_F(ElemGraphChangeOwnerMoveNeighborsToEnd, DISABLED_withAura)
{
    run_test_on_num_procs(3, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveNeighborsToEnd, DISABLED_withoutAura)
{
    run_test_on_num_procs(3, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerMoveNeighborsToDifferentProcs : public ElemGraphChangeOwner
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh_with_cyclic_decomp("1x2x2", auraOption);
        expect_mesh_created_correctly();
        expect_graph_correct_after_moving_neighbors_to_different_procs();
    }
    void expect_mesh_created_correctly()
    {
        expect_element_on_proc(1, 0);
        expect_element_on_proc(2, 1);
        expect_element_on_proc(3, 2);
        expect_element_on_proc(4, 3);
    }
    void expect_element_on_proc(stk::mesh::EntityId id, int proc)
    {
        if(get_bulk().parallel_rank() == proc)
            EXPECT_TRUE(is_owned_on_this_proc(get_bulk().get_entity(stk::topology::ELEMENT_RANK, id)));
    }
    void expect_graph_correct_after_moving_neighbors_to_different_procs()
    {
        create_elem_graph();
        check_initial_graph();
        move_elements({EntityIdProc(1, 1), EntityIdProc(3, 3)});

        ASSERT_TRUE(false) << "Need to add tests of graph data (expectations) if the code actually got this far.";
    }
    void check_initial_graph()
    {
        expect_element1_connected_correctly();
        expect_element2_connected_correctly();
    }

    void expect_element1_connected_correctly()
    {
        stk::mesh::Entity element1 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
        if(is_owned_on_this_proc(element1))
            expect_element1_connected_to_2_and_3(element1);
    }

    void expect_element1_connected_to_2_and_3(stk::mesh::Entity element1)
    {
        ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(element1));
        expect_connected_to_remote_elem_id_on_proc(element1, 0, 2, 1);
        expect_connected_to_remote_elem_id_on_proc(element1, 1, 3, 2);
    }

    void expect_element2_connected_correctly()
    {
        stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
        if(is_owned_on_this_proc(element2))
            expect_element2_connected_to_1_and_4(element2);
    }

    void expect_element2_connected_to_1_and_4(stk::mesh::Entity element2)
    {
        ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(element2));
        expect_connected_to_remote_elem_id_on_proc(element2, 0, 1, 0);
        expect_connected_to_remote_elem_id_on_proc(element2, 1, 4, 3);
    }

    void expect_connected_to_remote_elem_id_on_proc(stk::mesh::Entity elem,
                                                    size_t connectedIndex,
                                                    stk::mesh::EntityId connectedId,
                                                    int expectedProc)
    {
        expect_connected_to_remote_elem_id(elem, connectedIndex, connectedId);
        const stk::mesh::impl::parallel_info &parInfo = get_elem_graph().get_parallel_edge_info(elem, connectedId);
        EXPECT_EQ(expectedProc, parInfo.m_other_proc);
    }
};
TEST_F(ElemGraphChangeOwnerMoveNeighborsToDifferentProcs, DISABLED_withAura)
{
    run_test_on_num_procs(4, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveNeighborsToDifferentProcs, DISABLED_withoutAura)
{
    run_test_on_num_procs(4, stk::mesh::BulkData::NO_AUTO_AURA);
}



void change_entity_owner_hex_test_2_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);

    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(2, numLocallyOwnedElems);

        ElemElemGraphTester elem_graph(bulkData);

        // Create a vector of the elements to be moved
        std::vector <stk::mesh::Entity> elems_to_move;

        stk::mesh::EntityId elem_2_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

        if (proc == 0)
        {
            int side_from_elem2_to_elem3 = elem_graph.get_side_from_element1_to_remote_element2(elem_2, stk::mesh::EntityId(3));
            int side_from_elem2_to_elem1 = elem_graph.get_side_from_element1_to_locally_owned_element2(elem_2, bulkData.get_entity(stk::topology::ELEM_RANK,1));

            EXPECT_EQ(5, side_from_elem2_to_elem3);
            EXPECT_EQ(4, side_from_elem2_to_elem1);

            elems_to_move.push_back(elem_2);

            int other_proc = 1;
            for (unsigned i=0; i<elems_to_move.size(); i++)
            {
                EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
                EXPECT_EQ(0, bulkData.parallel_owner_rank(elems_to_move[i]));
                elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
            }
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move);

        elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

        if (proc == 0)
        {
            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            stk::mesh::impl::parallel_info &p_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));
            int other_side_ord = p_info.m_other_side_ord;
            EXPECT_EQ(4, other_side_ord);
            stk::mesh::EntityId chosen_face_id = 2;
            EXPECT_EQ(chosen_face_id, p_info.m_chosen_side_id);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(3)), std::logic_error);
        }

        if (proc == 1)
        {
            EXPECT_TRUE(bulkData.is_valid(elem_2));
            EXPECT_EQ(1, bulkData.parallel_owner_rank(elem_2));

            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_2));

            stk::mesh::Entity elem = elem_graph.get_connected_element(elem_2, 1);
            ASSERT_TRUE(elem_graph.is_connected_elem_locally_owned(elem_2, 1));
            EXPECT_EQ(3u, bulkData.identifier(elem));

            stk::mesh::EntityId connected_elem_global_id = elem_graph.get_entity_id_of_remote_element(elem_2, 0);
            ASSERT_FALSE(elem_graph.is_connected_elem_locally_owned(elem_2, 0));
            EXPECT_EQ(1u, connected_elem_global_id);

            stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            int side_from_elem2_to_elem3 = elem_graph.get_side_from_element1_to_locally_owned_element2(elem_2, elem_3);
            EXPECT_EQ(5, side_from_elem2_to_elem3);

            int side_from_elem2_to_elem1 = elem_graph.get_side_from_element1_to_remote_element2(elem_2, stk::mesh::EntityId(1));
            EXPECT_EQ(4, side_from_elem2_to_elem1);

            stk::mesh::impl::parallel_info &elem_2_to_1_p_info = elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1));
            int other_side_ord = elem_2_to_1_p_info.m_other_side_ord;
            EXPECT_EQ(5, other_side_ord);
            stk::mesh::EntityId chosen_face_id = 2;
            EXPECT_EQ(chosen_face_id, elem_2_to_1_p_info.m_chosen_side_id);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_3, stk::mesh::EntityId(2)), std::logic_error);
        }

        EXPECT_EQ(1u, elem_graph.num_parallel_edges());
    }
}

TEST(ElementGraph, test_change_entity_owner_2_procs_hex_mesh_with_aura)
{
    bool aura_on = true;
    change_entity_owner_hex_test_2_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_2_procs_hex_mesh_without_aura)
{
    bool aura_on = false;
    change_entity_owner_hex_test_2_procs(aura_on);
}

void change_entity_owner_then_death_hex_test_2_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);

    if(stk::parallel_machine_size(comm) == 2)
    {
        unsigned spatial_dim = 3;
        stk::mesh::MetaData meta(spatial_dim, stk::mesh::entity_rank_names());

        stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::mesh::Part& active = meta.declare_part("active", stk::topology::ELEMENT_RANK);
        stk::mesh::PartVector boundary_mesh_parts { &faces_part };

        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        stk::unit_test_util::put_mesh_into_part(bulkData, active);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(2, numLocallyOwnedElems);

        ElemElemGraphTester elem_graph(bulkData);

        // Create a vector of the elements to be moved
        std::vector <stk::mesh::Entity> elems_to_move;

        stk::mesh::EntityId elem_2_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

        if (proc == 0)
        {
            int side_from_elem2_to_elem3 = elem_graph.get_side_from_element1_to_remote_element2(elem_2, stk::mesh::EntityId(3));
            int side_from_elem2_to_elem1 = elem_graph.get_side_from_element1_to_locally_owned_element2(elem_2, bulkData.get_entity(stk::topology::ELEM_RANK,1));

            EXPECT_EQ(5, side_from_elem2_to_elem3);
            EXPECT_EQ(4, side_from_elem2_to_elem1);

            elems_to_move.push_back(elem_2);

            int other_proc = 1;
            for (unsigned i=0; i<elems_to_move.size(); i++)
            {
                EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
                EXPECT_EQ(0, bulkData.parallel_owner_rank(elems_to_move[i]));
                elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
            }
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move, &active);

        elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

        stk::mesh::EntityVector killedElements;
        std::vector<stk::mesh::PartVector> add_parts, remove_parts;
        if (proc == 1)
        {
            killedElements.push_back(elem_2);
            add_parts.push_back(stk::mesh::PartVector());
            remove_parts.push_back(stk::mesh::PartVector{&active});
        }

        bulkData.batch_change_entity_parts(killedElements, add_parts, remove_parts);
        boundary_mesh_parts.push_back(&active);

        process_killed_elements(bulkData, elem_graph, killedElements, active, boundary_mesh_parts, &boundary_mesh_parts);

        if (proc == 1)
        {
            EXPECT_TRUE(bulkData.is_valid(elem_2));
            EXPECT_EQ(1, bulkData.parallel_owner_rank(elem_2));
            EXPECT_FALSE(bulkData.bucket(elem_2).member(active));

            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_2));

            stk::mesh::Entity elem = elem_graph.get_connected_element(elem_2, 1);
            ASSERT_TRUE(elem_graph.is_connected_elem_locally_owned(elem_2, 1));
            EXPECT_EQ(3u, bulkData.identifier(elem));

            stk::mesh::EntityId connected_elem_global_id = elem_graph.get_entity_id_of_remote_element(elem_2, 0);
            ASSERT_FALSE(elem_graph.is_connected_elem_locally_owned(elem_2, 0));
            EXPECT_EQ(1u, connected_elem_global_id);

            stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            int side_from_elem2_to_elem3 = elem_graph.get_side_from_element1_to_locally_owned_element2(elem_2, elem_3);
            EXPECT_EQ(5, side_from_elem2_to_elem3);

            int side_from_elem2_to_elem1 = elem_graph.get_side_from_element1_to_remote_element2(elem_2, stk::mesh::EntityId(1));
            EXPECT_EQ(4, side_from_elem2_to_elem1);

            stk::mesh::impl::parallel_info &elem_2_to_1_p_info = elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1));
            int other_side_ord = elem_2_to_1_p_info.m_other_side_ord;
            EXPECT_EQ(5, other_side_ord);
        }
        if (proc == 0)
        {
            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            stk::mesh::impl::parallel_info &elem1_to_elem2_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));

            EXPECT_FALSE(elem1_to_elem2_info.m_in_body_to_be_skinned);

            EXPECT_EQ(4, elem1_to_elem2_info.m_other_side_ord);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(3)), std::logic_error);
       }

        EXPECT_EQ(1u, elem_graph.num_parallel_edges());
    }
}

TEST(ElementGraph, test_change_entity_owner_and_death_hex_mesh_2_procs_with_aura)
{
    bool aura_on = true;
    change_entity_owner_then_death_hex_test_2_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_and_death_hex_mesh_2_procs_without_aura)
{
    bool aura_on = false;
    change_entity_owner_then_death_hex_test_2_procs(aura_on);
}


void setup_hex_shell_hex_mesh(stk::mesh::BulkData& bulkData)
{
//
//                proc 0               proc 1           proc 2
//
//               block_1          |   block_2  |      block_3
//
//          3---------------7        7            7-------------11
//          /|             /|       /|           /|             /|
//         / |            / |      / |          / |            / |
//        /  |           /  |     /  |         /  |           /  |
//       4--------------8   |    8   |        8--------------12  |
//       |   |          |   |    |   |        |   |          |   |
//       |   |   1      |   |    | 2 |        |   |   3      |   |
//       |   |          |   |    |   |        |   |          |   |
//       |   2----------|---6    |   6        |   6----------|---10
//       |  /           |  /     |  /         |  /           |  /
//       | /            | /      | /          | /            | /
//       |/             |/       |/           |/             |/
//       1--------------5        5            5--------------9

    stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    unsigned spatial_dimension = 3;
    meta.initialize(spatial_dimension, stk::mesh::entity_rank_names());

    stk::mesh::Field<double>& field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
    stk::mesh::put_field(field, meta.universal_part());

    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::HEX_8);
    meta.commit();

    bulkData.modification_begin();

    stk::mesh::EntityIdVector elem1_nodes {1, 2, 3, 4, 5, 6, 7, 8};
    stk::mesh::EntityIdVector elem2_nodes {5, 6, 7, 8};
    stk::mesh::EntityIdVector elem3_nodes {5, 6, 7, 8, 9, 10, 11, 12};

    stk::mesh::EntityId elemId = 1;
    if (bulkData.parallel_rank() == 0) {
        stk::mesh::declare_element(bulkData, block_1, elemId, elem1_nodes);
        stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
        stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
        stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
        bulkData.add_node_sharing(node5, 1);
        bulkData.add_node_sharing(node6, 1);
        bulkData.add_node_sharing(node7, 1);
        bulkData.add_node_sharing(node8, 1);
        bulkData.add_node_sharing(node5, 2);
        bulkData.add_node_sharing(node6, 2);
        bulkData.add_node_sharing(node7, 2);
        bulkData.add_node_sharing(node8, 2);
    }
    else if (bulkData.parallel_rank() == 1) {
        elemId = 2;
        stk::mesh::declare_element(bulkData, block_2, elemId, elem2_nodes);
        stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
        stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
        stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
        bulkData.add_node_sharing(node5, 0);
        bulkData.add_node_sharing(node6, 0);
        bulkData.add_node_sharing(node7, 0);
        bulkData.add_node_sharing(node8, 0);
        bulkData.add_node_sharing(node5, 2);
        bulkData.add_node_sharing(node6, 2);
        bulkData.add_node_sharing(node7, 2);
        bulkData.add_node_sharing(node8, 2);
    }
    else if (bulkData.parallel_rank() == 2) {
        elemId = 3;
        stk::mesh::declare_element(bulkData, block_3, elemId, elem3_nodes);
        stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
        stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
        stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
        bulkData.add_node_sharing(node5, 0);
        bulkData.add_node_sharing(node6, 0);
        bulkData.add_node_sharing(node7, 0);
        bulkData.add_node_sharing(node8, 0);
        bulkData.add_node_sharing(node5, 1);
        bulkData.add_node_sharing(node6, 1);
        bulkData.add_node_sharing(node7, 1);
        bulkData.add_node_sharing(node8, 1);
    }

    bulkData.modification_end();
}

void change_entity_owner_hex_shell_hex_test_3_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);
    if(stk::parallel_machine_size(comm) == 3)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        setup_hex_shell_hex_mesh(bulkData);

        ElemElemGraphTester elem_graph(bulkData);

        const stk::mesh::Entity hex1   = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        const stk::mesh::Entity hex3   = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
        const stk::mesh::Entity shell2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);

        if (proc == 0) {
            // Connectivity for Hex Element 1
            EXPECT_EQ(1u, elem_graph.get_num_connected_elems(hex1));
            EXPECT_EQ(5,  elem_graph.get_side_id_to_connected_element(hex1, 0));
            EXPECT_EQ(2u, elem_graph.get_entity_id_of_remote_element(hex1, 0));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(hex1, 0));
            EXPECT_EQ(1u, elem_graph.num_edges());
            EXPECT_EQ(1u, elem_graph.num_parallel_edges());
        }
        else if (proc == 1) {
            // Connectivity for Shell Element 2
            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(shell2));
            EXPECT_EQ(0,  elem_graph.get_side_id_to_connected_element(shell2, 0));
            EXPECT_EQ(1,  elem_graph.get_side_id_to_connected_element(shell2, 1));
            EXPECT_EQ(3u, elem_graph.get_entity_id_of_remote_element(shell2, 0));
            EXPECT_EQ(1u, elem_graph.get_entity_id_of_remote_element(shell2, 1));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(shell2, 0));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(shell2, 1));
            EXPECT_EQ(2u, elem_graph.num_edges());
            EXPECT_EQ(2u, elem_graph.num_parallel_edges());
        }
        else if (proc == 2) {
            // Connectivity for Hex Element 3
            EXPECT_EQ(1u, elem_graph.get_num_connected_elems(hex3));
            EXPECT_EQ(4,  elem_graph.get_side_id_to_connected_element(hex3, 0));
            EXPECT_EQ(2u, elem_graph.get_entity_id_of_remote_element(hex3, 0));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(hex3, 0));
            EXPECT_EQ(1u, elem_graph.num_edges());
            EXPECT_EQ(1u, elem_graph.num_parallel_edges());
        }

        stk::mesh::EntityId elem_to_move_global_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_to_move_global_id);

        if (proc == 1)
        {
            int destination_proc = 2;
            elem_proc_pairs_to_move.push_back(std::make_pair(elem_to_move, destination_proc));
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElemsInMesh = counts[stk::topology::ELEM_RANK];

        size_t size_of_elem_graph = elem_graph.size();

        if (proc == 0)
        {
            EXPECT_EQ(1, numLocallyOwnedElemsInMesh);
            EXPECT_EQ(1u, size_of_elem_graph);

            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            stk::mesh::impl::parallel_info& elem1_to_elem2_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));
            EXPECT_EQ(2, elem1_to_elem2_info.m_other_proc);
            EXPECT_EQ(1u, elem_graph.num_edges());
            EXPECT_EQ(1u, elem_graph.num_parallel_edges());
        }
        if (proc == 1)
        {
            EXPECT_EQ(0, numLocallyOwnedElemsInMesh);
            EXPECT_EQ(0u, size_of_elem_graph);

            stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1)), std::logic_error);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(3)), std::logic_error);
            EXPECT_EQ(0u, elem_graph.num_edges());
            EXPECT_EQ(0u, elem_graph.num_parallel_edges());
        }
        if (proc == 2)
        {
            EXPECT_EQ(2, numLocallyOwnedElemsInMesh);
            EXPECT_EQ(2u, size_of_elem_graph);

            stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
            stk::mesh::impl::parallel_info& elem2_to_elem1_info = elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1));
            EXPECT_EQ(0, elem2_to_elem1_info.m_other_proc);
            EXPECT_EQ(3u, elem_graph.num_edges());
            EXPECT_EQ(1u, elem_graph.num_parallel_edges());
        }
    }
}

TEST(ElementGraph, test_change_entity_owner_3_procs_hex_shell_hex_mesh_with_aura)
{
    bool aura_on = true;
    change_entity_owner_hex_shell_hex_test_3_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_3_procs_hex_shell_hex_mesh_without_aura)
{
    bool aura_on = false;
    change_entity_owner_hex_shell_hex_test_3_procs(aura_on);
}

void change_entity_owner_hex_test_4_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);
    std::vector<double> wall_times;
    wall_times.reserve(10);
    std::vector<std::string> msgs;
    msgs.reserve(10);

    std::vector<size_t> mem_usage;

    wall_times.push_back(stk::wall_time());
    msgs.push_back("program-start");
    mem_usage.push_back(stk::get_memory_usage_now());

    if(stk::parallel_machine_size(comm) == 4)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(1, numLocallyOwnedElems);

        ElemElemGraphTester elem_graph(bulkData);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        // Create a vector of the elements to be moved
        std::vector <stk::mesh::Entity> elems_to_move;

        stk::mesh::EntityId elem_global_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_global_id);
        if (proc == 1)
        {
            elems_to_move.push_back(elem_to_move);

            int other_proc = 2;
            for (unsigned i=0; i<elems_to_move.size(); i++)
            {
                EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
                EXPECT_EQ(1, bulkData.parallel_owner_rank(elems_to_move[i]));
                elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
            }
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move);

        elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_global_id);

        if (proc == 2)
        {
            EXPECT_TRUE(bulkData.is_valid(elem_to_move));
            EXPECT_EQ(2, bulkData.parallel_owner_rank(elem_to_move));

            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_to_move));

            stk::mesh::Entity elem = elem_graph.get_connected_element(elem_to_move, 1);
            ASSERT_TRUE(elem_graph.is_connected_elem_locally_owned(elem_to_move, 1));
            EXPECT_EQ(3u, bulkData.identifier(elem));

            stk::mesh::EntityId connected_elem_global_id = elem_graph.get_entity_id_of_remote_element(elem_to_move, 0);
            ASSERT_FALSE(elem_graph.is_connected_elem_locally_owned(elem_to_move, 0));
            EXPECT_EQ(1u, connected_elem_global_id);

            stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_3, stk::mesh::EntityId(2)), std::logic_error);

            EXPECT_EQ(2u, elem_graph.num_parallel_edges());
        }
        else if (proc == 0)
        {
            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            stk::mesh::impl::parallel_info &elem_1_to_2_p_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));
            EXPECT_EQ(2, elem_1_to_2_p_info.m_other_proc);
            EXPECT_EQ(1u, elem_graph.num_parallel_edges());
        }
        else if (proc == 1)
        {
            EXPECT_EQ(0u, elem_graph.size());
            EXPECT_EQ(0u, elem_graph.num_edges());
            EXPECT_EQ(0u, elem_graph.num_parallel_edges());
        }
        else if (proc == 3)
        {
            EXPECT_EQ(1u, elem_graph.num_parallel_edges());
        }
    }
}

TEST(ElementGraph, test_change_entity_owner_4_procs_hex_mesh_with_aura)
{
    bool aura_on = true;
    change_entity_owner_hex_test_4_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_4_procs_hex_mesh_without_aura)
{
    bool aura_on = false;
    change_entity_owner_hex_test_4_procs(aura_on);
}

} //namespace
