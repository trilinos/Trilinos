#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>

namespace
{

class CommInfoObserver : public stk::mesh::ModificationObserver
{
public:
    CommInfoObserver()
    : commInfoWasChangedByRank(stk::topology::NUM_RANKS, false)
    {
    }

    virtual void local_entity_comm_info_changed_notification(stk::mesh::EntityRank rank)
    {
        commInfoWasChangedByRank[rank] = true;
    }

    bool was_comm_info_changed(stk::mesh::EntityRank rank) const
    {
        return commInfoWasChangedByRank[rank];
    }

    void reset_comm_info_status()
    {
        for(size_t i=0; i < commInfoWasChangedByRank.size(); i++)
        {
            commInfoWasChangedByRank[i] = false;
        }
    }
private:
    std::vector<bool> commInfoWasChangedByRank;
};

void ghost_element1_to_proc1(stk::mesh::BulkData &bulk, stk::mesh::Ghosting &ghost)
{
    stk::mesh::EntityProcVec entitiesToGhost;
    if(bulk.parallel_rank() == 0)
    {
        stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
        int destProc = 1;
        entitiesToGhost.push_back(stk::mesh::EntityProc(element1, destProc));
    }
    bulk.modification_begin();
    bulk.change_ghosting(ghost, entitiesToGhost);
    bulk.modification_end();
}

void expect_comm_info_changed_for_nodes_and_elements_only(const CommInfoObserver &observer)
{
    EXPECT_TRUE(observer.was_comm_info_changed(stk::topology::NODE_RANK));
    EXPECT_TRUE(observer.was_comm_info_changed(stk::topology::ELEM_RANK));

    EXPECT_TRUE(!observer.was_comm_info_changed(stk::topology::EDGE_RANK));
    EXPECT_TRUE(!observer.was_comm_info_changed(stk::topology::FACE_RANK));
    EXPECT_TRUE(!observer.was_comm_info_changed(stk::topology::CONSTRAINT_RANK));
}

void expect_comm_info_is_reset(CommInfoObserver &observer)
{
    for(stk::mesh::EntityRank rank = stk::topology::BEGIN_RANK; rank < stk::topology::END_RANK; rank++)
    {
        EXPECT_TRUE(!observer.was_comm_info_changed(rank));
    }
}

void stop_ghosting_of_element1_to_proc1(stk::mesh::BulkData &bulk, stk::mesh::Ghosting &ghost)
{
    std::vector<stk::mesh::EntityKey> entityKeysToStopGhosting;
    if(bulk.parallel_rank() == 1)
    {
        entityKeysToStopGhosting.push_back(stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1));
        for(int nodeId=1; nodeId <= 8; nodeId++)
        {
            entityKeysToStopGhosting.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, nodeId));
        }
    }
    bulk.modification_begin();
    bulk.change_ghosting(ghost, {}, entityKeysToStopGhosting);
    bulk.modification_end();
}

TEST(CommInfoObserver, addAndRemoveFromGhosting)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulk, comm);

        CommInfoObserver observer;
        bulk.register_observer(&observer);

        bulk.modification_begin();
        stk::mesh::Ghosting &ghost = bulk.create_ghosting("Clyde");
        bulk.modification_end();

        ghost_element1_to_proc1(bulk, ghost);
        expect_comm_info_changed_for_nodes_and_elements_only(observer);

        observer.reset_comm_info_status();
        expect_comm_info_is_reset(observer);

        stop_ghosting_of_element1_to_proc1(bulk, ghost);

        expect_comm_info_changed_for_nodes_and_elements_only(observer);
    }
}

TEST(CommInfoObserver, destroyGhosting)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulk, comm);

        CommInfoObserver observer;
        bulk.register_observer(&observer);

        bulk.modification_begin();
        stk::mesh::Ghosting &ghost = bulk.create_ghosting("Clyde");
        bulk.modification_end();

        ghost_element1_to_proc1(bulk, ghost);
        expect_comm_info_changed_for_nodes_and_elements_only(observer);

        observer.reset_comm_info_status();
        expect_comm_info_is_reset(observer);

        bulk.modification_begin();
        bulk.destroy_all_ghosting();
        bulk.modification_end();

        expect_comm_info_changed_for_nodes_and_elements_only(observer);
    }
}

}
