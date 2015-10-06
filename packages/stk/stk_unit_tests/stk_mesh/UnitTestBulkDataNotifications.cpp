#include <iostream>                     // for ostringstream, etc
#include <algorithm>                    // for sort
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>
#include <unit_tests/BulkDataTester.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/ModificationObserver.hpp>

namespace {

void reduce_max(std::vector<size_t> &localVector, stk::ParallelMachine comm)
{
    std::vector<size_t> globalVector(localVector.size());
    stk::all_reduce_max(comm, localVector.data(), globalVector.data(), globalVector.size());
    localVector.swap(globalVector);
}

class TestListener : public stk::mesh::ModificationObserver {
public:
    TestListener(stk::ParallelMachine communicator)
      : stk::mesh::ModificationObserver(),
        comm(communicator),
        buckets_changed(5, 0)
    {}

    virtual void local_buckets_changed_notification(stk::mesh::EntityRank rank)
    {
        buckets_changed[rank]++;
    }

    virtual void finished_modification_end_notification()
    {
        reduce_max(buckets_changed, comm);
    }

    size_t get_buckets_changed(stk::mesh::EntityRank rank) const
    { return buckets_changed[rank]; }

private:
    stk::ParallelMachine comm;
    std::vector<size_t> buckets_changed;
};

}

TEST(BulkDataNotifications, test_listener_buckets_changed)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(comm);
    if (numProcs==2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::unit_test::BulkDataTester mesh(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);

        const std::string generatedMeshSpec = "generated:1x1x2";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, mesh, comm);

        TestListener listener(comm);
        mesh.register_observer(&listener);

        int procId = stk::parallel_machine_rank(comm);

        stk::mesh::Part& newPart = meta.declare_part("new part", stk::topology::NODE_RANK);
        stk::mesh::Entity node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);

        if (procId == 0) {
            EXPECT_TRUE(mesh.is_valid(node1));
        }
        else {
            EXPECT_FALSE(mesh.is_valid(node1));
        }

        mesh.modification_begin();
        if (procId == 0) {
            mesh.change_entity_parts(node1, {&newPart});
        }
        mesh.modification_end();

        EXPECT_EQ(1u, listener.get_buckets_changed(stk::topology::NODE_RANK));
    }
}

