#include <gtest/gtest.h>                // for TEST
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for BulkData
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_util/parallel/ParallelComm.hpp>


namespace
{

void test_that_ids_are_unique(stk::mesh::BulkData &bulkData, stk::topology::rank_t rank, std::vector<stk::mesh::EntityId>& requestedIds)
{
    std::vector<stk::mesh::EntityId> ids_in_use;
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(rank, bulkData.mesh_meta_data().locally_owned_part());
    for (size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        for (size_t j=0;j<bucket.size();++j)
        {
            ids_in_use.push_back(bulkData.identifier(bucket[j]));
        }
    }

    // ========

    std::sort(requestedIds.begin(), requestedIds.end());
    std::vector<stk::mesh::EntityId>::iterator iter1 = std::unique(requestedIds.begin(), requestedIds.end());
    ThrowRequireMsg(iter1 == requestedIds.end(), "Oh no! " << __FILE__ << __LINE__);

    stk::CommAll comm(bulkData.parallel());

    for(int phase = 0; phase < 2; ++phase)
    {
        for(int i = 0; i < bulkData.parallel_size(); ++i)
        {
            if(i != bulkData.parallel_rank())
            {
                for(size_t j = 0; j < requestedIds.size(); ++j)
                {
                    bool is_id_unique = std::binary_search(ids_in_use.begin(), ids_in_use.end(), requestedIds[j]);
                    ThrowRequireMsg(is_id_unique == false, "Oh no! " <<  __FILE__<< __LINE__);
                    comm.send_buffer(i).pack<uint64_t>(requestedIds[j]);
                }
            }
        }

        if(phase == 0)
        {
            comm.allocate_buffers(bulkData.parallel_size() / 4);
        }
        else
        {
            comm.communicate();
        }
    }

    for(int i = 0; i < bulkData.parallel_size(); ++i)
    {
        if(i != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(i).remaining())
            {
                uint64_t key;
                comm.recv_buffer(i).unpack<uint64_t>(key);
                bool is_other_procs_id_on_this_proc = std::binary_search(requestedIds.begin(), requestedIds.end(), key);
                ThrowRequireMsg(is_other_procs_id_on_this_proc == false, "Oh no! " << __FILE__<< __LINE__);
                bool is_id_already_in_use = std::binary_search(ids_in_use.begin(), ids_in_use.end(), key);
                ThrowRequireMsg(is_id_already_in_use == false, "Oh no! " << __FILE__ << __LINE__);
            }
        }
    }
}

//BEGIN_TEST_1
TEST(StkMeshHowTo, use_generate_new_ids)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    int num_procs = -1;
    MPI_Comm_size(communicator, &num_procs);
    std::ostringstream os;
    os << "generated:1x1x" << num_procs;
    const std::string generatedMeshSpecification = os.str();

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    // Given a mesh, request 10 unique node ids

    std::vector<stk::mesh::EntityId> requestedIds;
    unsigned numRequested = 10;

    stkMeshBulkData.generate_new_ids(stk::topology::NODE_RANK, numRequested, requestedIds);

    test_that_ids_are_unique(stkMeshBulkData, stk::topology::NODE_RANK, requestedIds);
}
//END_TEST_1

}

