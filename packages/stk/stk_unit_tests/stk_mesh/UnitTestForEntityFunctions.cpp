#include <gtest/gtest.h>
#include <string>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_io/StkMeshIoBroker.hpp"

namespace
{

struct AlgorithmPerEntity
{
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity) = 0;
};

struct AlgorithmPerCommunicatedEntity
{
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity, const std::vector<int> &comm_procs) = 0;
};

class BulkDataForEntityTester : public stk::mesh::BulkData
{
public:
    BulkDataForEntityTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }
    virtual ~BulkDataForEntityTester()
    {
    }

    void for_each_entity_run(stk::topology::rank_t rank, AlgorithmPerEntity &functor)
    {
        const stk::mesh::BucketVector & buckets = this->buckets(rank);
        const size_t numBuckets = buckets.size();
        for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
        {
            stk::mesh::Bucket & bucket = *buckets[iBucket];
            const unsigned numEntitiesInBucket = bucket.size();
            for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
            {
                stk::mesh::Entity entity = bucket[iEntity];
                functor.run_on_entity(*this, entity);
            }
        }
    }

    void for_each_node_run(AlgorithmPerEntity &functor)
    {
        for_each_entity_run(stk::topology::NODE_RANK, functor);
    }

    void for_communicated_entities_run(AlgorithmPerCommunicatedEntity &algorithm) const
    {
        std::vector<int> comm_procs;
        stk::mesh::EntityCommListInfoVector::const_iterator iter = internal_comm_list().begin(),
                                                            iend = internal_comm_list().end();
        for(; iter != iend; ++iter)
        {
            stk::mesh::Entity entity = iter->entity;

            comm_procs.clear();
            const stk::mesh::EntityCommInfoVector& infovec = iter->entity_comm->comm_map;
            stk::mesh::PairIterEntityComm ec(infovec.begin(), infovec.end());
            for(; !ec.empty(); ++ec)
            {
                comm_procs.push_back(ec->proc);
            }

            algorithm.run_on_entity(*this, entity, comm_procs);
        }
    }
};

void generateMesh(const std::string &generatedMeshSpec, stk::mesh::BulkData &bulkData, MPI_Comm communicator)
{
    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.set_bulk_data(bulkData);
    exodusFileReader.add_mesh_database(generatedMeshSpec, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();
}

struct CountNumNodesAlgorithm : public AlgorithmPerEntity
{
    CountNumNodesAlgorithm(int &numNodes) :
            mNumNodes(numNodes)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
    {
        mNumNodes++;
    }
    int &mNumNodes;
};

TEST(ForEntityFunction, test_for_each_node_run)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        generateMesh(generatedMeshSpec, bulkData, communicator);

        int numNodes = 0;
        CountNumNodesAlgorithm countNumNodesAlgorithm(numNodes);
        bulkData.for_each_node_run(countNumNodesAlgorithm);

        EXPECT_EQ(16, numNodes);
    }
}

struct CountNumCommunicatedEntitiesAlgorithm : public AlgorithmPerCommunicatedEntity
{
    CountNumCommunicatedEntitiesAlgorithm(size_t &numCommunicatedEntities) :
            mNumCommunicatedEntities(numCommunicatedEntities)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity, const std::vector<int> &comm_procs)
    {
        mNumCommunicatedEntities++;
    }
    size_t &mNumCommunicatedEntities;
};

struct FillCommProcsAlgorithm : public AlgorithmPerCommunicatedEntity
{
    FillCommProcsAlgorithm(std::vector<int> &commProcs) :
            mCommProcs(commProcs)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity, const std::vector<int> &comm_procs)
    {
        for(size_t i = 0; i < comm_procs.size(); i++)
        {
            mCommProcs.push_back(comm_procs[i]);
        }
    }
    std::vector<int> &mCommProcs;
};

TEST(ForEntityFunction, test_for_communicated_entities_run)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        generateMesh(generatedMeshSpec, bulkData, communicator);

        size_t numCommunicatedEntities = 0;
        CountNumCommunicatedEntitiesAlgorithm countNumCommunicatedEntitiesAlgorithm(numCommunicatedEntities);
        bulkData.for_communicated_entities_run(countNumCommunicatedEntitiesAlgorithm);

        size_t numSharedNodes = 4;
        size_t numSendAuraNodes = 4;
        size_t numRecvAuraNodes = 4;
        size_t numSendAuraElems = 1;
        size_t numRecvAuraElems = 1;
        size_t numExpectedCommunicatedEntities = numSharedNodes + numSendAuraNodes + numRecvAuraNodes + numSendAuraElems
                                              + numRecvAuraElems;
        EXPECT_EQ(numExpectedCommunicatedEntities, numCommunicatedEntities);

        std::vector<int> commProcs;
        FillCommProcsAlgorithm fillCommProcs(commProcs);
        bulkData.for_communicated_entities_run(fillCommProcs);

        ASSERT_EQ(numExpectedCommunicatedEntities, commProcs.size());
        int otherProc = 1 - stk::parallel_machine_rank(communicator);
        for(size_t i = 0; i < commProcs.size(); i++)
        {
            EXPECT_EQ(otherProc, commProcs[i]);
        }
    }
}

}
