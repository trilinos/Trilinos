#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <mpi.h>
#include <stk_util/parallel/DistributedIndex.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

namespace {

class StkMeshCreator
{

public:
    StkMeshCreator(const std::string& generatedMeshSpec, MPI_Comm communicator) : m_exodusFilename("junk.par"),
    m_stkMeshMetaData(NULL), m_stkMeshBulkData(NULL)
    {
        this->writeExodusFile(generatedMeshSpec, m_exodusFilename, communicator);

        const int spatialDim = 3;
        m_stkMeshMetaData = new stk::mesh::MetaData(spatialDim);
        m_stkMeshBulkData = new stk::mesh::BulkData(*m_stkMeshMetaData, communicator);

        readExodusFileIntoStkMesh(m_exodusFilename, *m_stkMeshMetaData, *m_stkMeshBulkData, communicator);
    }

    ~StkMeshCreator()
    {
        delete m_stkMeshBulkData;
        delete m_stkMeshMetaData;
//        unlink(m_exodusFilename.c_str());
    }

    stk::mesh::MetaData* getMetaData() { return m_stkMeshMetaData; }
    stk::mesh::BulkData* getBulkData() { return m_stkMeshBulkData; }

private:

    void writeExodusFile(const std::string &generatedMeshSpecification, const std::string& exodusFileName, MPI_Comm communicator)
    {
        stk::io::StkMeshIoBroker stkIo(communicator);

        size_t index = stkIo.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
        stkIo.set_active_mesh(index);
        stkIo.create_input_mesh();
        stkIo.populate_bulk_data();

        size_t fh = stkIo.create_output_mesh(exodusFileName, stk::io::WRITE_RESULTS);
        stkIo.write_output_mesh(fh);
    }

    void readExodusFileIntoStkMesh(const std::string& exodusFileName, stk::mesh::MetaData &stkMeshMetaData, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm communicator)
    {
        stk::io::StkMeshIoBroker exodusFileReader(communicator);

        // Inform STK IO which STK Mesh objects to populate later
        exodusFileReader.set_bulk_data(stkMeshBulkData);

        exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

        // Populate the MetaData which has the descriptions of the Parts and Fields.
        exodusFileReader.create_input_mesh();

        // Populate entities in STK Mesh from Exodus file
        exodusFileReader.populate_bulk_data();
    }

private:
    std::string m_exodusFilename;
    stk::mesh::MetaData *m_stkMeshMetaData;
    stk::mesh::BulkData *m_stkMeshBulkData;
};

void updateDistributedIndexUsingStkMesh(stk::mesh::BulkData &stkMeshBulkData, const int myProc, stk::parallel::DistributedIndex& distributedIndex)
{
    stk::parallel::DistributedIndex::KeyTypeVector local_created_or_modified; // only store locally owned/shared entities

    // Iterate over all entities known to this process, putting
    // modified shared/owned entities in local_created_or_modified.
    size_t num_created_or_modified = 0;

    stk::mesh::impl::EntityRepository &m_entity_repo = stkMeshBulkData.get_entity_repository();

    for(stk::mesh::impl::EntityRepository::const_iterator i = m_entity_repo.begin(); i != m_entity_repo.end(); ++i)
    {
        stk::mesh::Entity entity = i->second;

        if(stk::mesh::in_owned_closure(stkMeshBulkData, entity, myProc))
        {
            // Has been changed and is in owned closure, may be shared
            ++num_created_or_modified;
        }
    }

    local_created_or_modified.reserve(num_created_or_modified);

    for(stk::mesh::impl::EntityRepository::const_iterator i = m_entity_repo.begin(); i != m_entity_repo.end(); ++i)
    {
        stk::mesh::Entity entity = i->second;

        //if(stkMeshBulkData.state(entity) != stk::mesh::Unchanged && stk::mesh::in_owned_closure(stkMeshBulkData, entity, myProc))
        if(stk::mesh::in_owned_closure(stkMeshBulkData, entity, myProc))
        {
            // Has been changed and is in owned closure, may be shared
            local_created_or_modified.push_back(stkMeshBulkData.entity_key(entity));
        }
    }

    // Update distributed index. Note that the DistributedIndex only
    // tracks ownership and sharing information.
    stk::parallel::DistributedIndex::KeyTypeVector::const_iterator begin = local_created_or_modified.begin();
    stk::parallel::DistributedIndex::KeyTypeVector::const_iterator end = local_created_or_modified.end();
    distributedIndex.update_keys(begin, end);
}

TEST( UnderstandingDistributedIndex, WithoutStkMeshBulkData)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int procCount = -1;
    int myProc = -1;
    MPI_Comm_size(communicator, &procCount);
    MPI_Comm_rank(communicator, &myProc);

    int outputProcId = 1;

    if(procCount >= 1)
    {
        const std::string generatedMeshSpec = "generated:2x2x2|sideset:xXyYzZ|nodeset:xXyYzZ";
        StkMeshCreator stkMesh(generatedMeshSpec, communicator);

        stk::mesh::MetaData &stkMeshMetaData = *stkMesh.getMetaData();
        stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();

        stk::parallel::DistributedIndex::KeySpanVector spans = stk::mesh::convert_entity_keys_to_spans(stkMeshMetaData);
        stk::parallel::DistributedIndex distributedIndex(communicator, spans);

        updateDistributedIndexUsingStkMesh(stkMeshBulkData, myProc, distributedIndex);

        const unsigned rankCount = stkMeshMetaData.entity_rank_count();
        EXPECT_EQ(4u, rankCount);

        for(size_t i = 0; i < spans.size(); i++)
        {
            stk::parallel::DistributedIndex::KeyType keyMin = spans[i].first;
            stk::parallel::DistributedIndex::KeyType keyMax = spans[i].second;
            stk::mesh::EntityKey key_1(static_cast<stk::mesh::EntityKey::entity_key_t>((keyMin)));
            stk::mesh::EntityKey key_2(static_cast<stk::mesh::EntityKey::entity_key_t>((keyMax)));

            stk::mesh::EntityRank eRank = static_cast<stk::mesh::EntityRank>(i);

            EXPECT_EQ(eRank, key_1.rank());
            EXPECT_EQ(eRank, key_2.rank());
        }

        MPI_Barrier(communicator);

        ///////////////////////////////////////////////////////////////////////////////

        std::vector<stk::parallel::DistributedIndex::KeyTypeVector> requested_key_types;
        std::vector<size_t> requests(rankCount, 0);
        if ( myProc == 0 )
        {
            requests[0] = 4;
            requests[1] = 5;
            requests[2] = 6;
            requests[3] = 7;
        }

        size_t totalCount = 0;
        for(size_t i = 0; i < requests.size(); i++)
        {
            totalCount += requests[i];
        }

        distributedIndex.generate_new_keys(requests, requested_key_types);

        for(size_t i = 0; i < requested_key_types.size(); i++)
        {
            stk::mesh::EntityRank eRank = static_cast<stk::mesh::EntityRank>(i);
            if(myProc == outputProcId)
            {
                for(size_t j = 0; j < requested_key_types[i].size(); j++)
                {
                    stk::mesh::EntityKey keyType(static_cast<stk::mesh::EntityKey::entity_key_t>(requested_key_types[i][j]));
                    std::cerr << "Key " << i << "," << j << " = " << keyType.id() << " with key type of " << requested_key_types[i][j] << std::endl;
                    EXPECT_EQ(eRank, keyType.rank());
                }
            }
        }

        MPI_Barrier(communicator);

        stk::parallel::DistributedIndex::KeyProcVector keys = distributedIndex.getKeys();
        for (int procIndex=0;procIndex<procCount;procIndex++)
        {
            if(myProc == procIndex)
            {
                std::cerr << "================ OUTPUT FOR PROC " << myProc << " ===================" << std::endl;
                for(size_t i = 0; i < keys.size(); i++)
                {
                    stk::mesh::EntityKey keyType(static_cast<stk::mesh::EntityKey::entity_key_t>(keys[i].first));
                    std::cerr << "for i = " << i << " key is " << keys[i].first << " with proc rank " << keys[i].second <<
                            " with a entity id of " << keyType << std::endl;
                }
            }
            MPI_Barrier(communicator);
        }
    }
}

TEST( UnderstandingDistributedIndex, ViaStkMeshBulkData)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int procCount = -1;
    int myProc = -1;
    MPI_Comm_size(communicator, &procCount);
    MPI_Comm_rank(communicator, &myProc);

    int outputProcId = 1;

    if(procCount >= 1)
    {
        const std::string generatedMeshSpec = "generated:3x3x3|sideset:xXyYzZ|nodeset:xXyYzZ";
        StkMeshCreator stkMesh(generatedMeshSpec, communicator);

        stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();

        MPI_Barrier(communicator);

        const int rankCount = 4; // nodes, edges, faces, and elements
        std::vector<size_t> requests(rankCount, 0);
        if ( myProc == 0 )
        {
            requests[0] = 4;
            requests[1] = 5;
            requests[2] = 6;
            requests[3] = 7;
        }

        size_t totalCount = 0;
        for(size_t i = 0; i < requests.size(); i++)
        {
            totalCount += requests[i];
        }

        std::vector<stk::mesh::Entity> requested_entities;

        stkMeshBulkData.modification_begin();
        stkMeshBulkData.generate_new_entities(requests, requested_entities);
        stkMeshBulkData.modification_end();

        EXPECT_EQ(totalCount, requested_entities.size());

        for(size_t i = 0; i < requested_entities.size(); i++)
        {
            if(myProc == outputProcId)
            {
                std::cerr << "Entity " << i << " has rank " << stkMeshBulkData.entity_rank(requested_entities[i])
                   << " and id " << stkMeshBulkData.identifier(requested_entities[i]) << std::endl;
            }
        }

        MPI_Barrier(communicator);
    }
}

TEST(UnderstandingMesh, TestSharedAndGhostedAndOwnedEntitiesWithoutAnyModification)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int procCount = -1;
    int myProc = -1;
    MPI_Comm_size(communicator, &procCount);
    MPI_Comm_rank(communicator, &myProc);

//    int outputProcId = 1;

    if(procCount == 2)
    {
        const std::string generatedMeshSpec = "generated:2x2x2|sideset:xXyYzZ|nodeset:xXyYzZ";
        StkMeshCreator stkMesh(generatedMeshSpec, communicator);

//        stk::mesh::MetaData &stkMeshMetaData = *stkMesh.getMetaData();
        stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();

        int sharedNodeIds[] = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

        int otherProcId = 1;
        if ( myProc == 1 )
        {
            otherProcId = 0;
        }

        size_t numSharedNodes=9;
        for (size_t i=0;i<numSharedNodes;i++)
        {
            stk::mesh::Entity sharedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, sharedNodeIds[i]);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(sharedNode);
            EXPECT_TRUE(bucket.shared());
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm_sharing(stkMeshBulkData.entity_key(sharedNode));
            EXPECT_TRUE(commStuff.size() == 1);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t sharedButNotGhostedId = 0;
            EXPECT_TRUE((*commStuff.first).ghost_id == sharedButNotGhostedId);
        }

        int ghostedNodeIds_onProc1[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        int ghostedNodeIds_onProc0[] = { 19, 20, 21, 22, 23, 24, 25, 26, 27 };

        int *ghostedNodeIds=ghostedNodeIds_onProc1;
        if ( myProc == 0 )
        {
            ghostedNodeIds = ghostedNodeIds_onProc0;
        }

        size_t numGhostedNodes=9;
        for (size_t i=0;i<numGhostedNodes;i++)
        {
            stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, ghostedNodeIds[i]);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
            bool isGhosted = !bucket.shared() && !bucket.owned();
            EXPECT_TRUE(isGhosted);
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm(stkMeshBulkData.entity_key(ghostedNode));
            size_t numProcsOtherThanMeOnWhichEntityIsAGhost = 1;
            EXPECT_TRUE(commStuff.size() == numProcsOtherThanMeOnWhichEntityIsAGhost);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t auraGhostedId = 1; // aura is not custom ghosted (automatically generated one layer of ghosting)
            EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostedId);
        }

        int elementsOnProc0[] = { 1, 2, 3, 4 }; // ghosted on proc 1
        int elementsOnProc1[] = { 5, 6, 7, 8 }; // ghosted on proc 0
        int *ghostedElements=elementsOnProc0;
        if ( myProc == 0 )
        {
            ghostedElements = elementsOnProc1;
        }

        size_t numGhostedElements=4;
        for (size_t i=0;i<numGhostedElements;i++)
        {
            stk::mesh::Entity ghostedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElements[i]);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedElement);
            bool isGhosted = !bucket.shared() && !bucket.owned();
            EXPECT_TRUE(isGhosted);
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm(stkMeshBulkData.entity_key(ghostedElement));
            size_t numProcsOtherThanMeOnWhichEntityIsAGhost = 1;
            EXPECT_TRUE(commStuff.size() == numProcsOtherThanMeOnWhichEntityIsAGhost);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t auraGhostedId = 1; // aura is not custom ghosted (automatically generated one layer of ghosting)
            EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostedId);
        }

        int *ownedElements = elementsOnProc1;
        if ( myProc == 0 )
        {
            ownedElements = elementsOnProc0;
        }

        size_t numOwnedElements=4;
        for (size_t i=0;i<numOwnedElements;i++)
        {
            stk::mesh::Entity ownedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ownedElements[i]);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ownedElement);
            EXPECT_TRUE(bucket.owned());
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm(stkMeshBulkData.entity_key(ownedElement));
            size_t numProcsOtherThanMeOnWhichEntityIsAGhost = 1;
            EXPECT_TRUE(commStuff.size() == numProcsOtherThanMeOnWhichEntityIsAGhost);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t auraGhostedId = 1; // aura is not custom ghosted (automatically generated one layer of ghosting)
            EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostedId);
        }
    }
}

TEST(UnderstandingMesh, GhostAnElement)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int procCount = -1;
    int myProc = -1;
    MPI_Comm_size(communicator, &procCount);
    MPI_Comm_rank(communicator, &myProc);

    if(procCount == 2)
    {
        const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
        StkMeshCreator stkMesh(generatedMeshSpec, communicator);

        stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();


        int sharedNodeIds[] = { 19, 20, 21, 22, 23, 24, 25, 26, 27 };

        int otherProcId = 1;
        if ( myProc == 1 )
        {
            otherProcId = 0;
        }

        size_t numSharedNodes=9;
        for (size_t i=0;i<numSharedNodes;i++)
        {
            stk::mesh::Entity sharedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, sharedNodeIds[i]);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(sharedNode);
            EXPECT_TRUE(bucket.shared());
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm_sharing(stkMeshBulkData.entity_key(sharedNode));
            EXPECT_TRUE(commStuff.size() == 1);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t sharedButNotGhostedId = 0;
            EXPECT_TRUE((*commStuff.first).ghost_id == sharedButNotGhostedId);
        }

        // elements 1-8 are proc 0
        // elements 9-16 are proc 1
        // 5-8 are ghosted on proc 1
        // 9-12 are ghosted on proc 0

        // proc 1 needs to ghost element 13 to proc 0
        int ghostedElementId = 13;
        int fromProc = 1;
        int toProc = 0;

        stkMeshBulkData.modification_begin();
        stk::mesh::Ghosting &ghosting = stkMeshBulkData.create_ghosting("Ghost Element 13");
        std::vector< std::pair<stk::mesh::Entity, int> > ghostingStruct;
        if ( myProc == fromProc )
        {
            stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElementId);
            ghostingStruct.push_back(std::make_pair(element, toProc));
        }
        stkMeshBulkData.change_ghosting(ghosting, ghostingStruct);
        stkMeshBulkData.modification_end();

        if ( myProc == toProc )
        {
            stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElementId);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(element);
            bool isGhosted = !bucket.shared() && !bucket.owned();
            EXPECT_TRUE(isGhosted);
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm(stkMeshBulkData.entity_key(element));
            EXPECT_TRUE(commStuff.size() == 1);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t customGhosting = ghosting.ordinal();
            EXPECT_TRUE((*commStuff.first).ghost_id == customGhosting);
        }
    }
}

TEST(UnderstandingMesh, KillAGhostedElement)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int procCount = -1;
    int myProc = -1;
    MPI_Comm_size(communicator, &procCount);
    MPI_Comm_rank(communicator, &myProc);

    if(procCount == 2)
    {
        const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
        StkMeshCreator stkMesh(generatedMeshSpec, communicator);

        stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();


        int sharedNodeIds[] = { 19, 20, 21, 22, 23, 24, 25, 26, 27 };

        int otherProcId = 1;
        if ( myProc == 1 )
        {
            otherProcId = 0;
        }

        size_t numSharedNodes=9;
        for (size_t i=0;i<numSharedNodes;i++)
        {
            stk::mesh::Entity sharedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, sharedNodeIds[i]);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(sharedNode);
            EXPECT_TRUE(bucket.shared());
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm_sharing(stkMeshBulkData.entity_key(sharedNode));
            EXPECT_TRUE(commStuff.size() == 1);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t sharedButNotGhostedId = 0;
            EXPECT_TRUE((*commStuff.first).ghost_id == sharedButNotGhostedId);
        }

        // elements 1-8 are proc 0
        // elements 9-16 are proc 1
        // 5-8 are ghosted on proc 1
        // 9-12 are ghosted on proc 0

        // proc 1 needs to ghost element 13 to proc 0
        int elementIdToKill = 8;
        int owningProc = 0;
        int ghostedToProc = 1;

        stk::mesh::Entity elementToKill = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdToKill);

        stkMeshBulkData.modification_begin();
        if ( myProc == ghostedToProc )
        {
            stkMeshBulkData.destroy_entity(elementToKill);
        }
        stkMeshBulkData.modification_end();

        if ( myProc == ghostedToProc )
        {
            stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdToKill);
            stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(element);
            bool isGhosted = !bucket.shared() && !bucket.owned();
            EXPECT_TRUE(isGhosted);
            stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.entity_comm(stkMeshBulkData.entity_key(element));
            EXPECT_TRUE(commStuff.size() == 1);
            EXPECT_TRUE((*commStuff.first).proc == otherProcId);
            size_t auraGhostingId = 1;
            EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostingId);
        }

        stkMeshBulkData.modification_begin();
        if ( myProc == owningProc )
        {
           stkMeshBulkData.destroy_entity(elementToKill);
        }
        stkMeshBulkData.modification_end();

        stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdToKill);
        EXPECT_FALSE(stkMeshBulkData.is_valid(element));
    }
}

}
