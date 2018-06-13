#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace
{

class StkToolsC : public stk::unit_test_util::MeshFixture
{};

TEST_F(StkToolsC, DeleteMeshExceptSpecifiedElems)
{
    const std::string unNamed = "mesh not specified";
    const std::string meshName = stk::unit_test_util::get_option("-i", unNamed);
    ThrowRequireMsg(meshName!=unNamed, "Please specify mesh with -i option.");
    setup_mesh(meshName, stk::mesh::BulkData::NO_AUTO_AURA);

    int invalidElemId = -1;
    int inputElemId = stk::unit_test_util::get_command_line_option("-e", invalidElemId);
    ThrowRequireMsg(inputElemId!=invalidElemId, "Please specify element with -e.");

    const stk::mesh::BucketVector &buckets = get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().locally_owned_part());

    std::set<stk::mesh::EntityId> elemIdsToKeep;
    elemIdsToKeep.insert(inputElemId);
//     elemIdsToKeep.insert(21081);
//     elemIdsToKeep.insert(21080);
//     elemIdsToKeep.insert(20832);
//     elemIdsToKeep.insert(20833);
//     elemIdsToKeep.insert(115099);
//     elemIdsToKeep.insert(105459);
//     elemIdsToKeep.insert(106976);
//     elemIdsToKeep.insert(108648);

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(get_bulk(), entityCounts);

    std::cout << "num entities = " << entityCounts[stk::topology::ELEMENT_RANK] << std::endl;

    stk::mesh::EntityVector elementsToDestroy;
    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket& bucket = *buckets[i];
        for(size_t j=0;j<bucket.size();j++)
        {
            stk::mesh::Entity element = bucket[j];
            stk::mesh::EntityId elemId = get_bulk().identifier(element);
            if (elemIdsToKeep.find(elemId) == elemIdsToKeep.end())
            {
                elementsToDestroy.push_back(element);
            }
            else
            {
                std::cout << "keeping element ID: " << elemId << std::endl;
            }
        }
    }
    destroy_elements(get_bulk(), elementsToDestroy);
    stk::io::write_mesh("modified.e", get_bulk());
}

TEST_F(StkToolsC, FlipElementConnectivity)
{
    const std::string unNamed = "mesh not specified";
    const std::string meshName = stk::unit_test_util::get_option("-i", unNamed);
    ThrowRequireMsg(meshName!=unNamed, "Please specify mesh with -i option.");
    setup_mesh(meshName, stk::mesh::BulkData::NO_AUTO_AURA);

    int invalidBlockId = -1;
    int inputBlockId = stk::unit_test_util::get_command_line_option("-b", invalidBlockId);
    ThrowRequireMsg(inputBlockId!=invalidBlockId, "Please specify block with -b.");

    std::ostringstream os;
    os << "block_" << inputBlockId;
    stk::mesh::Part* elemBlock = get_meta().get_part(os.str());

    stk::mesh::EntityVector elems;
    stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
    stk::mesh::get_selected_entities(stk::mesh::Selector(*elemBlock), buckets, elems);

    get_bulk().modification_begin();
    for (auto elem : elems)
    {
        stk::topology topology = get_bulk().bucket(elem).topology();
        ThrowRequireMsg(topology == stk::topology::HEX_8, "Input block must have HEX_8 topology but found topology " << topology);

        stk::mesh::EntityVector storedNodes;
        const stk::mesh::Entity* node = get_bulk().begin_nodes(elem);
        unsigned numNodes = get_bulk().num_nodes(elem);
        for (unsigned i = 0; i < numNodes; ++i)
        {
            storedNodes.push_back(node[i]);
        }
        for (unsigned i = 0; i < numNodes; ++i)
        {
            get_bulk().destroy_relation(elem, storedNodes[i], i);
        }
        for (unsigned i = 0; i < numNodes; ++i)
        {
            unsigned correctNodeId = (i+4)%8;
            get_bulk().declare_relation(elem, storedNodes[correctNodeId], i);
        }
    }
    get_bulk().modification_end();
    stk::io::write_mesh("flipped.e", get_bulk());
}
}
