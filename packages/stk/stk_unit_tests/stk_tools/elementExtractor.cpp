#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/DestroyElements.hpp>

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
    // elemIdsToKeep.insert(3403);
    // elemIdsToKeep.insert(3405);
    // elemIdsToKeep.insert(3402);
    // elemIdsToKeep.insert(2781);
    // elemIdsToKeep.insert(2783);
    // elemIdsToKeep.insert(2780);
    // elemIdsToKeep.insert(6417);

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
    stk::io::write_mesh("modifiedGapRemovalExample.e", get_bulk());
}
}
