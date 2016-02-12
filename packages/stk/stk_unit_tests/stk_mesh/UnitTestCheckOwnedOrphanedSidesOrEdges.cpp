#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "../../stk_mesh/stk_mesh/base/GetEntities.hpp"
#include "../../stk_util/stk_util/parallel/ParallelReduce.hpp"

namespace
{

class MeshCheckerOwnedOrphans : public stk::unit_test_util::MeshFixture
{
protected:
    MeshCheckerOwnedOrphans()
    {
        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    }
};

std::vector<stk::mesh::EntityKey> get_orphaned_owned_sides(const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()), sides);
    std::vector<stk::mesh::EntityKey> badSides;
    for(stk::mesh::Entity side : sides)
    {
        size_t num_elements = bulkData.num_elements(side);
        if(num_elements == 0)
            badSides.push_back(bulkData.entity_key(side));
    }
    return badSides;
}

void fill_mesh_with_orphaned_owned_sides(stk::mesh::BulkData& bulkData)
{
    bulkData.modification_begin();
    if(bulkData.parallel_rank()==1)
    {
        stk::mesh::EntityVector nodes(4);
        for(stk::mesh::EntityId id=1;id<=4;++id)
            nodes[id-1] = bulkData.declare_entity(stk::topology::NODE_RANK, id);

        stk::mesh::Entity face = bulkData.declare_entity(stk::topology::FACE_RANK, 1, bulkData.mesh_meta_data().get_topology_root_part(stk::topology::QUADRILATERAL_4));
        for(size_t i=0;i<nodes.size();++i)
            bulkData.declare_relation(face, nodes[i], i);
    }
    bulkData.modification_end();
}

TEST_F(MeshCheckerOwnedOrphans, check_mesh_without_orphaned_owned_sides)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        std::vector<stk::mesh::EntityKey> orphanedOwnedSides = get_orphaned_owned_sides(get_bulk());
        EXPECT_TRUE(orphanedOwnedSides.empty());
    }
}

TEST_F(MeshCheckerOwnedOrphans, check_mesh_with_orphaned_owned_sides)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        std::vector<std::vector<stk::mesh::EntityKey> > badSidesPerProc = {
                {},
                {stk::mesh::EntityKey(stk::topology::FACE_RANK, 1)}
        };

        fill_mesh_with_orphaned_owned_sides(get_bulk());
        std::vector<stk::mesh::EntityKey> orphanedOwnedSides = get_orphaned_owned_sides(get_bulk());
        EXPECT_EQ(badSidesPerProc[get_bulk().parallel_rank()].size(), orphanedOwnedSides.size());
        EXPECT_EQ(badSidesPerProc[get_bulk().parallel_rank()], orphanedOwnedSides);
    }
}

}
