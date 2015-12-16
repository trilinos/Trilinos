#include "gtest/gtest.h"
#include <mpi.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

class BulkDataIdMapperTest : public stk::unit_test_util::MeshFixture
{
protected:
    BulkDataIdMapperTest()
    {
        if(stk::parallel_machine_size(get_comm()) == 1)
        {
            setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);
            stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEMENT_RANK), elements);
            build_elem_to_local_id_vector();
        }
    }
    void build_elem_to_local_id_vector()
    {
        elemToLocalIds.resize(get_num_entities());
        for(size_t i=0; i<elements.size(); i++)
            elemToLocalIds[elements[i].local_offset()] = i;
    }
    unsigned get_num_entities()
    {
        const unsigned numInvalidEntities = 1;
        unsigned numEntities = numInvalidEntities;
        std::vector<unsigned> countPerRank;
        stk::mesh::count_entities(get_meta().universal_part(), get_bulk(), countPerRank);
        for(unsigned count : countPerRank)
            numEntities += count;
        return numEntities;
    }
    stk::mesh::EntityVector elements;
    std::vector<stk::mesh::impl::LocalId> elemToLocalIds;
};

TEST_F(BulkDataIdMapperTest, localToGlobal)
{
    if(stk::parallel_machine_size(get_comm()) == 1)
    {
        stk::mesh::impl::BulkDataIdMapper idMapper(get_bulk(), elements, elemToLocalIds);
        EXPECT_EQ(1u, idMapper.localToGlobal(0));
    }
}

TEST_F(BulkDataIdMapperTest, globalToLocal)
{
    if(stk::parallel_machine_size(get_comm()) == 1)
    {
        stk::mesh::impl::BulkDataIdMapper idMapper(get_bulk(), elements, elemToLocalIds);
        EXPECT_EQ(0, idMapper.globalToLocal(1));
    }
}

}
