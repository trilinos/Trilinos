#include "MetaDataTester.hpp"
#include <gtest/gtest.h>

namespace stk { namespace mesh { namespace unit_test {

MetaDataTester::MetaDataTester(size_t spatial_dimension,
                               RootTopologiesInduceOption root_topos_induce,
                               const std::vector<std::string>& rank_names)
    : MetaData(),
      m_rootTopologiesInduce(root_topos_induce)
{
    initialize(spatial_dimension);
}


Part & MetaDataTester::declare_internal_part( const std::string & p_name, EntityRank rank)
{
    std::string internal_name = impl::convert_to_internal_name(p_name);

    const bool force_dont_induce = !m_rootTopologiesInduce;
    return declare_part(internal_name, rank, force_dont_induce);
}

}}}

namespace {

TEST ( UnitTestMetaDataTester, rootTopologiesInduceOptionWorks)
{
    const size_t spatialDim = 3;

    {
        stk::mesh::unit_test::MetaDataTester
            meta_no_induce(spatialDim, stk::mesh::unit_test::MetaDataTester::RootTopologiesInduceOption::NO_INDUCE);

        stk::mesh::PartVector parts = meta_no_induce.get_parts();
        for (size_t i = 0; i < parts.size(); ++ i)
        {
            stk::mesh::Part &part_i = *parts[i];

            if (part_i.topology() != stk::topology::INVALID_TOPOLOGY)
            {
                EXPECT_TRUE(part_i.force_no_induce());
            }
        }
    }

    {
        stk::mesh::unit_test::MetaDataTester
            meta_induce(spatialDim, stk::mesh::unit_test::MetaDataTester::RootTopologiesInduceOption::INDUCE);

        stk::mesh::PartVector parts = meta_induce.get_parts();
        for (size_t i = 0; i < parts.size(); ++ i)
        {
            stk::mesh::Part &part_i = *parts[i];

            if (part_i.topology() != stk::topology::INVALID_TOPOLOGY)
            {
                EXPECT_FALSE(part_i.force_no_induce());
            }
        }
    }
}

}
