#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <vector>

namespace {

//EquivalentElements
TEST(stk_topology_understanding, equivalent_elements)
{
    std::pair<bool, unsigned> areElementsEquivalent;

    {
        if (stk::topology::topology_type<stk::topology::HEX_8>::num_permutations > 1) {
            unsigned hex1[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
            unsigned hex2[8] = { 4, 7, 6, 5, 0, 3, 2, 1 };
            unsigned hex3[8] = { 4, 5, 6, 7, 0, 1, 2, 3 };

            stk::topology hex8 = stk::topology::HEX_8;

            areElementsEquivalent = hex8.equivalent(hex1, hex2);
            EXPECT_TRUE(areElementsEquivalent.first);
            areElementsEquivalent = hex8.equivalent(hex1, hex3);
            EXPECT_FALSE(areElementsEquivalent.first);
        }
    }

    {
        unsigned triangle_1[3] = {0, 1, 2};
        unsigned triangle_2[3] = {0, 2, 1};

        stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;

        areElementsEquivalent = triangular_shell.equivalent(triangle_1, triangle_2);

        EXPECT_TRUE(areElementsEquivalent.first);

        unsigned permutation_index = areElementsEquivalent.second;
        unsigned goldValue = 3;
        EXPECT_EQ(goldValue, permutation_index); // From previous unit test, this is the 4th permutation
    }

}
//EndEquivalentElements
}

