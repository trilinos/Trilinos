#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

namespace {
TEST(stk_topology_how_to, runtime_vs_compiletime_topology )
{
    stk::topology runtime_hex8 = stk::topology::HEX_8;

    typedef stk::topology::topology_type<stk::topology::HEX_8> compiletime_hex8;

    const unsigned compiletime_num_nodes = compiletime_hex8::num_nodes;

    EXPECT_EQ( runtime_hex8.num_nodes(), compiletime_num_nodes );

    //declare a static array with length given by compile-time num-nodes
    double compile_time_sized_array[compiletime_num_nodes];
    EXPECT_EQ(sizeof(compile_time_sized_array), sizeof(double)*compiletime_num_nodes);
}
}

