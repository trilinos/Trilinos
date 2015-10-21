#include <gtest/gtest.h>
#include <vector>
#include <map>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include "ElementGraphTester.hpp"

namespace {

class ElemGraphAddElementsToEmptyGraphTester : public stk::unit_test_util::MeshTestFixture
{
protected:
    ElemGraphAddElementsToEmptyGraphTester() : elementGraph(nullptr) { }
    ~ElemGraphAddElementsToEmptyGraphTester()
    {
        delete elementGraph;
    }

    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh("generated:1x1x4", auraOption);
    }

    void create_elem_graph()
    {
        elementGraph = new ElemElemGraphTester(get_bulk());
        EXPECT_EQ(0u, elementGraph->size());
    }

    void expect_graph_correct_after_moving_everything_to_proc0()
    {
        create_elem_graph();
        add_elements_to_graph();
        test_edges();
    }

    void add_elements_to_graph()
    {
        stk::mesh::EntityVector elements_to_add;
        get_bulk().get_entities(stk::topology::ELEMENT_RANK, get_bulk().mesh_meta_data().locally_owned_part(), elements_to_add);
        elementGraph->add_elements(elements_to_add);
        EXPECT_EQ(elements_to_add.size(), elementGraph->size());
    }

    void test_num_edges()
    {
        std::vector<size_t> numGraphEdgesPerProc = {1, 2, 2, 1};
        size_t numGraphEdges = numGraphEdgesPerProc[get_bulk().parallel_rank()];
        EXPECT_EQ(numGraphEdges, elementGraph->num_edges());
        EXPECT_EQ(numGraphEdges, elementGraph->num_parallel_edges());
    }

    void test_edges()
    {
        test_num_edges();

        GraphEdges graphEdgesThisProc = goldElement1ToElement2SideOrdinalsPerOwningProc[get_bulk().parallel_rank()];
        for(GraphEdgeMock const &graphEdge : graphEdgesThisProc)
            EXPECT_EQ(graphEdge.sideOrdinalConnectingElement1ToElement2, elementGraph->get_side_from_element1_to_element2(graphEdge.element1, graphEdge.element2));
    }
protected:
    ElemElemGraphTester *elementGraph;

    struct GraphEdgeMock
    {
        int element1;
        int element2;
        int sideOrdinalConnectingElement1ToElement2;
    };

    typedef std::vector<GraphEdgeMock> GraphEdges;
    std::vector<GraphEdges> goldElement1ToElement2SideOrdinalsPerOwningProc =
    {
        {{1, 2, 5}, {1, 3, -1}, {1, 4, -1}},
        {{2, 1, 4}, {2, 3,  5}, {2, 4, -1}},
        {{3, 2, 4}, {3, 4,  5}, {3, 1, -1}},
        {{4, 3, 4}, {4, 1, -1}, {4, 2, -1}}
    };
};
TEST_F(ElemGraphAddElementsToEmptyGraphTester, withAura)
{
    run_test_on_num_procs(4, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphAddElementsToEmptyGraphTester, withoutAura)
{
    run_test_on_num_procs(4, stk::mesh::BulkData::NO_AUTO_AURA);
}

}
