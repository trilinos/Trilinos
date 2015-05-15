#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <stk_topology/topology.hpp>

namespace
{

int check_connectivity(const std::vector<std::vector<int64_t> >& elem_graph, const std::vector<std::vector<int64_t> > &via_side, int64_t element_id1, int64_t element_id2)
{
    int side=-1;

    if(element_id1 >=1 && element_id1 <=3 && element_id2 >=1 && element_id2 <=3)
    {
        const std::vector<int64_t>& conn_elements = elem_graph[element_id1-1];

        std::vector<int64_t>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), element_id2);
        if ( iter != conn_elements.end() )
        {
            int64_t index = iter - conn_elements.begin();
            side = via_side[element_id1-1][index];
        }
    }

    return side;
}

TEST(ElementGraph, check_graph_connectivity)
{
    // element1 --> element2 --> element3
    std::vector<std::vector<int64_t> > elem_graph = {
            {2},
            {1,3},
            {2}
    };

    std::vector<std::vector<int64_t> > via_side = {
            {5},
            {2,6},
            {4}
    };

    EXPECT_EQ(5, check_connectivity(elem_graph, via_side, 1, 2));
    EXPECT_EQ(2, check_connectivity(elem_graph, via_side, 2, 1));
    EXPECT_EQ(6, check_connectivity(elem_graph, via_side, 2, 3));
    EXPECT_EQ(4, check_connectivity(elem_graph, via_side, 3, 2));

    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 1, 3));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 3, 1));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 4, 1));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 1, 4));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 1, 1));
}

std::vector<std::pair<int64_t, int64_t> > skin_mesh(const std::vector<std::vector<int64_t> >& elem_graph, const std::vector<std::vector<int64_t> > &via_side,
   const std::vector<stk::topology> &element_topologies)
{
    std::vector<std::pair<int64_t, int64_t> > element_side_pairs;

    std::vector<int64_t> elem_sides;

    size_t num_elems = via_side.size();
    for(size_t i=0; i<num_elems; ++i)
    {
        const std::vector<int64_t>& internal_sides = via_side[i];
        size_t num_sides = element_topologies[i].num_sides();

        elem_sides.assign(num_sides, -1);
        for(size_t j=0; j<internal_sides.size(); ++j)
        {
            int64_t zeroBasedSideId = internal_sides[j] - 1;
            elem_sides[zeroBasedSideId] = internal_sides[j];
        }

        int64_t elementId = i+1;
        for(size_t j=0; j<num_sides; ++j)
        {
            if (elem_sides[j] == -1)
            {
                int64_t sideId = j+1;
                element_side_pairs.push_back(std::make_pair(elementId, sideId));
            }
        }
    }
    return element_side_pairs;
}

TEST(ElementGraph, skin_mesh_using_graph)
{
    // element1 --> element2 --> element3
    std::vector<std::vector<int64_t> > elem_graph = {
            {2},
            {1,3},
            {2}
    };

    std::vector<std::vector<int64_t> > via_side = {
            {5},
            {2,6},
            {4}
    };

    std::vector<stk::topology> element_topologies{
        stk::topology::HEXAHEDRON_8,
        stk::topology::HEXAHEDRON_8,
        stk::topology::HEXAHEDRON_8
    };

    std::vector<std::pair<int64_t, int64_t> > element_side_pairs = skin_mesh(elem_graph, via_side, element_topologies);

    std::vector<std::pair<int64_t,int64_t> >gold_element_side_pairs{
        {1,1},
        {1,2},
        {1,3},
        {1,4},
        {1,6},
        {2,1},
        {2,3},
        {2,4},
        {2,5},
        {3,1},
        {3,2},
        {3,3},
        {3,5},
        {3,6}
    };

    ASSERT_EQ(gold_element_side_pairs.size(), element_side_pairs.size());

    for (size_t i=0;i<gold_element_side_pairs.size();++i)
    {
        std::vector<std::pair<int64_t, int64_t> >::iterator iter = std::find(element_side_pairs.begin(), element_side_pairs.end(), gold_element_side_pairs[i]);
        EXPECT_TRUE(iter != element_side_pairs.end()) << "gold elem-side-pair=" << gold_element_side_pairs[i].first << ", " << gold_element_side_pairs[i].second;
    }
}

}
