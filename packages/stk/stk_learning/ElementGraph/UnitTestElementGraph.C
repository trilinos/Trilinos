#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

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

void set_local_ids_and_fill_element_entities(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    size_t local_id = 0;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            local_id_to_element_entity[local_id] = bucket[j];
            bulkData.set_local_id(bucket[j], local_id);
            local_id++;
        }
    }
}

void fill_graph(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        stk::topology topology = bucket.topology();
        unsigned num_sides = topology.num_sides();
        for(size_t j=0; j<bucket.size(); ++j)
        {
            size_t local_elem_id = bulkData.local_id(bucket[j]);
            const stk::mesh::Entity* elem_nodes = bucket.begin_nodes(j);
            std::vector<std::pair<int64_t,int64_t> > elem_side_pairs;
            for(unsigned side_index=0; side_index<num_sides; ++side_index)
            {
                unsigned num_side_nodes = topology.side_topology(side_index).num_nodes();
                stk::mesh::EntityVector side_nodes(num_side_nodes);
                topology.side_nodes(elem_nodes, side_index, side_nodes.begin());
                stk::mesh::EntityVector connected_elements;
                stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, num_side_nodes, side_nodes.data(), connected_elements);
                for(size_t elem_index=0; elem_index<connected_elements.size(); ++elem_index)
                {
                    if (connected_elements[elem_index] != bucket[j])
                    {
                        elem_side_pairs.push_back(std::make_pair(bulkData.local_id(connected_elements[elem_index]),side_index));
                    }
                }
            }

            std::sort(elem_side_pairs.begin(), elem_side_pairs.end());
            std::vector<std::pair<int64_t,int64_t> >::iterator new_end = std::unique(elem_side_pairs.begin(), elem_side_pairs.end());
            elem_side_pairs.resize(new_end - elem_side_pairs.begin());
            for(size_t index=0; index<elem_side_pairs.size(); ++index)
            {
                elem_graph[local_elem_id].push_back(elem_side_pairs[index].first);
                via_sides[local_elem_id].push_back(elem_side_pairs[index].second);
            }
        }
    }
}

TEST(ElementGraph, test_graph_creation_using_stk_mesh)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) == 1)
    {
        std::string dimension = unitTestUtils::getOption("--zdim", "3");
        const int zdim = std::atoi(dimension.c_str());

        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        std::ostringstream os;
        os << "generated:1x1x" << zdim;
        std::string filename = os.str();

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(zdim, numElems);

        stk::mesh::EntityVector local_id_to_element_entity(numElems, 0);
        set_local_ids_and_fill_element_entities(bulkData, local_id_to_element_entity);
        size_t expectedNumElems = counts[stk::topology::ELEM_RANK];
        ASSERT_EQ(expectedNumElems, local_id_to_element_entity.size());
        std::vector<stk::mesh::EntityId> element_global_ids(numElems);
        for(int i=0; i<numElems; ++i)
        {
            stk::mesh::EntityId expectedId = static_cast<stk::mesh::EntityId>(i+1);
            EXPECT_EQ(expectedId, bulkData.identifier(local_id_to_element_entity[i])) << "elem id for i=" << i << ": "<< bulkData.identifier(local_id_to_element_entity[i]);
        }

        std::vector<std::vector<int64_t> > elem_graph(numElems);
        std::vector<std::vector<int64_t> > via_sides(numElems);

        fill_graph(bulkData, elem_graph, via_sides);

        int64_t left_side_id = 4;
        int64_t right_side_id = 5;

        for(size_t i=0; i<elem_graph.size(); ++i)
        {
            const std::vector<int64_t>& conn_elements = elem_graph[i];
            if (i == 0)
            {
                ASSERT_EQ(1u, conn_elements.size());
                EXPECT_EQ(1, conn_elements[0]);
                EXPECT_EQ(right_side_id, via_sides[i][0]);
            }
            else if (i == elem_graph.size() - 1)
            {
                int64_t second_to_last_element_index = elem_graph.size() - 2;
                ASSERT_EQ(1u, conn_elements.size());
                EXPECT_EQ(second_to_last_element_index, conn_elements[0]);
                EXPECT_EQ(left_side_id, via_sides[i][0]);
            }
            else
            {
                ASSERT_EQ(2u, conn_elements.size());
                int64_t element_to_the_left = i-1;
                int64_t element_to_the_right = i+1;
                EXPECT_EQ(element_to_the_left, conn_elements[0]);
                EXPECT_EQ(element_to_the_right, conn_elements[1]);
                EXPECT_EQ(left_side_id, via_sides[i][0]);
                EXPECT_EQ(right_side_id, via_sides[i][1]);
            }
        }
    }
}

}
