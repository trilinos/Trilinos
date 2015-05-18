#include <gtest/gtest.h>

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace
{

int check_connectivity(const std::vector<std::vector<int64_t> >& elem_graph, const std::vector<std::vector<int64_t> > &via_side, int64_t element_id1, int64_t element_id2)
{
    int side=-1;

    if(element_id1 >=0 && element_id1 <=2 && element_id2 >=0 && element_id2 <=2)
    {
        const std::vector<int64_t>& conn_elements = elem_graph[element_id1];

        std::vector<int64_t>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), element_id2);
        if ( iter != conn_elements.end() )
        {
            int64_t index = iter - conn_elements.begin();
            side = via_side[element_id1][index];
        }
    }

    return side;
}

TEST(ElementGraph, check_graph_connectivity)
{
    // element0 --> element1 --> element2
    std::vector<std::vector<int64_t> > elem_graph = {
            {1},
            {0,2},
            {1}
    };

    std::vector<std::vector<int64_t> > via_side = {
            {4},
            {1,5},
            {3}
    };

    EXPECT_EQ(4, check_connectivity(elem_graph, via_side, 0, 1));
    EXPECT_EQ(1, check_connectivity(elem_graph, via_side, 1, 0));
    EXPECT_EQ(5, check_connectivity(elem_graph, via_side, 1, 2));
    EXPECT_EQ(3, check_connectivity(elem_graph, via_side, 2, 1));

    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 0, 2));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 2, 0));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 3, 0));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 0, 3));
    EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 0, 0));
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

        if (internal_sides.size() < num_sides)
        {
            elem_sides.assign(num_sides, -1);
            for(size_t j=0; j<internal_sides.size(); ++j)
            {
                int64_t zeroBasedSideId = internal_sides[j];
                elem_sides[zeroBasedSideId] = internal_sides[j];
            }

            int64_t elementId = i;
            for(size_t j=0; j<num_sides; ++j)
            {
                if (elem_sides[j] == -1)
                {
                    int64_t sideId = j;
                    element_side_pairs.push_back(std::make_pair(elementId, sideId));
                }
            }
        }
    }
    return element_side_pairs;
}

TEST(ElementGraph, skin_mesh_using_graph)
{
    // element0 --> element1 --> element2
    std::vector<std::vector<int64_t> > elem_graph = {
            {1},
            {0,2},
            {1}
    };

    std::vector<std::vector<int64_t> > via_side = {
            {4},
            {1,5},
            {3}
    };

    std::vector<stk::topology> element_topologies{
        stk::topology::HEXAHEDRON_8,
        stk::topology::HEXAHEDRON_8,
        stk::topology::HEXAHEDRON_8
    };

    std::vector<std::pair<int64_t, int64_t> > element_side_pairs = skin_mesh(elem_graph, via_side, element_topologies);

    std::vector<std::pair<int64_t,int64_t> >gold_element_side_pairs{
        {0,0},
        {0,1},
        {0,2},
        {0,3},
        {0,5},
        {1,0},
        {1,2},
        {1,3},
        {1,4},
        {2,0},
        {2,1},
        {2,2},
        {2,4},
        {2,5}
    };

    ASSERT_EQ(gold_element_side_pairs.size(), element_side_pairs.size());

    for (size_t i=0;i<gold_element_side_pairs.size();++i)
    {
        std::vector<std::pair<int64_t, int64_t> >::iterator iter = std::find(element_side_pairs.begin(), element_side_pairs.end(), gold_element_side_pairs[i]);
        EXPECT_TRUE(iter != element_side_pairs.end()) << "gold elem-side-pair=" << gold_element_side_pairs[i].first << ", " << gold_element_side_pairs[i].second;
    }
}

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    size_t local_id = 0;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            local_id_to_element_entity[local_id] = bucket[j];
            element_topologies[local_id] = bucket.topology();
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
        std::vector<std::pair<int64_t,int64_t> > elem_side_pairs;
        stk::mesh::EntityVector side_nodes;
        stk::mesh::EntityVector connected_elements;
        for(size_t j=0; j<bucket.size(); ++j)
        {
            size_t local_elem_id = bulkData.local_id(bucket[j]);
            const stk::mesh::Entity* elem_nodes = bucket.begin_nodes(j);
            elem_side_pairs.clear();
            for(unsigned side_index=0; side_index<num_sides; ++side_index)
            {
                unsigned num_side_nodes = topology.side_topology(side_index).num_nodes();
                side_nodes.resize(num_side_nodes);
                topology.side_nodes(elem_nodes, side_index, side_nodes.begin());
                connected_elements.clear();
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
    std::vector<double> wall_times;
    wall_times.reserve(10);
    std::vector<std::string> msgs;
    msgs.reserve(10);

    std::vector<size_t> mem_usage;

    wall_times.push_back(stk::wall_time());
    msgs.push_back("program-start");
    mem_usage.push_back(stk::get_memory_usage_now());

    if(stk::parallel_machine_size(comm) == 1)
    {
        std::string dimension = unitTestUtils::getOption("--zdim", "3");
        const int zdim = std::atoi(dimension.c_str());

        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        std::ostringstream os;
        os << "generated:1x1x" << zdim;
        bool check_results = true;

        std::string filename = os.str();

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(zdim, numElems);

        stk::mesh::EntityVector local_id_to_element_entity(numElems, 0);
        std::vector<stk::topology> element_topologies(numElems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);
        size_t expectedNumElems = counts[stk::topology::ELEM_RANK];

        if ( check_results )
        {
            ASSERT_EQ(expectedNumElems, local_id_to_element_entity.size());

            std::vector<stk::mesh::EntityId> element_global_ids(numElems);
            for(int i=0; i<numElems; ++i)
            {
                stk::mesh::EntityId expectedId = static_cast<stk::mesh::EntityId>(i+1);
                EXPECT_EQ(expectedId, bulkData.identifier(local_id_to_element_entity[i])) << "elem id for i=" << i << ": "<< bulkData.identifier(local_id_to_element_entity[i]);
            }
        }

        std::vector<std::vector<int64_t> > elem_graph(numElems);
        std::vector<std::vector<int64_t> > via_sides(numElems);

        fill_graph(bulkData, elem_graph, via_sides);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        if ( check_results )
        {
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

    for(size_t i=0;i<wall_times.size();++i)
    {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
    }

    for(size_t i=0;i<mem_usage.size();++i)
     {
         std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
     }
}

TEST(ElementGraph, test_create_skin_for_stk_mesh)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<double> wall_times;
    wall_times.reserve(10);
    std::vector<std::string> msgs;
    msgs.reserve(10);

    std::vector<size_t> mem_usage;

    wall_times.push_back(stk::wall_time());
    msgs.push_back("program-start");
    mem_usage.push_back(stk::get_memory_usage_now());

    if(stk::parallel_machine_size(comm) == 1)
    {
        std::string dimension = unitTestUtils::getOption("--zdim", "3");
        const int zdim = std::atoi(dimension.c_str());

        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& skin_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(skin_part);
        stk::mesh::BulkData bulkData(meta, comm);
        std::ostringstream os;
        os << "generated:" << zdim << "x" << zdim << "x" << zdim;
        bool check_results = false;

        std::string filename = os.str();

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numElems = counts[stk::topology::ELEM_RANK];
        if ( check_results )
        {
            EXPECT_EQ(zdim, numElems);
        }

        stk::mesh::EntityVector local_id_to_element_entity(numElems, 0);
        std::vector<stk::topology> element_topologies(numElems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);
        size_t expectedNumElems = counts[stk::topology::ELEM_RANK];

        if ( check_results )
        {
            ASSERT_EQ(expectedNumElems, local_id_to_element_entity.size());

            std::vector<stk::mesh::EntityId> element_global_ids(numElems);
            for(int i=0; i<numElems; ++i)
            {
                stk::mesh::EntityId expectedId = static_cast<stk::mesh::EntityId>(i+1);
                EXPECT_EQ(expectedId, bulkData.identifier(local_id_to_element_entity[i])) << "elem id for i=" << i << ": "<< bulkData.identifier(local_id_to_element_entity[i]);
            }
        }

        std::vector<std::vector<int64_t> > elem_graph(numElems);
        std::vector<std::vector<int64_t> > via_sides(numElems);

        fill_graph(bulkData, elem_graph, via_sides);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<std::pair<int64_t, int64_t> > elem_side_pairs = skin_mesh(elem_graph, via_sides, element_topologies);

        bulkData.modification_begin();

        for(size_t face_index=0; face_index<elem_side_pairs.size(); ++face_index)
        {
            stk::mesh::Entity element = local_id_to_element_entity[elem_side_pairs[face_index].first];
            stk::mesh::EntityId face_global_id = face_index + 1;
            stk::mesh::declare_element_side(bulkData, face_global_id, element, elem_side_pairs[face_index].second, &skin_part);
        }

        bulkData.modification_end_for_entity_creation(stk::topology::FACE_RANK);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after create-faces");
        mem_usage.push_back(stk::get_memory_usage_now());

        stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());

        if ( check_results )
        {
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

    for(size_t i=0;i<wall_times.size();++i)
    {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
    }

    for(size_t i=0;i<mem_usage.size();++i)
     {
         std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
     }
}

TEST(ElementGraph, test_skin_mesh)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<double> wall_times;
    wall_times.reserve(10);
    std::vector<std::string> msgs;
    msgs.reserve(10);

    std::vector<size_t> mem_usage;

    wall_times.push_back(stk::wall_time());
    msgs.push_back("program-start");
    mem_usage.push_back(stk::get_memory_usage_now());

    if(stk::parallel_machine_size(comm) == 1)
    {
        std::string dimension = unitTestUtils::getOption("--zdim", "3");
        const int zdim = std::atoi(dimension.c_str());

        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& skin_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(skin_part);
        stk::mesh::BulkData bulkData(meta, comm);
        std::ostringstream os;

        os << "generated:" << zdim << "x" << zdim << "x" << zdim;

        std::string filename = os.str();
        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        stk::mesh::PartVector parts(1, &skin_part);
        stk::mesh::skin_mesh(bulkData, parts);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after stk::skin_mesh");
        mem_usage.push_back(stk::get_memory_usage_now());
    }

    for(size_t i = 0; i < wall_times.size(); ++i)
    {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
    }

    for(size_t i = 0; i < mem_usage.size(); ++i)
    {
        std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
    }
}

}
