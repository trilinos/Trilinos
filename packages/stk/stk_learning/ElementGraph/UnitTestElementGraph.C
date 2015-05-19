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
#include <stk_util/parallel/CommSparse.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace
{

class BulkDataElementGraphTester : public stk::mesh::BulkData
{

public:

    BulkDataElementGraphTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }

    ~BulkDataElementGraphTester(){}

    bool my_internal_modification_end_for_skin_mesh(stk::mesh::EntityRank entity_rank, modification_optimization opt, stk::mesh::Selector selectedToSkin,
            const stk::mesh::Selector * only_consider_second_element_from_this_selector = 0)
    {
        return this->internal_modification_end_for_skin_mesh(entity_rank, opt, selectedToSkin, only_consider_second_element_from_this_selector);
    }
};

int check_connectivity(const std::vector<std::vector<int64_t> >& elem_graph, const std::vector<std::vector<int64_t> > &via_side, int64_t element_id1, int64_t element_id2);

std::vector<std::pair<int64_t, int64_t> >
skin_mesh(const std::vector<std::vector<int64_t> > &via_side, const std::vector<stk::topology> &element_topologies);

void fill_parallel_graph_with_test(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides);
void fill_parallel_graph(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides);

void add_possibly_connected_elements_to_graph_using_side_nodes(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides,
                                                                               const stk::mesh::EntityVector& side_nodes);

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements_to_communicate);

stk::mesh::EntityVector get_elements_to_communicate(stk::mesh::BulkData& bulkData);

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies);

void add_element_side_pair_for_unused_sides(int64_t elementId, stk::topology topology, const std::vector<int64_t> &internal_sides, std::vector<std::pair<int64_t, int64_t> >& element_side_pairs);

void fill_graph(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides);

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

    std::vector<std::pair<int64_t, int64_t> > element_side_pairs = skin_mesh(via_side, element_topologies);

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

TEST(ElementGraph, create_element_graph_serial)
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
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);
        std::ostringstream os;

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(4, numElems);

        stk::mesh::EntityVector local_id_to_element_entity(numElems, 0);
        std::vector<stk::topology> element_topologies(numElems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);
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

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

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

        if (stk::parallel_machine_rank(comm) == 0)
        {
            for(size_t i=0;i<wall_times.size();++i)
            {
                std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
            }

            for(size_t i=0;i<mem_usage.size();++i)
            {
                std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
            }
        }
    }
}

TEST(ElementGraph, create_element_graph_parallel)
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

    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(2, numLocallyOwnedElems);

        stk::mesh::EntityVector local_id_to_element_entity(numLocallyOwnedElems, 0);
        std::vector<stk::topology> element_topologies(numLocallyOwnedElems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        std::vector<std::vector<int64_t> > elem_graph(numLocallyOwnedElems);
        std::vector<std::vector<int64_t> > via_sides(numLocallyOwnedElems);

        fill_parallel_graph_with_test(bulkData, elem_graph, via_sides);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        unsigned element_to_test_local_id = std::numeric_limits<unsigned>::max();
        int64_t side_id = -1;
        int64_t left_side_id = 4;
        int64_t right_side_id = 5;

        if (stk::parallel_machine_rank(comm) == 0)
        {
            element_to_test_local_id = 1;
            side_id = right_side_id;
        }
        else
        {
            element_to_test_local_id = 0;
            side_id = left_side_id;
        }

        for(size_t i=0;i<elem_graph.size();++i)
        {
            if (i != element_to_test_local_id)
            {
                EXPECT_TRUE(elem_graph[i].empty());
                EXPECT_TRUE(via_sides[i].empty());
            }
            else
            {
                ASSERT_TRUE(elem_graph[i].size()==1);
                ASSERT_TRUE(via_sides[i].size()==1);
                EXPECT_EQ(-1, elem_graph[i][0]);
                EXPECT_EQ(side_id, via_sides[i][0]);
            }
        }

        if (stk::parallel_machine_rank(comm) == 0)
        {
            for(size_t i=0;i<wall_times.size();++i)
            {
                std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
            }

            for(size_t i=0;i<mem_usage.size();++i)
            {
                std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
            }
        }
    }
}

TEST(ElementGraph, skin_mesh_using_element_graph_serial)
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
        BulkDataElementGraphTester bulkData(meta, comm);
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

        std::vector<std::pair<int64_t, int64_t> > elem_side_pairs = skin_mesh(via_sides, element_topologies);

        bulkData.modification_begin();

        for(size_t face_index=0; face_index<elem_side_pairs.size(); ++face_index)
        {
            stk::mesh::Entity element = local_id_to_element_entity[elem_side_pairs[face_index].first];
            stk::mesh::EntityId face_global_id = face_index + 1;
            stk::mesh::declare_element_side(bulkData, face_global_id, element, elem_side_pairs[face_index].second, &skin_part);
        }

        stk::mesh::Selector element_selector = bulkData.mesh_meta_data().locally_owned_part();
        bulkData.my_internal_modification_end_for_skin_mesh(stk::topology::FACE_RANK, stk::mesh::BulkData::MOD_END_SORT, element_selector, NULL);

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

        if (stk::parallel_machine_rank(comm) == 0)
        {
            for(size_t i=0;i<wall_times.size();++i)
            {
                std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
            }

            for(size_t i=0;i<mem_usage.size();++i)
            {
                std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
            }
        }
    }
}

TEST(ElementGraph, skin_mesh_using_element_graph_parallel)
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

    if(stk::parallel_machine_size(comm) <= 2)
    {
        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& skin_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(skin_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        unsigned num_locally_owned_elems = stk::mesh::count_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK));

        stk::mesh::EntityVector local_id_to_element_entity(num_locally_owned_elems, 0);
        std::vector<stk::topology> element_topologies(num_locally_owned_elems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        std::vector<std::vector<int64_t> > elem_graph(num_locally_owned_elems);
        std::vector<std::vector<int64_t> > via_sides(num_locally_owned_elems);

        fill_graph(bulkData, elem_graph, via_sides);
        if (stk::parallel_machine_size(comm) > 1)
        {
            fill_parallel_graph_with_test(bulkData, elem_graph, via_sides);
        }

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<std::pair<int64_t, int64_t> > elem_side_pairs = skin_mesh(via_sides, element_topologies);

        bulkData.modification_begin();

        int offset_per_proc = 40*bulkData.parallel_rank();

        for(size_t face_index=0; face_index<elem_side_pairs.size(); ++face_index)
        {
            stk::mesh::Entity element = local_id_to_element_entity[elem_side_pairs[face_index].first];
            stk::mesh::EntityId face_global_id = face_index + 1 + offset_per_proc;
            stk::mesh::declare_element_side(bulkData, face_global_id, element, elem_side_pairs[face_index].second, &skin_part);
        }

        stk::mesh::Selector element_selector = bulkData.mesh_meta_data().locally_owned_part();
        bulkData.my_internal_modification_end_for_skin_mesh(stk::topology::FACE_RANK, stk::mesh::BulkData::MOD_END_SORT, element_selector, NULL);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after create-faces");
        mem_usage.push_back(stk::get_memory_usage_now());

        stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());

        std::vector<size_t> counts;
        stk::mesh::Selector skin = skin_part;
        stk::mesh::comm_mesh_counts(bulkData, counts, &skin);

        size_t num_faces = counts[stk::topology::FACE_RANK];
        EXPECT_EQ(18u, num_faces);

        if (stk::parallel_machine_rank(comm) == 0)
        {
            for(size_t i=0;i<wall_times.size();++i)
            {
                std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
            }

            for(size_t i=0;i<mem_usage.size();++i)
            {
                std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
            }
        }
    }
}

TEST(ElementGraph, create_faces_using_element_graph_serial)
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
        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(new_faces_part);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x3", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numElems = counts[stk::topology::ELEM_RANK];

        stk::mesh::EntityVector local_id_to_element_entity(numElems, 0);
        std::vector<stk::topology> element_topologies(numElems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        std::vector<std::vector<int64_t> > elem_graph(numElems);
        std::vector<std::vector<int64_t> > via_sides(numElems);

        fill_graph(bulkData, elem_graph, via_sides);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        bulkData.modification_begin();

        stk::mesh::EntityId face_global_id = 1;
        for(size_t i=0; i<elem_graph.size(); ++i)
        {
            const std::vector<int64_t>& connected_elements = elem_graph[i];
            stk::mesh::Entity element1 = local_id_to_element_entity[i];

            int64_t signed_i = i;

            for(size_t j=0; j<connected_elements.size(); ++j)
            {
                if (signed_i > connected_elements[j])
                {
                    stk::mesh::Entity face = stk::mesh::declare_element_side(bulkData, face_global_id, element1, via_sides[i][j], &new_faces_part);

                    const stk::mesh::Entity* side_nodes = bulkData.begin_nodes(face);
                    unsigned num_side_nodes = bulkData.num_nodes(face);
                    stk::mesh::EntityVector side_nodes_vec(side_nodes, side_nodes+num_side_nodes);

                    stk::mesh::Entity element2 = local_id_to_element_entity[connected_elements[j]];
                    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm = stk::mesh::get_ordinal_and_permutation(bulkData, element2, stk::topology::FACE_RANK, side_nodes_vec);
                    bulkData.declare_relation(element2, face, ord_and_perm.first, ord_and_perm.second);

                    face_global_id++;
                }
            }

            std::vector<std::pair<int64_t, int64_t> > element_side_pairs;
            add_element_side_pair_for_unused_sides(i, element_topologies[i], via_sides[i], element_side_pairs);

            for(size_t j=0;j<element_side_pairs.size();j++)
            {
                stk::mesh::declare_element_side(bulkData, face_global_id, element1, element_side_pairs[j].second, &new_faces_part);
                face_global_id++;
            }
        }

        bulkData.modification_end_for_entity_creation(stk::topology::FACE_RANK);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after create-faces");
        mem_usage.push_back(stk::get_memory_usage_now());

        unsigned num_faces = stk::mesh::count_selected_entities(new_faces_part, bulkData.buckets(stk::topology::FACE_RANK));

        EXPECT_EQ(16u, num_faces);
        EXPECT_EQ(num_faces, face_global_id-1);

        stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());

        if (stk::parallel_machine_rank(comm) == 0)
        {
            for(size_t i=0;i<wall_times.size();++i)
            {
                std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
            }

            for(size_t i=0;i<mem_usage.size();++i)
            {
                std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
            }
        }
    }
}

std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim)
{
    std::ostringstream os;
    os << "generated:" << xdim << "x" << ydim << "x" << zdim;
    return os.str();
}

TEST(ElementGraph, compare_performance_skin_mesh)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    //wall_times.push_back(stk::wall_time());
    // mem_usage.push_back(stk::get_memory_usage_now());

    std::string dimension = unitTestUtils::getOption("--zdim", "none");

    int xdim = 3;
    if ( dimension != "none")
    {
        xdim = std::atoi(dimension.c_str());
    }

    int ydim = xdim;
    int zdim = xdim * stk::parallel_machine_size(comm);

    std::string filename = get_name_of_generated_mesh(xdim, ydim, zdim);

    {
        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& skin_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(skin_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        {
            double wall_time_start = stk::wall_time();

            stk::mesh::PartVector parts(1, &skin_part);
            stk::mesh::skin_mesh(bulkData, parts);

            double elapsed_time = stk::wall_time() - wall_time_start;

            if (stk::parallel_machine_rank(comm) == 0)
            {
                std::cerr << "STK time: " << elapsed_time << std::endl;
            }
        }
    }

    {
        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& skin_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(skin_part);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        {
            double wall_time_start = stk::wall_time();

            unsigned num_locally_owned_elems = stk::mesh::count_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK));

            stk::mesh::EntityVector local_id_to_element_entity(num_locally_owned_elems, 0);
            std::vector<stk::topology> element_topologies(num_locally_owned_elems);
            set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

            std::vector<std::vector<int64_t> > elem_graph(num_locally_owned_elems);
            std::vector<std::vector<int64_t> > via_sides(num_locally_owned_elems);

            fill_graph(bulkData, elem_graph, via_sides);
            fill_parallel_graph(bulkData, elem_graph, via_sides);

            std::vector<std::pair<int64_t, int64_t> > elem_side_pairs = skin_mesh(via_sides, element_topologies);

            bulkData.modification_begin();

            int offset_per_proc = 100000*bulkData.parallel_rank();

            for(size_t face_index=0; face_index<elem_side_pairs.size(); ++face_index)
            {
                stk::mesh::Entity element = local_id_to_element_entity[elem_side_pairs[face_index].first];
                stk::mesh::EntityId face_global_id = face_index + 1 + offset_per_proc;
                stk::mesh::declare_element_side(bulkData, face_global_id, element, elem_side_pairs[face_index].second, &skin_part);
            }

            stk::mesh::Selector element_selector = bulkData.mesh_meta_data().locally_owned_part();
            bulkData.my_internal_modification_end_for_skin_mesh(stk::topology::FACE_RANK, stk::mesh::BulkData::MOD_END_SORT, element_selector, NULL);

            double elapsed_time = stk::wall_time() - wall_time_start;

            stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());

            std::vector<size_t> counts;
            stk::mesh::comm_mesh_counts(bulkData, counts);

            if (stk::parallel_machine_rank(comm) == 0)
            {
                std::cerr << "Element graph time: " << elapsed_time << std::endl;
                std::cerr << "Total # of elements: " << counts[stk::topology::ELEM_RANK] << std::endl;
            }
        }
    }
}

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

void add_element_side_pair_for_unused_sides(int64_t elementId, stk::topology topology, const std::vector<int64_t> &internal_sides, std::vector<std::pair<int64_t, int64_t> >& element_side_pairs)
{
    size_t num_sides = topology.num_sides();
    std::vector<int64_t> elem_sides;

    if (internal_sides.size() < num_sides)
    {
        elem_sides.assign(num_sides, -1);
        for(size_t j=0; j<internal_sides.size(); ++j)
        {
            int64_t zeroBasedSideId = internal_sides[j];
            elem_sides[zeroBasedSideId] = internal_sides[j];
        }

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

std::vector<std::pair<int64_t, int64_t> > skin_mesh(const std::vector<std::vector<int64_t> > &via_side,
   const std::vector<stk::topology> &element_topologies)
{
    std::vector<std::pair<int64_t, int64_t> > element_side_pairs;

    size_t num_elems = via_side.size();
    for(size_t i=0; i<num_elems; ++i)
    {
        const std::vector<int64_t>& internal_sides = via_side[i];
        add_element_side_pair_for_unused_sides(i, element_topologies[i], internal_sides, element_side_pairs);
    }
    return element_side_pairs;
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

stk::mesh::EntityVector get_elements_to_communicate(stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector elements_to_communicate;
    std::set<stk::mesh::Entity> element_set;
    const stk::mesh::BucketVector& shared_node_buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().globally_shared_part());
    for(size_t i=0; i<shared_node_buckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *shared_node_buckets[i];
        for(size_t node_index=0; node_index<bucket.size(); ++node_index)
        {
            stk::mesh::Entity node = bucket[node_index];
            const stk::mesh::Entity* elements = bulkData.begin_elements(node);
            unsigned num_elements = bulkData.num_elements(node);
            for(unsigned element_index=0; element_index<num_elements; ++element_index)
            {
                if (bulkData.bucket(elements[element_index]).owned())
                {
                    element_set.insert(elements[element_index]);
                }
            }
        }
    }
    elements_to_communicate.assign(element_set.begin(), element_set.end());
    return elements_to_communicate;
}

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements_to_communicate)
{
    for(size_t element_index=0; element_index<elements_to_communicate.size(); ++element_index)
    {
        stk::mesh::Entity elem = elements_to_communicate[element_index];
        stk::topology topology = bulkData.bucket(elem).topology();
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);
        unsigned num_sides = topology.num_sides();
        for(unsigned side_index=0; side_index<num_sides; ++side_index)
        {
            unsigned num_nodes_this_side = topology.side_topology(side_index).num_nodes();
            stk::mesh::EntityVector side_nodes(num_nodes_this_side);
            topology.side_nodes(elem_nodes, side_index, side_nodes.begin());

            std::vector<stk::mesh::EntityKey> side_node_entity_keys(num_nodes_this_side);
            for(size_t i=0; i<num_nodes_this_side; ++i)
            {
                side_node_entity_keys[i] = bulkData.entity_key(side_nodes[i]);
            }

            std::vector<int> sharing_procs;
            bulkData.shared_procs_intersection(side_node_entity_keys, sharing_procs);

            for(size_t proc_index=0; proc_index<sharing_procs.size(); ++proc_index)
            {
                comm.send_buffer(sharing_procs[proc_index]).pack<unsigned>(num_nodes_this_side);
                for(size_t i=0; i<num_nodes_this_side; ++i)
                {
                    comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityKey>(side_node_entity_keys[i]);
                }
            }
        }
    }
}

void add_possibly_connected_elements_to_graph_using_side_nodes(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides,
                                                                               const stk::mesh::EntityVector& side_nodes)
{
    stk::mesh::EntityVector elements;
    unsigned num_side_nodes = side_nodes.size();
    stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, num_side_nodes, side_nodes.data(), elements);
    for(size_t element_index=0; element_index<elements.size(); ++element_index)
    {
        stk::mesh::Entity elem = elements[element_index];
        stk::topology topology = bulkData.bucket(elem).topology();
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);
        unsigned num_sides = topology.num_sides();
        for(unsigned side_index=0; side_index<num_sides; ++side_index)
        {
            unsigned num_nodes_this_side = topology.side_topology(side_index).num_nodes();
            if (num_nodes_this_side == num_side_nodes)
            {
                stk::mesh::EntityVector side_nodes_this_side(num_nodes_this_side);
                topology.side_nodes(elem_nodes, side_index, side_nodes_this_side.begin());

                if (topology.side_topology(side_index).equivalent(side_nodes_this_side, side_nodes).first == true)
                {
                    unsigned local_elem_id = bulkData.local_id(elem);
                    elem_graph[local_elem_id].push_back(-1);
                    via_sides[local_elem_id].push_back(side_index);
                    break;
                }
            }
        }
    }
}

void fill_parallel_graph_with_test(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides)
{
    stk::mesh::EntityVector elements_to_communicate = get_elements_to_communicate(bulkData);
    ASSERT_EQ(1u, elements_to_communicate.size());

    if (bulkData.parallel_rank() == 0)
    {
        EXPECT_EQ(2u, bulkData.identifier(elements_to_communicate[0]));
    }
    else
    {
        EXPECT_EQ(3u, bulkData.identifier(elements_to_communicate[0]));
    }

    stk::CommSparse comm(bulkData.parallel());

    for(int phase=0; phase<2; ++phase)
    {
        pack_shared_side_nodes_of_elements(comm, bulkData, elements_to_communicate);

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    for(int proc_id=0; proc_id<bulkData.parallel_size(); ++proc_id)
    {
        if (proc_id != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                unsigned num_side_nodes = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(num_side_nodes);
                stk::mesh::EntityVector side_nodes(num_side_nodes);
                for(unsigned i=0; i<num_side_nodes; ++i)
                {
                    stk::mesh::EntityKey key;
                    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(key);
                    side_nodes[i] = bulkData.get_entity(key);
                }

                add_possibly_connected_elements_to_graph_using_side_nodes(bulkData, elem_graph, via_sides, side_nodes);
            }
        }
    }
}

void fill_parallel_graph(stk::mesh::BulkData& bulkData, std::vector<std::vector<int64_t> >& elem_graph, std::vector<std::vector<int64_t> >& via_sides)
{
    stk::mesh::EntityVector elements_to_communicate = get_elements_to_communicate(bulkData);

    stk::CommSparse comm(bulkData.parallel());

    for(int phase=0; phase<2; ++phase)
    {
        pack_shared_side_nodes_of_elements(comm, bulkData, elements_to_communicate);

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    for(int proc_id=0; proc_id<bulkData.parallel_size(); ++proc_id)
    {
        if (proc_id != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                unsigned num_side_nodes = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(num_side_nodes);
                stk::mesh::EntityVector side_nodes(num_side_nodes);
                for(unsigned i=0; i<num_side_nodes; ++i)
                {
                    stk::mesh::EntityKey key;
                    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(key);
                    side_nodes[i] = bulkData.get_entity(key);
                }

                add_possibly_connected_elements_to_graph_using_side_nodes(bulkData, elem_graph, via_sides, side_nodes);
            }
        }
    }
}

}
