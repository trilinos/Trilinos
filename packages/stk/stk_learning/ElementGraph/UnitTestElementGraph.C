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
#include <stk_mesh/base/CreateFaces.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace
{

struct sharing_info
{
    stk::mesh::Entity m_entity;
    int m_sharing_proc;
    int m_owner;
    sharing_info(stk::mesh::Entity entity, int sharing_proc, int owner) :
        m_entity(entity), m_sharing_proc(sharing_proc), m_owner(owner) {}
};

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

    bool my_modification_end_for_entity_creation(const std::vector<sharing_info>& shared_modified, modification_optimization opt = MOD_END_SORT)
    {
        if ( this->in_synchronized_state() ) { return false ; }

        ThrowAssertMsg(stk::mesh::impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

        if (parallel_size() > 1)
        {
            stk::mesh::PartVector shared_part, owned_part, empty;
            shared_part.push_back(&m_mesh_meta_data.globally_shared_part());
            owned_part.push_back(&m_mesh_meta_data.locally_owned_part());

            stk::mesh::EntityVector modified_entities(shared_modified.size());
            for(size_t i = 0; i < shared_modified.size(); ++i)
            {
                stk::mesh::Entity entity = shared_modified[i].m_entity;
                int sharing_proc = shared_modified[i].m_sharing_proc;
                entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
                int owning_proc = shared_modified[i].m_owner;
                const bool am_not_owner = this->internal_set_parallel_owner_rank_but_not_comm_lists(entity, owning_proc);
                if (am_not_owner)
                {
                    stk::mesh::EntityKey key = this->entity_key(entity);
                    internal_change_owner_in_comm_data(key, owning_proc);
                    internal_change_entity_parts(entity, shared_part /*add*/, owned_part /*remove*/);
                }
                else
                {
                    internal_change_entity_parts(entity, shared_part /*add*/, empty /*remove*/);
                }
                modified_entities[i] = entity;
            }

            std::sort(modified_entities.begin(), modified_entities.end(), stk::mesh::EntityLess(*this));
            stk::mesh::EntityVector::iterator iter = std::unique(modified_entities.begin(), modified_entities.end());
            modified_entities.resize(iter-modified_entities.begin());

            add_comm_list_entries_for_entities( modified_entities );

            internal_resolve_shared_membership();

            if ( this->get_automatic_aura_option() == AUTO_AURA)
            {
              this->resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(mesh_meta_data().side_rank(), mesh_meta_data().universal_part());
            }

            check_mesh_consistency();
        }
        else
        {
            std::vector<stk::mesh::Entity> modified_entities ;
            internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(mesh_meta_data().side_rank(), modified_entities);
        }

        this->internal_finish_modification_end(opt);

        return true ;
    }
};

typedef int64_t LocalId;
typedef int ProcId;
typedef int SideId;

struct parallel_info
{
    ProcId m_other_proc;
    SideId m_other_side_ord;
    int m_permutation;
    parallel_info(ProcId proc, SideId side_ord, int perm) :
        m_other_proc(proc), m_other_side_ord(side_ord), m_permutation(perm) {}
};

typedef std::pair<LocalId,SideId> ElementSidePair;
typedef std::map<std::pair<LocalId,stk::mesh::EntityId>, parallel_info > ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<SideId> > SidesForElementGraph;

int check_connectivity(const ElementGraph& elem_graph, const SidesForElementGraph &via_side,
        LocalId element_id1, LocalId element_id2);

std::vector<ElementSidePair>
skin_mesh(const SidesForElementGraph &via_side, const std::vector<stk::topology> &element_topologies);

void fill_parallel_graph(stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info);

void add_possibly_connected_elements_to_graph_using_side_nodes(stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, const stk::mesh::EntityVector& side_nodes, ParallelGraphInfo& parallel_graph_info,
        LocalId other_element, SideId other_side, ProcId other_proc);

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements_to_communicate);

stk::mesh::EntityVector get_elements_to_communicate(stk::mesh::BulkData& bulkData);

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity,
        std::vector<stk::topology>& element_topologies);

void add_element_side_pairs_for_unused_sides(LocalId elementId, stk::topology topology, const std::vector<SideId> &internal_sides,
        std::vector<ElementSidePair>& element_side_pairs);

void fill_graph(stk::mesh::BulkData& bulkData, ElementGraph& elem_graph, SidesForElementGraph& via_sides);

void create_faces_using_graph(BulkDataElementGraphTester& bulkData, stk::mesh::Part& part);

TEST(ElementGraph, check_graph_connectivity)
{
    // element0 --> element1 --> element2
    ElementGraph elem_graph = {
            {1},
            {0,2},
            {1}
    };

    SidesForElementGraph via_side = {
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
    ElementGraph elem_graph = {
            {1},
            {0,2},
            {1}
    };

    SidesForElementGraph via_side = {
            {4},
            {1,5},
            {3}
    };

    std::vector<stk::topology> element_topologies{
        stk::topology::HEXAHEDRON_8,
        stk::topology::HEXAHEDRON_8,
        stk::topology::HEXAHEDRON_8
    };

    std::vector<ElementSidePair> element_side_pairs = skin_mesh(via_side, element_topologies);

    std::vector<ElementSidePair>gold_element_side_pairs{
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
        std::vector<ElementSidePair >::iterator iter = std::find(element_side_pairs.begin(), element_side_pairs.end(), gold_element_side_pairs[i]);
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

        ElementGraph elem_graph(numElems);
        SidesForElementGraph via_sides(numElems);

        fill_graph(bulkData, elem_graph, via_sides);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        SideId left_side_id = 4;
        SideId right_side_id = 5;

        for(size_t i=0; i<elem_graph.size(); ++i)
        {
            const std::vector<LocalId>& conn_elements = elem_graph[i];
            if (i == 0)
            {
                ASSERT_EQ(1u, conn_elements.size());
                EXPECT_EQ(1, conn_elements[0]);
                EXPECT_EQ(right_side_id, via_sides[i][0]);
            }
            else if (i == elem_graph.size() - 1)
            {
                LocalId second_to_last_element_index = elem_graph.size() - 2;
                ASSERT_EQ(1u, conn_elements.size());
                EXPECT_EQ(second_to_last_element_index, conn_elements[0]);
                EXPECT_EQ(left_side_id, via_sides[i][0]);
            }
            else
            {
                ASSERT_EQ(2u, conn_elements.size());
                LocalId element_to_the_left = i-1;
                LocalId element_to_the_right = i+1;
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

        ElementGraph elem_graph(numLocallyOwnedElems);
        SidesForElementGraph via_sides(numLocallyOwnedElems);

        ParallelGraphInfo parallel_graph_info;
        fill_parallel_graph(bulkData, elem_graph, via_sides, parallel_graph_info);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        LocalId element_to_test_local_id = std::numeric_limits<LocalId>::max();
        SideId side_id = -1;
        SideId left_side_id = 4;
        SideId right_side_id = 5;

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
            if (static_cast<LocalId>(i) != element_to_test_local_id)
            {
                EXPECT_TRUE(elem_graph[i].empty());
                EXPECT_TRUE(via_sides[i].empty());
            }
            else
            {
                ASSERT_TRUE(elem_graph[i].size()==1);
                ASSERT_TRUE(via_sides[i].size()==1);
                EXPECT_GE(-1, elem_graph[i][0]);
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

        ElementGraph elem_graph(numElems);
        SidesForElementGraph via_sides(numElems);

        fill_graph(bulkData, elem_graph, via_sides);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<ElementSidePair> elem_side_pairs = skin_mesh(via_sides, element_topologies);

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
            SideId left_side_id = 4;
            SideId right_side_id = 5;

            for(size_t i=0; i<elem_graph.size(); ++i)
            {
                const std::vector<LocalId>& conn_elements = elem_graph[i];
                if (i == 0)
                {
                    ASSERT_EQ(1u, conn_elements.size());
                    EXPECT_EQ(1, conn_elements[0]);
                    EXPECT_EQ(right_side_id, via_sides[i][0]);
                }
                else if (i == elem_graph.size() - 1)
                {
                    LocalId second_to_last_element_index = elem_graph.size() - 2;
                    ASSERT_EQ(1u, conn_elements.size());
                    EXPECT_EQ(second_to_last_element_index, conn_elements[0]);
                    EXPECT_EQ(left_side_id, via_sides[i][0]);
                }
                else
                {
                    ASSERT_EQ(2u, conn_elements.size());
                    LocalId element_to_the_left = i-1;
                    LocalId element_to_the_right = i+1;
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

        ElementGraph elem_graph(num_locally_owned_elems);
        SidesForElementGraph via_sides(num_locally_owned_elems);

        fill_graph(bulkData, elem_graph, via_sides);
        if (stk::parallel_machine_size(comm) > 1)
        {
            ParallelGraphInfo parallel_graph_info;
            fill_parallel_graph(bulkData, elem_graph, via_sides, parallel_graph_info);
        }

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<ElementSidePair> elem_side_pairs = skin_mesh(via_sides, element_topologies);

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

stk::mesh::Entity get_element_side(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal)
{
    stk::mesh::Entity side = stk::mesh::Entity();
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    unsigned elem_num_faces = bulkData.num_connectivity(element, side_rank);
    const stk::mesh::Entity * elem_sides = bulkData.begin(element, side_rank);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulkData.begin_ordinals(element, side_rank);
    for (unsigned i=0 ; i<elem_num_faces ; ++i)
    {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal))
        {
            side = elem_sides[i];
            break;
        }
    }

    return side;
}

bool does_element_side_exist(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal)
{
    stk::mesh::Entity side = stk::mesh::Entity();
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    unsigned elem_num_faces = bulkData.num_connectivity(element, side_rank);
    const stk::mesh::Entity * elem_sides = bulkData.begin(element, side_rank);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulkData.begin_ordinals(element, side_rank);
    for (unsigned i=0 ; i<elem_num_faces ; ++i)
    {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal))
        {
            side = elem_sides[i];
            break;
        }
    }

    return bulkData.is_valid(side);
}

bool does_side_exist_with_different_permutation(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::ConnectivityOrdinal side_ordinal, stk::mesh::Permutation perm)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    unsigned elem_num_faces = bulkData.num_connectivity(element, side_rank);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulkData.begin_ordinals(element, side_rank);
    const stk::mesh::Permutation * elem_perm_it = bulkData.begin_permutations(element, side_rank);

    for (unsigned i=0 ; i<elem_num_faces ; ++i)
    {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal))
        {
            if (perm != elem_perm_it[i])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    return false;
}

bool is_id_already_in_use_locally(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
{
    stk::mesh::Entity entity = bulkData.get_entity(rank, id);
    return bulkData.is_valid(entity);
}

stk::mesh::Entity connect_face_to_element(stk::mesh::BulkData& bulkData, stk::mesh::Entity element,
        stk::mesh::EntityId side_global_id, stk::mesh::ConnectivityOrdinal side_ordinal,
        stk::mesh::Permutation side_permutation, stk::mesh::Part& part)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    stk::mesh::Entity side = bulkData.declare_entity(side_rank, side_global_id, part);

    // connect element to side
    bulkData.declare_relation(element, side, side_ordinal, side_permutation);

    // connect side to nodes
    const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
    stk::topology side_top = bulkData.bucket(element).topology().side_topology(side_ordinal);
    stk::mesh::EntityVector side_nodes(side_top.num_nodes());
    bulkData.bucket(element).topology().side_nodes(nodes, side_ordinal, side_nodes.begin());
    stk::mesh::EntityVector permuted_side_nodes(side_top.num_nodes());
    side_top.permutation_nodes(side_nodes, side_permutation, permuted_side_nodes.begin());
    for(size_t i=0;i<permuted_side_nodes.size();++i)
    {
        bulkData.declare_relation(side, permuted_side_nodes[i], i);
    }

    return side;
}

TEST(ElementSide, get_or_create_element_side_with_permutation)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) == 1)
    {
        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(new_faces_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x3", bulkData, comm);

        //////////////// Make first side

        bulkData.modification_begin();

        //get_or_create_element_side_with_permutation(bulkData, element);
        stk::mesh::Entity element1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        stk::mesh::EntityId side_global_id = 11;
        stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
        stk::mesh::Permutation side_permutation = static_cast<stk::mesh::Permutation>(4);
        stk::mesh::ConnectivityOrdinal side_ordinal = static_cast<stk::mesh::ConnectivityOrdinal>(1);

        EXPECT_FALSE(is_id_already_in_use_locally(bulkData, side_rank, side_global_id));
        EXPECT_FALSE(does_side_exist_with_different_permutation(bulkData, element1, side_ordinal, side_permutation));
        EXPECT_FALSE(does_element_side_exist(bulkData, element1, side_ordinal));

        connect_face_to_element(bulkData, element1, side_global_id, side_ordinal, side_permutation, new_faces_part);

        bulkData.modification_end();

        stk::mesh::Entity side1 = get_element_side(bulkData, element1, side_ordinal);
        EXPECT_TRUE(bulkData.is_valid(side1));

        EXPECT_TRUE(is_id_already_in_use_locally(bulkData, side_rank, side_global_id));

        stk::mesh::Permutation side_permutation1 = static_cast<stk::mesh::Permutation>(0);
        EXPECT_TRUE(does_side_exist_with_different_permutation(bulkData, element1, side_ordinal, side_permutation1));

        size_t num_sides = stk::mesh::count_selected_entities(new_faces_part, bulkData.buckets(side_rank));
        EXPECT_EQ(1u, num_sides);
    }
}

void test_parallel_graph_info(const ElementGraph& elem_graph, const ParallelGraphInfo& parallel_graph_info,
        LocalId this_element, LocalId other_element, ProcId other_proc, SideId other_side_ord, int permutation)
{
    ParallelGraphInfo::const_iterator iter1 = parallel_graph_info.find(std::make_pair(this_element, other_element));
    ASSERT_TRUE(iter1 != parallel_graph_info.end());

    for(size_t i=0;i<elem_graph.size();++i)
    {
        const std::vector<LocalId>& conn_elements = elem_graph[i];
        for(size_t j=0;j<conn_elements.size();++j)
        {
            if(conn_elements[j]==-1*other_element && static_cast<LocalId>(i) == this_element)
            {
                ParallelGraphInfo::const_iterator iter = parallel_graph_info.find(std::make_pair(this_element, other_element));

                ASSERT_TRUE(iter != parallel_graph_info.end());
                EXPECT_EQ(other_proc, iter->second.m_other_proc);
                EXPECT_EQ(other_side_ord, iter->second.m_other_side_ord);
                EXPECT_EQ(permutation, iter->second.m_permutation);
            }
        }
    }
}

TEST(ElementGraph, test_parallel_graph_info_data_structure)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        ElementGraph elem_graph {
                {1},
                {0,-3},
        };

        SidesForElementGraph via_side {
                {4},
                {1,5},
        };


        ParallelGraphInfo parallel_graph_info;
        ProcId other_proc = 1;
        SideId other_side_ord = 2;
        LocalId local_element = 1;
        LocalId other_element = 3;
        int permutation = 0;

        parallel_graph_info.insert(std::make_pair(std::make_pair(local_element, other_element), parallel_info(other_proc, other_side_ord, permutation)));

        size_t num_elems_this_proc = elem_graph.size();
        EXPECT_EQ(2u, num_elems_this_proc);

        test_parallel_graph_info(elem_graph, parallel_graph_info, local_element, other_element, other_proc, other_side_ord, permutation);
    }
}

TEST(ElementGraph, test_parallel_graph_info_with_parallel_element_graph)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];

        stk::mesh::EntityVector local_id_to_element_entity(numLocallyOwnedElems, 0);
        std::vector<stk::topology> element_topologies(numLocallyOwnedElems);
        set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        ElementGraph elem_graph(numLocallyOwnedElems);
        SidesForElementGraph via_sides(numLocallyOwnedElems);

        ParallelGraphInfo parallel_graph_info;
        fill_parallel_graph(bulkData, elem_graph, via_sides, parallel_graph_info);

        if(stk::parallel_machine_rank(comm)==0)
        {
            LocalId local_element = 1;
            LocalId other_element = 3;
            ProcId other_proc = 1;
            SideId other_side_ord = 4; // 4 left, 5 right
            int permutation = 4;

            test_parallel_graph_info(elem_graph, parallel_graph_info, local_element, other_element, other_proc, other_side_ord, permutation);
        }
        else
        {
            LocalId local_element = 0;
            LocalId other_element = 2;
            ProcId other_proc = 0;
            SideId other_side_ord = 5; // 4 left, 5 right
            int permutation = 4;

            test_parallel_graph_info(elem_graph, parallel_graph_info, local_element, other_element, other_proc, other_side_ord, permutation);
        }
    }
}

int get_element_face_multiplier()
{
    return 10;
}

template<class T>
void print_graph(const std::string &title, int proc_id, T& elem_graph)
{
    std::ostringstream os;

    os << title << " for processor " << proc_id << std::endl;
    for (size_t i=0;i<elem_graph.size();++i)
    {
        os << "Element " << i << "::\t";
        for (size_t j=0;j<elem_graph[i].size();++j)
        {
            os << elem_graph[i][j] << "\t";
        }
        os << std::endl;
    }
    std::cerr << os.str();
}



TEST(ElementGraph, create_faces_using_element_graph_parallel)
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
        stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(new_faces_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        create_faces_using_graph(bulkData, new_faces_part);

        const stk::mesh::BucketVector& sharedNodeBuckets = bulkData.get_buckets(stk::topology::NODE_RANK, meta.globally_shared_part());
        for(size_t bucket_index=0; bucket_index<sharedNodeBuckets.size(); ++bucket_index)
        {
            const stk::mesh::Bucket& bucket = *sharedNodeBuckets[bucket_index];
            EXPECT_TRUE(bucket.member(new_faces_part));
        }

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after create-faces");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        stk::mesh::EntityRank side_rank = meta.side_rank();
        unsigned num_faces = entity_counts[side_rank];
        EXPECT_EQ(21u, num_faces);

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

TEST(ElementGraph, create_faces_using_element_graph_parallel_block_membership)
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
        stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);

        stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
        EXPECT_EQ(stk::topology::ELEM_RANK, block_2.primary_entity_rank());

        stk::io::put_io_part_attribute(new_faces_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        stk::mesh::Part& block_1 = *meta.get_part("block_1");

        bulkData.modification_begin();

        if (bulkData.parallel_rank() == 1)
        {
            //Move elem 3 from block_1 to block_2 so that the boundary between block_1 and block_2
            //will coincide with the proc boundary, and the shared face between elems 2 & 3
            //should have both block_1 and block_2.
            stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            ASSERT_TRUE(bulkData.is_valid(elem3));
            stk::mesh::PartVector add_parts(1, &block_2);
            stk::mesh::PartVector rem_parts(1, &block_1);
            bulkData.change_entity_parts(elem3, add_parts, rem_parts);
        }

        bulkData.modification_end();

        create_faces_using_graph(bulkData, new_faces_part);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after create-faces");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
        unsigned num_faces = entity_counts[side_rank];
        EXPECT_EQ(21u, num_faces);

        const stk::mesh::BucketVector& sharedFaceBuckets = bulkData.get_buckets(stk::topology::FACE_RANK, meta.globally_shared_part());
        if (bulkData.parallel_size() == 2)
        {
            ASSERT_EQ(1u, sharedFaceBuckets.size());

            const stk::mesh::Bucket& bucket = *sharedFaceBuckets[0];
            ASSERT_EQ(1u, bucket.size());
            EXPECT_TRUE(bucket.member(block_2));
            EXPECT_TRUE(bucket.member(block_1));
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
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        {
            double wall_time_start = stk::wall_time();

            unsigned num_locally_owned_elems = stk::mesh::count_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK));

            stk::mesh::EntityVector local_id_to_element_entity(num_locally_owned_elems, 0);
            std::vector<stk::topology> element_topologies(num_locally_owned_elems);
            set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

            ElementGraph elem_graph(num_locally_owned_elems);
            SidesForElementGraph via_sides(num_locally_owned_elems);

            fill_graph(bulkData, elem_graph, via_sides);

            ParallelGraphInfo parallel_graph_info;
            fill_parallel_graph(bulkData, elem_graph, via_sides, parallel_graph_info);

            std::vector<ElementSidePair> elem_side_pairs = skin_mesh(via_sides, element_topologies);

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

void create_faces_using_graph(BulkDataElementGraphTester& bulkData, stk::mesh::Part& part)
{
    double wall_time_start = stk::wall_time();

    std::vector<unsigned> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numElems = counts[stk::topology::ELEM_RANK];

    stk::mesh::EntityVector local_id_to_element_entity(numElems, 0);
    std::vector<stk::topology> element_topologies(numElems);
    set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

    ElementGraph elem_graph(numElems);
    SidesForElementGraph via_sides(numElems);

    fill_graph(bulkData, elem_graph, via_sides);

    ParallelGraphInfo parallel_graph_info;
    fill_parallel_graph(bulkData, elem_graph, via_sides, parallel_graph_info);

    double graph_time = stk::wall_time() - wall_time_start;
    wall_time_start = stk::wall_time();

    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    bulkData.modification_begin();

    std::vector<sharing_info> shared_modified;

    for(size_t i = 0; i < elem_graph.size(); ++i)
    {
        const std::vector<LocalId>& connected_elements = elem_graph[i];
        stk::mesh::Entity element1 = local_id_to_element_entity[i];

        LocalId this_element = i;

        for(size_t j = 0; j < connected_elements.size(); ++j)
        {
            if(this_element < connected_elements[j] && connected_elements[j] > 0)
            {
                stk::mesh::EntityId face_global_id = get_element_face_multiplier() * bulkData.identifier(element1) + via_sides[i][j];
                if ( is_id_already_in_use_locally(bulkData, side_rank, face_global_id) )
                {

                }
                stk::mesh::Entity face = stk::mesh::impl::get_or_create_face_at_element_side(bulkData, element1, via_sides[i][j],
                        face_global_id, part);

                const stk::mesh::Entity* side_nodes = bulkData.begin_nodes(face);
                unsigned num_side_nodes = bulkData.num_nodes(face);
                stk::mesh::EntityVector side_nodes_vec(side_nodes, side_nodes + num_side_nodes);

                stk::mesh::Entity element2 = local_id_to_element_entity[connected_elements[j]];
                std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm = stk::mesh::get_ordinal_and_permutation(bulkData, element2, stk::topology::FACE_RANK, side_nodes_vec);
                bulkData.declare_relation(element2, face, ord_and_perm.first, ord_and_perm.second);
            }
            else if(connected_elements[j] < 0)
            {
                LocalId other_element = -1 * connected_elements[j];
                ParallelGraphInfo::const_iterator iter = parallel_graph_info.find(std::make_pair(this_element, other_element));
                ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
                ProcId other_proc = iter->second.m_other_proc;
                SideId other_side = iter->second.m_other_side_ord;

                ProcId this_proc = bulkData.parallel_rank();
                ProcId owning_proc = this_proc < other_proc ? this_proc : other_proc;

                stk::mesh::EntityId face_global_id = 0;
                stk::mesh::Permutation perm;
                if(owning_proc == this_proc)
                {
                    stk::mesh::EntityId id = bulkData.identifier(element1);
                    face_global_id = get_element_face_multiplier() * id + via_sides[i][j];
                    perm = static_cast<stk::mesh::Permutation>(0);
                }
                else
                {
                    face_global_id = get_element_face_multiplier() * other_element + other_side;
                    perm = static_cast<stk::mesh::Permutation>(iter->second.m_permutation);
                }

                stk::mesh::ConnectivityOrdinal side_ord = static_cast<stk::mesh::ConnectivityOrdinal>(via_sides[i][j]);

                std::string msg = "Program error. Contact sierra-help@sandia.gov for support.";

                ThrowRequireMsg(!is_id_already_in_use_locally(bulkData, side_rank, face_global_id), msg);
                ThrowRequireMsg(!does_side_exist_with_different_permutation(bulkData, element1, side_ord, perm), msg);
                ThrowRequireMsg(!does_element_side_exist(bulkData, element1, side_ord), msg);

                stk::mesh::Entity face = connect_face_to_element(bulkData, element1, face_global_id, side_ord, perm, part);

                shared_modified.push_back(sharing_info(face, other_proc, owning_proc));
            }
        }

        std::vector<ElementSidePair> element_side_pairs;
        add_element_side_pairs_for_unused_sides(i, element_topologies[i], via_sides[i], element_side_pairs);

        for(size_t j = 0; j < element_side_pairs.size(); j++)
        {
            stk::mesh::EntityId face_global_id = get_element_face_multiplier() * bulkData.identifier(element1) + element_side_pairs[j].second;
            stk::mesh::impl::get_or_create_face_at_element_side(bulkData, element1, element_side_pairs[j].second,
                    face_global_id, part);

        }
    }

    double start_mod_end = stk::wall_time();
    bulkData.my_modification_end_for_entity_creation(shared_modified);
    double mod_end_time = stk::wall_time() - start_mod_end;

    double create_faces_time = stk::wall_time() - wall_time_start;

    if(bulkData.parallel_rank() == 0)
    {
        std::cerr << "graph time: " << graph_time << std::endl;
        std::cerr << "create faces time: " << create_faces_time << std::endl;
        std::cerr << "mod end time: " << mod_end_time << std::endl;
    }
}

TEST(ElementGraph, compare_performance_create_faces)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

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
        stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(faces_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        {
            double wall_time_start = stk::wall_time();

            stk::mesh::PartVector parts(1, &faces_part);
            stk::mesh::create_faces(bulkData);

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
        bool force_no_induce = true;
        stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4, force_no_induce);
        stk::io::put_io_part_attribute(faces_part);
        BulkDataElementGraphTester bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);

        {
            double wall_time_start = stk::wall_time();

            create_faces_using_graph(bulkData, faces_part);

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

int check_connectivity(const ElementGraph& elem_graph, const SidesForElementGraph &via_side,
        LocalId element_id1, LocalId element_id2)
{
    int side=-1;

    if(element_id1 >=0 && element_id1 <=2 && element_id2 >=0 && element_id2 <=2)
    {
        const std::vector<LocalId>& conn_elements = elem_graph[element_id1];

        std::vector<LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), element_id2);
        if ( iter != conn_elements.end() )
        {
            int64_t index = iter - conn_elements.begin();
            side = via_side[element_id1][index];
        }
    }

    return side;
}

void add_element_side_pairs_for_unused_sides(LocalId elementId, stk::topology topology, const std::vector<SideId> &internal_sides,
        std::vector<ElementSidePair>& element_side_pairs)
{
    size_t num_sides = topology.num_sides();
    std::vector<SideId> elem_sides;

    if (internal_sides.size() < num_sides)
    {
        elem_sides.assign(num_sides, -1);
        for(size_t j=0; j<internal_sides.size(); ++j)
        {
            SideId sideId = internal_sides[j];
            elem_sides[sideId] = internal_sides[j];
        }

        for(size_t j=0; j<num_sides; ++j)
        {
            if (elem_sides[j] == -1)
            {
                SideId sideId = j;
                element_side_pairs.push_back(std::make_pair(elementId, sideId));
            }
        }
    }
}

std::vector<ElementSidePair>
skin_mesh(const SidesForElementGraph &via_side, const std::vector<stk::topology> &element_topologies)
{
    std::vector<ElementSidePair> element_side_pairs;

    size_t num_elems = via_side.size();
    for(size_t i=0; i<num_elems; ++i)
    {
        const std::vector<SideId>& internal_sides = via_side[i];
        add_element_side_pairs_for_unused_sides(i, element_topologies[i], internal_sides, element_side_pairs);
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

void fill_graph(stk::mesh::BulkData& bulkData, ElementGraph& elem_graph, SidesForElementGraph& via_sides)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        stk::topology topology = bucket.topology();
        int num_sides = topology.num_sides();
        std::vector<ElementSidePair> elem_side_pairs;
        stk::mesh::EntityVector side_nodes;
        stk::mesh::EntityVector connected_elements;
        for(size_t j=0; j<bucket.size(); ++j)
        {
            LocalId local_elem_id = bulkData.local_id(bucket[j]);
            const stk::mesh::Entity* elem_nodes = bucket.begin_nodes(j);
            elem_side_pairs.clear();
            for(SideId side_index=0; side_index<num_sides; ++side_index)
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
            std::vector<ElementSidePair>::iterator new_end = std::unique(elem_side_pairs.begin(), elem_side_pairs.end());
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
        stk::mesh::EntityId element_id = bulkData.identifier(elem);

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
                comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityId>(element_id);
                comm.send_buffer(sharing_procs[proc_index]).pack<unsigned>(side_index);
                comm.send_buffer(sharing_procs[proc_index]).pack<unsigned>(num_nodes_this_side);
                for(size_t i=0; i<num_nodes_this_side; ++i)
                {
                    comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityKey>(side_node_entity_keys[i]);
                }
            }
        }
    }
}

void add_possibly_connected_elements_to_graph_using_side_nodes(stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, const stk::mesh::EntityVector& side_nodes, ParallelGraphInfo& parallel_graph_info,
        LocalId other_element, SideId other_side, ProcId other_proc)
{
    stk::mesh::EntityVector elements;
    unsigned num_side_nodes = side_nodes.size();
    stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, num_side_nodes, side_nodes.data(), elements);
    for(size_t element_index=0; element_index<elements.size(); ++element_index)
    {
        stk::mesh::Entity elem = elements[element_index];
        stk::topology topology = bulkData.bucket(elem).topology();
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);
        int num_sides = topology.num_sides();
        for(SideId side_index=0; side_index<num_sides; ++side_index)
        {
            unsigned num_nodes_this_side = topology.side_topology(side_index).num_nodes();
            if (num_nodes_this_side == num_side_nodes)
            {
                stk::mesh::EntityVector side_nodes_this_side(num_nodes_this_side);
                topology.side_nodes(elem_nodes, side_index, side_nodes_this_side.begin());

                std::pair<bool,unsigned> result = topology.side_topology(side_index).equivalent(side_nodes_this_side, side_nodes);
                if (result.first == true)
                {
                    LocalId local_elem_id = bulkData.local_id(elem);
                    elem_graph[local_elem_id].push_back(-1*other_element);
                    via_sides[local_elem_id].push_back(side_index);

                    parallel_graph_info.insert(std::make_pair(std::make_pair(local_elem_id, other_element), parallel_info(other_proc, other_side, result.second)));

                    break;
                }
            }
        }
    }
}

void fill_parallel_graph(stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info)
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
                stk::mesh::EntityId element_id;
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(element_id);
                unsigned side_index = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(side_index);
                unsigned num_side_nodes = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(num_side_nodes);
                stk::mesh::EntityVector side_nodes(num_side_nodes);
                for(unsigned i=0; i<num_side_nodes; ++i)
                {
                    stk::mesh::EntityKey key;
                    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(key);
                    side_nodes[i] = bulkData.get_entity(key);
                }

                add_possibly_connected_elements_to_graph_using_side_nodes(bulkData, elem_graph, via_sides, side_nodes,
                        parallel_graph_info, element_id, side_index, proc_id);
            }
        }
    }
}

}
