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

    size_t num_entity_keys() const
    {
        return m_entity_keys.size();
    }

};

typedef int64_t LocalId;
typedef int ProcId;
typedef int SideId;

//BeginDocExample3
struct parallel_info
{
    ProcId m_other_proc;
    SideId m_other_side_ord;
    int m_permutation;
    bool m_in_part;
    parallel_info(ProcId proc, SideId side_ord, int perm) :
        m_other_proc(proc), m_other_side_ord(side_ord), m_permutation(perm), m_in_part(true) {}
};
//EndDocExample3

typedef std::pair<LocalId,SideId> ElementSidePair;
typedef std::map<std::pair<LocalId,stk::mesh::EntityId>, parallel_info > ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<SideId> > SidesForElementGraph;

int check_connectivity(const ElementGraph& elem_graph, const SidesForElementGraph &via_side,
        LocalId element_id1, LocalId element_id2);

std::vector<ElementSidePair>
skin_mesh(const SidesForElementGraph &via_side, const std::vector<stk::topology> &element_topologies);

void fill_parallel_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
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

void fill_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph, SidesForElementGraph& via_sides);

void create_faces_using_graph(BulkDataElementGraphTester& bulkData, stk::mesh::Part& part);

SideId get_side_from_element1_to_element2(const ElementGraph& elem_graph, const SidesForElementGraph &via_side, LocalId element_id1, LocalId element_id2);

class ElemElemGraph
{
public:
    ElemElemGraph(stk::mesh::BulkData& bulkData) : m_bulk_data(bulkData)
    {
        size_data_members();

        set_local_ids_and_fill_element_entities_and_topologies(m_bulk_data, local_id_to_element_entity, element_topologies);
        fill_graph(m_bulk_data, elem_graph, via_sides);
        fill_parallel_graph(m_bulk_data, elem_graph, via_sides, parallel_graph_info);
    }

    ~ElemElemGraph() {}

    size_t get_num_connected_elems(stk::mesh::Entity local_element) const
    {
        LocalId local_id = get_local_element_id(local_element);
        return elem_graph[local_id].size();
    }

    bool is_connected_elem_locally_owned(stk::mesh::Entity local_element, size_t index_conn_elem) const
    {
        LocalId local_id = get_local_element_id(local_element);
        return elem_graph[local_id][index_conn_elem] >= 0;
    }

    SideId get_side_id_to_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
    {
        LocalId local_id = get_local_element_id(local_element);
        return via_sides[local_id][index_conn_elem];
    }

    stk::mesh::Entity get_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
    {
        LocalId local_id = get_local_element_id(local_element);
        LocalId other_element_id = elem_graph[local_id][index_conn_elem];
        return local_id_to_element_entity[other_element_id];
    }

    stk::mesh::EntityId get_entity_id_of_remote_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
    {
        ThrowRequireMsg(!is_connected_elem_locally_owned(local_element, index_conn_elem) , "Program error. Contact sierra-help@sandia.gov for support.");
        LocalId local_id = get_local_element_id(local_element);
        stk::mesh::EntityId id = -elem_graph[local_id][index_conn_elem];
        return id;
    }

    ProcId get_owning_proc_id_of_remote_element(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
    {
        LocalId local_id = get_local_element_id(local_element);
        ParallelGraphInfo::const_iterator iter = parallel_graph_info.find(std::make_pair(local_id, other_element_id));
        ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
        ProcId other_proc = iter->second.m_other_proc;
        return other_proc;
    }

    SideId get_side_from_element1_to_remote_element2(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
    {
        LocalId remote_element_local_id = -other_element_id;
        LocalId element1_local_id = get_local_element_id(local_element);

        SideId side = -1;
        const std::vector<LocalId>& conn_elements = elem_graph[element1_local_id];

        std::vector<LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), remote_element_local_id);
        if ( iter != conn_elements.end() )
        {
            int64_t index = iter - conn_elements.begin();
            side = via_sides[element1_local_id][index];
        }
        return side;
    }

    SideId get_side_from_element1_to_locally_owned_element2(stk::mesh::Entity local_element, stk::mesh::Entity other_element) const
    {
        LocalId other_element_id = get_local_element_id(other_element);
        LocalId element1_local_id = get_local_element_id(local_element);

        SideId side = -1;
        const std::vector<LocalId>& conn_elements = elem_graph[element1_local_id];

        std::vector<LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), other_element_id);
        if ( iter != conn_elements.end() )
        {
            int64_t index = iter - conn_elements.begin();
            side = via_sides[element1_local_id][index];
        }
        return side;
    }

    parallel_info& get_parallel_edge_info(stk::mesh::Entity element, stk::mesh::EntityId remote_id)
    {
        LocalId this_elem_local_id = get_local_element_id(element);

        ParallelGraphInfo::iterator iter = parallel_graph_info.find(std::make_pair(this_elem_local_id, remote_id));
        ThrowRequireMsg( iter != parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
        return iter->second;
    }

private:

    LocalId get_local_element_id(stk::mesh::Entity local_element) const
    {
        ThrowRequireMsg(m_bulk_data.is_valid(local_element), "Program error. Contact sierra-help@sandia.gov for support.");
        LocalId local_id = m_bulk_data.local_id(local_element);
        return local_id;
    }

    void size_data_members()
    {
        std::vector<unsigned> counts;
        stk::mesh::count_entities(m_bulk_data.mesh_meta_data().locally_owned_part(), m_bulk_data, counts);
        int numElems = counts[stk::topology::ELEM_RANK];

        elem_graph.resize(numElems);
        via_sides.resize(numElems);
        local_id_to_element_entity.resize(numElems, 0);
        element_topologies.resize(numElems);
    }

    stk::mesh::BulkData &m_bulk_data;
    ElementGraph elem_graph;
    SidesForElementGraph via_sides;
    ParallelGraphInfo parallel_graph_info;
    stk::mesh::EntityVector local_id_to_element_entity;
    std::vector<stk::topology> element_topologies;
};

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

    //BeginDocExample4

    bulkData.modification_begin();

    std::vector<sharing_info> shared_modified;

    for(size_t i = 0; i < elem_graph.size(); ++i)
    {
        const std::vector<LocalId>& connected_elements = elem_graph[i];
        stk::mesh::Entity element1 = local_id_to_element_entity[i];

        LocalId this_element = i;

        for(size_t j = 0; j < connected_elements.size(); ++j)
        {
            if(this_element < connected_elements[j] && connected_elements[j] >= 0)
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

    //EndDocExample4

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

stk::mesh::EntityVector get_killed_elements(stk::mesh::BulkData& bulkData, const int killValue, const stk::mesh::Part& inactive)
{
    stk::mesh::EntityVector killedElements;
    const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);
    for(size_t b = 0; b < buckets.size(); ++b)
    {
        const stk::mesh::Bucket &bucket = *buckets[b];
        if(bucket.owned() && !bucket.member(inactive))
        {
            for(size_t e = 0; e < bucket.size(); ++e)
            {
                stk::mesh::Entity entity = bucket[e];
                bool should_element_be_killed = bulkData.identifier(entity) < static_cast<stk::mesh::EntityId>(killValue);
                if(!bulkData.bucket(entity).member(inactive) && should_element_be_killed == true)
                {
                    killedElements.push_back(bucket[e]);
                }
            }
        }
    }
    return killedElements;
}

void move_killled_elements_to_part(stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& killedElements, stk::mesh::Part& block_1, stk::mesh::Part& inactive)
{
    std::vector<stk::mesh::PartVector> add_parts(killedElements.size());
    std::vector<stk::mesh::PartVector> rm_parts(killedElements.size());

    stk::mesh::PartVector add_part_vec(1, &inactive);
    stk::mesh::PartVector rm_part_vec(1, &block_1);

    for (size_t j=0;j<killedElements.size();++j)
    {
        add_parts[j] = add_part_vec;
        rm_parts[j] = rm_part_vec;
    }

    bulkData.batch_change_entity_parts(killedElements, add_parts, rm_parts);
}

struct graphEdgeProc
{
    stk::mesh::EntityId m_localElementId;
    stk::mesh::EntityId m_remoteElementId;
    int m_proc_id;
    graphEdgeProc(const stk::mesh::EntityId& localElementId, const stk::mesh::EntityId &remoteElementId, int proc_id) :
        m_localElementId(localElementId), m_remoteElementId(remoteElementId), m_proc_id(proc_id) {}
};

std::vector<graphEdgeProc> get_elements_to_communicate(stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &killedElements,
        const ElemElemGraph& elem_graph)
{
    std::vector<graphEdgeProc> elements_to_comm;

    for(size_t i=0;i<killedElements.size();++i)
    {
        stk::mesh::Entity this_elem_entity = killedElements[i];
        for(size_t j=0;j<elem_graph.get_num_connected_elems(this_elem_entity);++j)
        {
            if(!elem_graph.is_connected_elem_locally_owned(this_elem_entity, j))
            {
                stk::mesh::EntityId other_element_id = elem_graph.get_entity_id_of_remote_element(this_elem_entity,j);
                ProcId other_proc = elem_graph.get_owning_proc_id_of_remote_element(this_elem_entity, other_element_id);
                elements_to_comm.push_back(graphEdgeProc(bulkData.identifier(this_elem_entity), other_element_id, other_proc));
            }
        }
    }

    return elements_to_comm;
}

stk::mesh::EntityId get_face_global_id(stk::mesh::BulkData& bulkData, const ElementGraph& elem_graph, const SidesForElementGraph &via_sides,
        const stk::mesh::EntityVector local_id_to_element_entity, LocalId element1_local_id, LocalId element2_local_id)
{
    stk::mesh::EntityId face_global_id;

    stk::mesh::Entity element1_entity = local_id_to_element_entity[element1_local_id];
    if(element1_local_id < element2_local_id)
    {
        SideId side_id = get_side_from_element1_to_element2(elem_graph, via_sides, element1_local_id, element2_local_id);
        EXPECT_TRUE(side_id != -1);
        face_global_id = get_element_face_multiplier() * bulkData.identifier(element1_entity) + side_id;
    }
    else
    {
        SideId side_id = get_side_from_element1_to_element2(elem_graph, via_sides, element2_local_id, element1_local_id);
        EXPECT_TRUE(side_id != -1);
        stk::mesh::Entity element2_entity = local_id_to_element_entity[element2_local_id];
        face_global_id = get_element_face_multiplier() * bulkData.identifier(element2_entity) + side_id;
    }

    return face_global_id;
}

void communicate_killed_entities(stk::mesh::BulkData& bulkData, const std::vector<graphEdgeProc>& elements_to_comm,
        std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges)
{
    stk::CommSparse comm(bulkData.parallel());
    for(int iphase=0;iphase<2;++iphase)
    {
        for(size_t i=0;i<elements_to_comm.size();++i)
        {
            int remote_proc = elements_to_comm[i].m_proc_id;
            stk::mesh::EntityId localId = elements_to_comm[i].m_localElementId;
            stk::mesh::EntityId remoteId = elements_to_comm[i].m_remoteElementId;

            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(localId);
            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(remoteId);
        }
        if(iphase==0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    for(int i=0;i<bulkData.parallel_size();++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            stk::mesh::EntityId remoteId;
            stk::mesh::EntityId localId;
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(remoteId);
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(localId);
            remote_edges.push_back(std::make_pair(localId, remoteId));
        }
    }
}

stk::mesh::EntityId get_face_global_id(const stk::mesh::BulkData &bulkData, const ElemElemGraph& elementGraph, stk::mesh::Entity element1, stk::mesh::Entity element2,
        SideId element1_side_id)
{
    stk::mesh::EntityId element1_global_id = bulkData.identifier(element1);
    stk::mesh::EntityId element2_global_id = bulkData.identifier(element2);
    stk::mesh::EntityId face_global_id = 0;

    if(element1_global_id < element2_global_id)
    {
        face_global_id = get_element_face_multiplier() * element1_global_id + element1_side_id;
    }
    else
    {
        SideId side_id = elementGraph.get_side_from_element1_to_locally_owned_element2(element2, element1);
        EXPECT_TRUE(side_id != -1);
        face_global_id = get_element_face_multiplier() * element2_global_id + side_id;
    }

    return face_global_id;
}

stk::mesh::EntityId get_face_id_for_remotely_connected_element(stk::mesh::EntityId local_element_id, SideId local_side_id, stk::mesh::EntityId remote_element_id, SideId remote_side_id)
{
    stk::mesh::EntityId face_global_id = 0;
    if(local_element_id < remote_element_id)
    {
        face_global_id = get_element_face_multiplier() * local_element_id + local_side_id;
    }
    else
    {
        face_global_id = get_element_face_multiplier() * remote_element_id + remote_side_id;
    }

    return face_global_id;
}

stk::mesh::Permutation get_permutation_for_new_face(const parallel_info& parallel_edge_info, stk::mesh::EntityId local_element_id, stk::mesh::EntityId remote_element_id)
{
    stk::mesh::Permutation perm;
    if(local_element_id < remote_element_id)
    {
        perm = static_cast<stk::mesh::Permutation>(0);
    }
    else
    {
        perm = static_cast<stk::mesh::Permutation>(parallel_edge_info.m_permutation);
    }
    return perm;
}

void create_or_delete_shared_face(stk::mesh::BulkData& bulkData, const parallel_info& parallel_edge_info, const ElemElemGraph& elementGraph,
        stk::mesh::EntityId local_id, stk::mesh::EntityId remote_id, bool create_face, stk::mesh::Part& faces_part,
        std::vector<sharing_info> &shared_modified, stk::mesh::EntityVector &deletedEntities)
{
    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, local_id);
    SideId side_id = elementGraph.get_side_from_element1_to_remote_element2(element, remote_id);
    ASSERT_TRUE(side_id != -1);

    ProcId this_proc_id = bulkData.parallel_rank();
    ProcId other_proc = parallel_edge_info.m_other_proc;
    ProcId owning_proc = this_proc_id < other_proc ? this_proc_id : other_proc;

    SideId remote_side_id = parallel_edge_info.m_other_side_ord;
    stk::mesh::EntityId face_global_id = get_face_id_for_remotely_connected_element(local_id, side_id, remote_id, remote_side_id);

    stk::mesh::ConnectivityOrdinal side_ord = static_cast<stk::mesh::ConnectivityOrdinal>(side_id);
    std::string msg = "Program error. Contact sierra-help@sandia.gov for support.";

    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
    if(create_face)
    {
        stk::mesh::Permutation perm = get_permutation_for_new_face(parallel_edge_info, local_id, remote_id);

        ThrowRequireMsg(!is_id_already_in_use_locally(bulkData, side_rank, face_global_id), msg);
        ThrowRequireMsg(!does_side_exist_with_different_permutation(bulkData, element, side_ord, perm), msg);
        ThrowRequireMsg(!does_element_side_exist(bulkData, element, side_ord), msg);

        stk::mesh::Entity face = connect_face_to_element(bulkData, element, face_global_id, side_ord, perm, faces_part);
        shared_modified.push_back(sharing_info(face, other_proc, owning_proc));
    }
    else
    {
        ThrowRequireMsg(is_id_already_in_use_locally(bulkData, side_rank, face_global_id), msg);
        ThrowRequireMsg(does_element_side_exist(bulkData, element, side_ord), msg);

        stk::mesh::Entity face = bulkData.get_entity(stk::topology::FACE_RANK, face_global_id);
        ThrowRequireMsg(bulkData.is_valid(face), msg);
        deletedEntities.push_back(face);
    }
}

TEST(ElementGraph, test_element_death)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) <= 2)
    {
        std::string dimension = unitTestUtils::getOption("--zdim", "none");

        int xdim = 3;
        if(dimension != "none")
        {
            xdim = std::atoi(dimension.c_str());
        }

        int ydim = xdim;
        int zdim = xdim; //  * stk::parallel_machine_size(comm);

        std::string filename = get_name_of_generated_mesh(xdim, ydim, zdim);

        {
            unsigned spatialDim = 3;
            stk::mesh::MetaData meta(spatialDim);
            stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
            stk::io::put_io_part_attribute(faces_part);
            BulkDataElementGraphTester bulkData(meta, comm);

            stk::mesh::Part& inactive = meta.declare_part("inactive", stk::topology::ELEMENT_RANK);
            stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);
            stk::unit_test_util::write_mesh_using_stk_io("orig.exo", bulkData, bulkData.parallel());

            double start_graph = stk::wall_time();

            ASSERT_TRUE(meta.get_part("block_1") != NULL);

            stk::mesh::Part& block_1 = *meta.get_part("block_1");

            ElemElemGraph elementGraph(bulkData);

            double elapsed_graph_time = stk::wall_time() - start_graph;
            std::cerr << "Time to create graph: " << elapsed_graph_time << std::endl;

            stk::mesh::EntityRank side_rank = meta.side_rank();

            double total_mod_ed = 0;
            double total_mod_ec = 0;

            {
                double wall_time_start = stk::wall_time();

                int num_time_steps = xdim * ydim * zdim;

                for(int i = 0; i < num_time_steps; ++i)
                {
                    stk::mesh::EntityVector killedElements = get_killed_elements(bulkData, i, inactive);
                    move_killled_elements_to_part(bulkData, killedElements, block_1, inactive);

                    std::vector<sharing_info> shared_modified;
                    stk::mesh::EntityVector deletedEntities;

                    std::vector<graphEdgeProc> elements_to_comm = get_elements_to_communicate(bulkData, killedElements, elementGraph);
                    std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> > remote_edges;

                    communicate_killed_entities(bulkData, elements_to_comm, remote_edges);

                    bulkData.modification_begin();

                    for(size_t re=0;re<remote_edges.size();++re)
                    {
                        stk::mesh::EntityId local_id = remote_edges[re].first;
                        stk::mesh::EntityId remote_id = remote_edges[re].second;

                        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, local_id);
                        bool create_face = true;
                        if(bulkData.bucket(element).member(inactive))
                        {
                            create_face = false;
                        }

                        parallel_info &parallel_edge_info = elementGraph.get_parallel_edge_info(element, remote_id);
                        create_or_delete_shared_face(bulkData, parallel_edge_info, elementGraph, local_id, remote_id, create_face, faces_part,
                                shared_modified, deletedEntities);
                        parallel_edge_info.m_in_part = false;
                    }

                    std::vector<ElementSidePair> element_side_pairs;
                    element_side_pairs.reserve(get_element_face_multiplier() * killedElements.size());

                    for(size_t k = 0; k < killedElements.size(); ++k)
                    {
                        stk::mesh::Entity this_elem_entity = killedElements[k];

                        for(size_t j = 0; j < elementGraph.get_num_connected_elems(this_elem_entity); ++j)
                        {
                            if(elementGraph.is_connected_elem_locally_owned(this_elem_entity, j))
                            {
                                stk::mesh::Entity other_element = elementGraph.get_connected_element(this_elem_entity, j);
                                SideId side_id = elementGraph.get_side_id_to_connected_element(this_elem_entity, j);

                                bool is_other_element_alive = !bulkData.bucket(other_element).member(inactive);
                                if(is_other_element_alive)
                                {
                                    stk::mesh::EntityId face_global_id = get_face_global_id(bulkData, elementGraph, this_elem_entity, other_element, side_id);

                                    stk::mesh::Entity face = stk::mesh::impl::get_or_create_face_at_element_side(bulkData, this_elem_entity, side_id,
                                            face_global_id, faces_part);

                                    const stk::mesh::Entity* side_nodes = bulkData.begin_nodes(face);
                                    unsigned num_side_nodes = bulkData.num_nodes(face);
                                    stk::mesh::EntityVector side_nodes_vec(side_nodes, side_nodes + num_side_nodes);


                                    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm = stk::mesh::get_ordinal_and_permutation(bulkData, other_element,
                                            stk::topology::FACE_RANK,
                                            side_nodes_vec);
                                    bulkData.declare_relation(other_element, face, ord_and_perm.first, ord_and_perm.second);
                                }
                                else
                                {
                                    unsigned num_faces = bulkData.num_faces(this_elem_entity);
                                    const stk::mesh::Entity *faces = bulkData.begin_faces(this_elem_entity);
                                    const stk::mesh::ConnectivityOrdinal *ordinals = bulkData.begin_face_ordinals(this_elem_entity);

                                    stk::mesh::Entity face;

                                    for (unsigned ii=0;ii<num_faces;++ii)
                                    {
                                        if (ordinals[ii] == static_cast<stk::mesh::ConnectivityOrdinal>(side_id))
                                        {
                                            face = faces[ii];
                                            break;
                                        }
                                    }

                                    EXPECT_TRUE(bulkData.is_valid(face));
                                    deletedEntities.push_back(face);
                                }
                            }
                            else
                            {
                                stk::mesh::EntityId local_id = bulkData.identifier(this_elem_entity);
                                stk::mesh::EntityId remote_id = elementGraph.get_entity_id_of_remote_element(this_elem_entity, j);

                                parallel_info &parallel_edge_info = elementGraph.get_parallel_edge_info(this_elem_entity, remote_id);
                                bool other_element_active = parallel_edge_info.m_in_part;
                                bool create_face = false;
                                if(other_element_active)
                                {
                                    create_face = true;
                                }

                                create_or_delete_shared_face(bulkData, parallel_edge_info, elementGraph, local_id, remote_id, create_face, faces_part,
                                        shared_modified, deletedEntities);
                            }
                        }
                    }

                    double start_mod1 = stk::wall_time();
                    bulkData.my_modification_end_for_entity_creation(shared_modified);
                    double elapsed_mod1 = stk::wall_time() - start_mod1;
                    total_mod_ec += elapsed_mod1;

                    double start_mod2 = stk::wall_time();
                    bulkData.modification_begin();
                    stk::mesh::impl::delete_entities_and_upward_relations(bulkData, deletedEntities);
                    bulkData.modification_end();
                    double elapsed_mod2 = stk::wall_time() - start_mod2;
                    total_mod_ed += elapsed_mod2;
                }

                double elapsed_time = stk::wall_time() - wall_time_start;

                stk::mesh::Selector sel = block_1;
                std::vector<size_t> counts1;
                stk::mesh::comm_mesh_counts(bulkData, counts1, &sel);

                size_t num_active = counts1[stk::topology::ELEM_RANK];

                stk::mesh::Selector sel2 = faces_part;
                stk::mesh::comm_mesh_counts(bulkData, counts1, &sel2);

                size_t num_faces = counts1[side_rank];

                EXPECT_EQ(2u, num_active);
                EXPECT_EQ(5u, num_faces);

                if(stk::parallel_machine_rank(comm) == 0)
                {
                    std::cerr << "Total time: " << elapsed_time << std::endl;
                    std::cerr << "\tTotal entity creation mod end time: " << total_mod_ec << std::endl;
                    std::cerr << "\tTotal entity deletion mod end time: " << total_mod_ed << std::endl;
                    std::cerr << "Total # of alive elements: " << num_active << std::endl;
                    std::cerr << "Total # of faces: " << num_faces << std::endl;

//                    {
//                        const stk::mesh::BucketVector &face_buckets = bulkData.buckets(stk::topology::FACE_RANK);
//                        for(size_t i=0;i<face_buckets.size();++i)
//                        {
//                            const stk::mesh::Bucket &bucket = *face_buckets[i];
//                            for(size_t j=0;j<bucket.size();++j)
//                            {
//                                std::cerr << "Face " << bulkData.identifier(bucket[j]) << " exists.\n";
//                            }
//                        }
//                    }
//
//                    {
//                        const stk::mesh::BucketVector &buckets = bulkData.buckets(stk::topology::ELEM_RANK);
//                        for(size_t i=0;i<buckets.size();++i)
//                        {
//                            const stk::mesh::Bucket &bucket = *buckets[i];
//                            if(!bucket.member(inactive))
//                            {
//                                for(size_t j=0;j<bucket.size();++j)
//                                {
//                                    std::cerr << "Element " << bulkData.identifier(bucket[j]) << " exists.\n";
//                                }
//                            }
//                        }
//                    }
                }

                stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());
            }
        }
    }
}

SideId get_side_from_element1_to_element2(const ElementGraph& elem_graph, const SidesForElementGraph &via_side, LocalId element_id1, LocalId element_id2)
{
    SideId side = -1;
    const std::vector<LocalId>& conn_elements = elem_graph[element_id1];

    std::vector<LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), element_id2);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = via_side[element_id1][index];
    }
    return side;
}

int check_connectivity(const ElementGraph& elem_graph, const SidesForElementGraph &via_side,
        LocalId element_id1, LocalId element_id2)
{
    int side=-1;
    if(element_id1 >=0 && element_id1 <=2 && element_id2 >=0 && element_id2 <=2)
    {
        side = get_side_from_element1_to_element2(elem_graph, via_side, element_id1, element_id2);
    }
    return side;
}

//BeginDocExample1
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
//EndDocExample1

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

void fill_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph, SidesForElementGraph& via_sides)
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

stk::mesh::EntityVector get_elements_to_communicate(const stk::mesh::BulkData& bulkData)
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

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements_to_communicate)
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

//BeginDocExample2
void add_possibly_connected_elements_to_graph_using_side_nodes(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
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
//EndDocExample2


void fill_parallel_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info)
{
    stk::mesh::EntityVector elements_to_communicate = get_elements_to_communicate(bulkData);

    stk::CommSparse comm(bulkData.parallel());

    pack_shared_side_nodes_of_elements(comm, bulkData, elements_to_communicate);
    comm.allocate_buffers();

    pack_shared_side_nodes_of_elements(comm, bulkData, elements_to_communicate);
    comm.communicate();

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
