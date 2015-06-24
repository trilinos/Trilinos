#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <stk_mesh/base/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "UnitTestElementDeathUtils.hpp"

namespace
{

class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
    ElemElemGraphTester(stk::mesh::BulkData& bulkData)
      : ElemElemGraph(bulkData) {};

    virtual ~ElemElemGraphTester() {};

    void fill_graph() { ElemElemGraph::fill_graph(); }

    void fill_parallel_graph(stk::mesh::impl::ElemSideToProcAndFaceId& elem_side_comm) { ElemElemGraph::fill_parallel_graph(elem_side_comm); }

    stk::mesh::impl::ElementGraph & get_element_graph() { return m_elem_graph; }
    stk::mesh::impl::SidesForElementGraph & get_via_sides() { return m_via_sides; }
    stk::mesh::impl::ParallelGraphInfo & get_parallel_graph_info() { return m_parallel_graph_info; }
    size_t get_graph_size() { return m_local_id_to_element_entity.size(); }

    int check_local_connectivity(stk::mesh::Entity elem1, stk::mesh::Entity elem2)
    {
        int side=-1;
        if (is_valid_graph_element(elem1) && is_valid_graph_element(elem2)) {
            side = get_side_from_element1_to_locally_owned_element2(elem1, elem2);
        }
        return side;
    }

    int check_remote_connectivity(stk::mesh::Entity elem, stk::mesh::EntityId other_elem_id)
    {
        int side=-1;
        if (is_valid_graph_element(elem)) {
            side = get_side_from_element1_to_remote_element2(elem, other_elem_id);
        }
        return side;
    }
};

class BulkDataElementGraphTester : public stk::mesh::BulkData
{

public:

    BulkDataElementGraphTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }

    ~BulkDataElementGraphTester(){}

    bool my_internal_modification_end_for_skin_mesh(stk::mesh::EntityRank entity_rank, stk::mesh::impl::MeshModification::modification_optimization opt, stk::mesh::Selector selectedToSkin,
            const stk::mesh::Selector * only_consider_second_element_from_this_selector = 0)
    {
        return this->internal_modification_end_for_skin_mesh(entity_rank, opt, selectedToSkin, only_consider_second_element_from_this_selector);
    }

    bool my_modification_end_for_entity_creation(const std::vector<stk::mesh::sharing_info>& shared_modified, stk::mesh::impl::MeshModification::modification_optimization opt = stk::mesh::impl::MeshModification::MOD_END_SORT)
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

using namespace stk::mesh::impl;
using namespace stk::mesh;

// Not to be used with ElemElemGraph or ElemElemGraphTester class.
bool is_valid_graph_element(const impl::ElementGraph &elem_graph, stk::mesh::impl::LocalId elem_id);

// Not to be used with ElemElemGraph or ElemElemGraphTester class.
int check_connectivity(const impl::ElementGraph &elem_graph, const impl::SidesForElementGraph &via_sides,
                       stk::mesh::impl::LocalId element_id1, stk::mesh::impl::LocalId element_id2);

// Not to be used with ElemElemGraph or ElemElemGraphTester class.
int get_side_from_element1_to_element2(const impl::ElementGraph &elem_graph,
                                       const impl::SidesForElementGraph &via_sides,
                                       stk::mesh::impl::LocalId element1_local_id,
                                       stk::mesh::impl::LocalId other_element_id);

void add_element_side_pairs_for_unused_sides(LocalId elementId, stk::topology topology, const std::vector<int> &internal_sides,
        std::vector<ElementSidePair>& element_side_pairs);

void create_faces_using_graph(BulkDataElementGraphTester& bulkData, stk::mesh::Part& part);


TEST(ElementGraph, add_elements_to_graph_serial)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    std::vector<size_t> mem_usage;

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);

        ElemElemGraphTester elem_graph(bulkData);

        EXPECT_EQ(0u, elem_graph.size());

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(4, numLocallyOwnedElems);

        stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
        stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
        stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);

        stk::mesh::EntityVector elements_to_add;
        elements_to_add.push_back(elem1);
        elements_to_add.push_back(elem2);
        elements_to_add.push_back(elem3);
        elements_to_add.push_back(elem4);

        for (unsigned i=0; i<elements_to_add.size(); i++)
        {
        	EXPECT_TRUE(bulkData.is_valid(elements_to_add[i]));
        	EXPECT_EQ(0, bulkData.parallel_owner_rank(elements_to_add[i]));
        }

        elem_graph.add_elements_to_graph(elements_to_add);

        EXPECT_EQ(4u, elem_graph.size());
        EXPECT_EQ(6u, elem_graph.num_edges());

        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem2, elem3));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem3, elem2));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));
    }
}

TEST(ElementGraph, DISABLED_add_elements_to_graph_parallel)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int p_rank = stk::parallel_machine_rank(comm);

    std::vector<size_t> mem_usage;

    if(stk::parallel_machine_size(comm) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);

        ElemElemGraphTester elem_graph(bulkData);

        EXPECT_EQ(0u, elem_graph.size());

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(2, numLocallyOwnedElems);

        stk::mesh::Entity elem1 = stk::mesh::Entity();
        stk::mesh::Entity elem2 = stk::mesh::Entity();
        stk::mesh::Entity elem3 = stk::mesh::Entity();
        stk::mesh::Entity elem4 = stk::mesh::Entity();

        stk::mesh::EntityVector elements_to_add;

        if (p_rank == 0)
        {
            elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);

            elements_to_add.push_back(elem1);
            elements_to_add.push_back(elem2);
        }
        else
        {
            elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);

            elements_to_add.push_back(elem3);
            elements_to_add.push_back(elem4);
        }

        if (p_rank == 0)
        {
            EXPECT_TRUE(bulkData.is_valid(elem1));
            EXPECT_EQ(0, bulkData.parallel_owner_rank(elem1));
            EXPECT_TRUE(bulkData.is_valid(elem2));
            EXPECT_EQ(0, bulkData.parallel_owner_rank(elem2));
        }
        else
        {
            EXPECT_TRUE(bulkData.is_valid(elem3));
            EXPECT_EQ(1, bulkData.parallel_owner_rank(elem3));
            EXPECT_TRUE(bulkData.is_valid(elem4));
            EXPECT_EQ(1, bulkData.parallel_owner_rank(elem4));
        }

        elem_graph.add_elements_to_graph(elements_to_add);

        EXPECT_EQ(2u, elem_graph.size());
        EXPECT_EQ(3u, elem_graph.num_edges());

        if (p_rank ==0)
        {
            const stk::mesh::EntityId elem3_global_id = 3;
            EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
            EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
            EXPECT_EQ(5, elem_graph.check_remote_connectivity(elem2, elem3_global_id));
            EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        }
        else
        {
            const stk::mesh::EntityId elem2_global_id = 2;
            EXPECT_EQ(4, elem_graph.check_remote_connectivity(elem3, elem2_global_id));
            EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
            EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
            EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));
        }
    }
}

TEST(ElementGraph, delete_elements_from_graph_serial)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(4, numLocallyOwnedElems);

        stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
        stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
        stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);

        ElemElemGraphTester elem_graph(bulkData);

        EXPECT_EQ(4u, elem_graph.size());
        EXPECT_EQ(6u, elem_graph.num_edges());

        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem2, elem3));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem3, elem2));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));

        std::vector <stk::mesh::Entity> elems_to_delete;
        elems_to_delete.push_back(elem2);
        elems_to_delete.push_back(elem3);

        for (unsigned i=0; i<elems_to_delete.size(); i++)
        {
        	EXPECT_TRUE(bulkData.is_valid(elems_to_delete[i]));
        	EXPECT_EQ(0, bulkData.parallel_owner_rank(elems_to_delete[i]));
        }

        elem_graph.delete_elements_from_graph(elems_to_delete);

        EXPECT_EQ(2u, elem_graph.size());
        EXPECT_EQ(0u, elem_graph.num_edges());

        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem2));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem1));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem3, elem2));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem3, elem4));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem4, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));
    }
}

TEST(ElementGraph, add_and_delete_elements_from_graph_serial)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) == 1)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulkData(meta, comm);

        ElemElemGraphTester elem_graph(bulkData);

        EXPECT_EQ(0u, elem_graph.size());
        EXPECT_EQ(0u, elem_graph.num_edges());

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(4, numLocallyOwnedElems);

        stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
        stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
        stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);

        stk::mesh::EntityVector elements_to_add;
        elements_to_add.push_back(elem1);
        elements_to_add.push_back(elem2);
        elements_to_add.push_back(elem3);
        elements_to_add.push_back(elem4);

        for (unsigned i=0; i<elements_to_add.size(); i++)
        {
        	EXPECT_TRUE(bulkData.is_valid(elements_to_add[i]));
        	EXPECT_EQ(0, bulkData.parallel_owner_rank(elements_to_add[i]));
        }

        elem_graph.add_elements_to_graph(elements_to_add);

        EXPECT_EQ(4u, elem_graph.size());
        EXPECT_EQ(6u, elem_graph.num_edges());

        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem2, elem3));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem3, elem2));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));

        std::vector <stk::mesh::Entity> elems_to_delete;
        elems_to_delete.push_back(elem2);
        elems_to_delete.push_back(elem3);

        elem_graph.delete_elements_from_graph(elems_to_delete);

        EXPECT_EQ(2u, elem_graph.size());
        EXPECT_EQ(0u, elem_graph.num_edges());

        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem2));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem1));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem3, elem2));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem3, elem4));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem4, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));

        elements_to_add.clear();
        elements_to_add.push_back(elem2);
        elements_to_add.push_back(elem3);

        elem_graph.add_elements_to_graph(elements_to_add);

        EXPECT_EQ(4u, elem_graph.size());
        EXPECT_EQ(6u, elem_graph.num_edges());

        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem2, elem3));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem3, elem2));
        EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
        EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
        EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));
    }
}

TEST(ElementGraph, HexAddShellSerial)
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |   |
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Added single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size != 1) {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2 };

    // Build the base hex mesh
    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    mesh.modification_end();

    // Initialize the graph based on the existing hex mesh
    ElemElemGraphTester elem_graph(mesh);
    EXPECT_EQ(1u, elem_graph.size());

    // Add a shell
    mesh.modification_begin();
    stk::mesh::EntityVector added_shells;
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
    }
    mesh.modification_end();

    elem_graph.add_elements_to_graph(added_shells);

    EXPECT_EQ(2u, elem_graph.size());

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elem_graph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elem_graph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(shell2, elem_graph.get_connected_element(hex1, 0));
    EXPECT_TRUE(elem_graph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,             elem_graph.get_num_connected_elems(shell2));
    EXPECT_EQ(1,              elem_graph.get_side_id_to_connected_element(shell2, 0));
    EXPECT_EQ(hex1, elem_graph.get_connected_element(shell2, 0));
    EXPECT_TRUE(elem_graph.is_connected_elem_locally_owned(shell2, 0));

}

TEST( ElementGraph, DISABLED_HexAddShellAddShellSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |3.0|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2, 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    mesh.modification_begin();
    stk::mesh::EntityVector added_shells;
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
    }
    mesh.modification_end();

    elemElemGraph.add_elements_to_graph(added_shells);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
    EXPECT_EQ(shell2, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell2, 0));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
}

TEST( ElementGraph, HexAddShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Added single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    mesh.modification_begin();
    stk::mesh::EntityVector added_shells;
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
    }
    mesh.modification_end();

    elemElemGraph.add_elements_to_graph(added_shells);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Hex Element 2
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell3, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
}

TEST( ElementGraph, DISABLED_HexAddShellAddShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Added two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    mesh.modification_begin();
    stk::mesh::EntityVector added_shells;
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
    }
    mesh.modification_end();

    elemElemGraph.add_elements_to_graph(added_shells);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex1, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 1));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex2, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell3, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell4, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell4, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell4, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
}


void add_element_side_pairs_for_unused_sides(LocalId elementId, stk::topology topology, const std::vector<int> &internal_sides,
        std::vector<ElementSidePair>& element_side_pairs)
{
    size_t num_sides = topology.num_sides();
    std::vector<int> elem_sides;

    if (internal_sides.size() < num_sides)
    {
        elem_sides.assign(num_sides, -1);
        for(size_t j=0; j<internal_sides.size(); ++j)
        {
            int sideId = internal_sides[j];
            elem_sides[sideId] = internal_sides[j];
        }

        for(size_t j=0; j<num_sides; ++j)
        {
            if (elem_sides[j] == -1)
            {
                int sideId = j;
                element_side_pairs.push_back(std::make_pair(elementId, sideId));
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

void change_local_id_to_negative_global_id(ElementGraph &elem_graph, LocalId elem_local_id, stk::mesh::EntityId elem_global_id)
{
    ThrowRequire(is_valid_graph_element(elem_graph, elem_local_id));
    for(unsigned id = 0; id < elem_graph.size(); ++id)
    {
        std::vector <LocalId>::iterator pos_of_move_elem_in_current_id = std::find(elem_graph[id].begin(), elem_graph[id].end(), elem_local_id);
        if (pos_of_move_elem_in_current_id != elem_graph[id].end())
        {
            int index_of_move_elem_in_current_id = pos_of_move_elem_in_current_id - elem_graph[id].begin();
            elem_graph[id][index_of_move_elem_in_current_id] = -elem_global_id;
        }
    }
}

void change_negative_global_id_to_local_id(ElementGraph &elem_graph, stk::mesh::EntityId elem_global_id, LocalId elem_local_id)
{
    ThrowRequire(is_valid_graph_element(elem_graph, elem_local_id));
    for(unsigned id = 0; id < elem_graph.size(); ++id)
    {
        LocalId negative_elem_global_id = -elem_global_id;
        std::vector <LocalId>::iterator pos_of_negative_global_id = std::find(elem_graph[id].begin(), elem_graph[id].end(), negative_elem_global_id);
        if (pos_of_negative_global_id != elem_graph[id].end())
        {
            int index_of_move_elem_in_current_id = pos_of_negative_global_id - elem_graph[id].begin();
            elem_graph[id][index_of_move_elem_in_current_id] = elem_local_id;
        }
    }
}

void test_parallel_graph_info(const ElementGraph& elem_graph, const ParallelGraphInfo& parallel_graph_info,
        LocalId this_element, LocalId other_element, int other_proc, int other_side_ord, int permutation)
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


void create_faces_using_graph(BulkDataElementGraphTester& bulkData, stk::mesh::Part& part)
{
    double wall_time_start = stk::wall_time();

    std::vector<unsigned> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numElems = counts[stk::topology::ELEM_RANK];

    stk::mesh::EntityVector local_id_to_element_entity(numElems, Entity());
    std::vector<stk::topology> element_topologies(numElems);
    impl::set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

    ElemElemGraphTester elemElemGraph(bulkData);

    ElementGraph & elem_graph = elemElemGraph.get_element_graph();
    SidesForElementGraph & via_sides = elemElemGraph.get_via_sides();
    stk::mesh::impl::ParallelGraphInfo & parallel_graph_info = elemElemGraph.get_parallel_graph_info();

    double graph_time = stk::wall_time() - wall_time_start;
    wall_time_start = stk::wall_time();

    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
    stk::mesh::PartVector parts {&part};
    //BeginDocExample4

    bulkData.modification_begin();

    std::vector<stk::mesh::sharing_info> shared_modified;

    for(size_t i = 0; i < elem_graph.size(); ++i)
    {
        const std::vector<LocalId>& connected_elements = elem_graph[i];
        stk::mesh::Entity element1 = local_id_to_element_entity[i];

        LocalId this_element = i;

        for(size_t j = 0; j < connected_elements.size(); ++j)
        {
            if(this_element < connected_elements[j] && connected_elements[j] >= 0)
            {
                stk::mesh::EntityId face_global_id = impl::get_element_side_multiplier() * bulkData.identifier(element1) + via_sides[i][j];
                if ( impl::is_id_already_in_use_locally(bulkData, side_rank, face_global_id) )
                {

                }
                stk::mesh::Entity face = stk::mesh::impl::get_or_create_face_at_element_side(bulkData, element1, via_sides[i][j],
                        face_global_id, stk::mesh::PartVector(1,&part));

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
                int other_proc = iter->second.m_other_proc;
                int other_side = iter->second.m_other_side_ord;

                int this_proc = bulkData.parallel_rank();
                int owning_proc = this_proc < other_proc ? this_proc : other_proc;

                stk::mesh::EntityId face_global_id = 0;
                stk::mesh::Permutation perm;
                if(owning_proc == this_proc)
                {
                    stk::mesh::EntityId id = bulkData.identifier(element1);
                    face_global_id = impl::get_element_side_multiplier() * id + via_sides[i][j];
                    perm = static_cast<stk::mesh::Permutation>(0);
                }
                else
                {
                    face_global_id = impl::get_element_side_multiplier() * other_element + other_side;
                    perm = static_cast<stk::mesh::Permutation>(iter->second.m_permutation);
                }

                stk::mesh::ConnectivityOrdinal side_ord = static_cast<stk::mesh::ConnectivityOrdinal>(via_sides[i][j]);

                std::string msg = "Program error. Contact sierra-help@sandia.gov for support.";

                ThrowRequireMsg(!impl::is_id_already_in_use_locally(bulkData, side_rank, face_global_id), msg);
                ThrowRequireMsg(!impl::does_side_exist_with_different_permutation(bulkData, element1, side_ord, perm), msg);
                ThrowRequireMsg(!impl::does_element_side_exist(bulkData, element1, side_ord), msg);

                stk::mesh::Entity face = impl::connect_side_to_element(bulkData, element1, face_global_id, side_ord, perm, parts);

                shared_modified.push_back(stk::mesh::sharing_info(face, other_proc, owning_proc));
            }
        }

        std::vector<ElementSidePair> element_side_pairs;
        add_element_side_pairs_for_unused_sides(i, element_topologies[i], via_sides[i], element_side_pairs);

        for(size_t j = 0; j < element_side_pairs.size(); j++)
        {
            stk::mesh::EntityId face_global_id = impl::get_element_side_multiplier() * bulkData.identifier(element1) + element_side_pairs[j].second;
            stk::mesh::impl::get_or_create_face_at_element_side(bulkData, element1, element_side_pairs[j].second,
                    face_global_id, stk::mesh::PartVector(1,&part));

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

//BeginDocExample1
std::vector<ElementSidePair>
skin_mesh(const SidesForElementGraph &via_side, const std::vector<stk::topology> &element_topologies)
{
    std::vector<ElementSidePair> element_side_pairs;

    size_t num_elems = via_side.size();
    for(size_t i=0; i<num_elems; ++i)
    {
        const std::vector<int>& internal_sides = via_side[i];
        add_element_side_pairs_for_unused_sides(i, element_topologies[i], internal_sides, element_side_pairs);
    }
    return element_side_pairs;
}

std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim)
{
    std::ostringstream os;
    os << "generated:" << xdim << "x" << ydim << "x" << zdim << "|nodeset:zZ";
    return os.str();
}

stk::mesh::EntityVector get_killed_elements(stk::mesh::BulkData& bulkData, const int killValue, const stk::mesh::Part& active)
{
    stk::mesh::EntityVector killedElements;
    const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);
    for(size_t b = 0; b < buckets.size(); ++b)
    {
        const stk::mesh::Bucket &bucket = *buckets[b];
        if(bucket.owned() && bucket.member(active))
        {
            for(size_t e = 0; e < bucket.size(); ++e)
            {
                stk::mesh::Entity entity = bucket[e];
                bool should_element_be_killed = bulkData.identifier(entity) < static_cast<stk::mesh::EntityId>(killValue);
                if(bulkData.bucket(entity).member(active) && should_element_be_killed == true)
                {
                    killedElements.push_back(bucket[e]);
                }
            }
        }
    }
    return killedElements;
}

void move_killled_elements_to_part(stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& killedElements, stk::mesh::Part& block_1, stk::mesh::Part& active)
{
    std::vector<stk::mesh::PartVector> add_parts(killedElements.size());
    std::vector<stk::mesh::PartVector> rm_parts(killedElements.size());

    stk::mesh::PartVector add_part_vec;
    stk::mesh::PartVector rm_part_vec = { &block_1, &active };

    for (size_t j=0;j<killedElements.size();++j)
    {
        add_parts[j] = add_part_vec;
        rm_parts[j] = rm_part_vec;
    }

    bulkData.batch_change_entity_parts(killedElements, add_parts, rm_parts);
}

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

        ElemElemGraphTester elemElemGraph(bulkData);

        size_t expectedNumElems = counts[stk::topology::ELEM_RANK];
        ASSERT_EQ(expectedNumElems, elemElemGraph.get_graph_size());

        ElementGraph & elem_graph = elemElemGraph.get_element_graph();
        SidesForElementGraph & via_sides = elemElemGraph.get_via_sides();

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after ElemElemGraphTester constructor");
        mem_usage.push_back(stk::get_memory_usage_now());

        int left_side_id = 4;
        int right_side_id = 5;

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

        ElemElemGraphTester elemElemGraph(bulkData);

        size_t expectedNumElems = counts[stk::topology::ELEM_RANK];
        ASSERT_EQ(expectedNumElems, elemElemGraph.get_graph_size());

        ElementGraph & elem_graph = elemElemGraph.get_element_graph();
        SidesForElementGraph & via_sides = elemElemGraph.get_via_sides();

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after ElemElemGraphTester constructor");
        mem_usage.push_back(stk::get_memory_usage_now());

        LocalId element_to_test_local_id = std::numeric_limits<LocalId>::max();
        int side_id = -1;
        int left_side_id = 4;
        int right_side_id = 5;

        EXPECT_EQ(2u, elem_graph.size());

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
            if (static_cast<LocalId>(i) == element_to_test_local_id)
            {
                // Element on parallel boundary
                ASSERT_EQ(2u, elem_graph[i].size());
                ASSERT_EQ(2u, via_sides[i].size());
                ASSERT_GE(-1, elem_graph[i][1]);
                ASSERT_EQ(side_id, via_sides[i][1]);
            }
            else
            {
                EXPECT_EQ(1u, elem_graph[i].size());
                EXPECT_EQ(1u, via_sides[i].size());
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

        ElemElemGraphTester elemElemGraph(bulkData);

        size_t expectedNumElems = counts[stk::topology::ELEM_RANK];

        ElementGraph & elem_graph = elemElemGraph.get_element_graph();
        SidesForElementGraph & via_sides = elemElemGraph.get_via_sides();

        if ( check_results )
        {
            ASSERT_EQ(expectedNumElems, elemElemGraph.get_graph_size());
        }

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after ElemElemGraphTester constructor");
        mem_usage.push_back(stk::get_memory_usage_now());

        stk::mesh::EntityVector local_id_to_element_entity(numElems, Entity());
        std::vector<stk::topology> element_topologies(numElems);
        impl::set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        std::vector<ElementSidePair> elem_side_pairs = skin_mesh(via_sides, element_topologies);

        bulkData.modification_begin();

        for(size_t face_index=0; face_index<elem_side_pairs.size(); ++face_index)
        {
            stk::mesh::Entity element = local_id_to_element_entity[elem_side_pairs[face_index].first];
            stk::mesh::EntityId face_global_id = face_index + 1;
            stk::mesh::declare_element_side(bulkData, face_global_id, element, elem_side_pairs[face_index].second, &skin_part);
        }

        stk::mesh::Selector element_selector = bulkData.mesh_meta_data().locally_owned_part();
        bulkData.my_internal_modification_end_for_skin_mesh(stk::topology::FACE_RANK, stk::mesh::impl::MeshModification::MOD_END_SORT, element_selector, NULL);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after create-faces");
        mem_usage.push_back(stk::get_memory_usage_now());

        stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());

        if ( check_results )
        {
            int left_side_id = 4;
            int right_side_id = 5;

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

void change_entity_owner(stk::mesh::BulkData &bulkData, stk::mesh::ElemElemGraph &elem_graph, std::vector< std::pair< stk::mesh::Entity, int > > &elem_proc_pairs_to_move, stk::mesh::Part *active_part=NULL)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
    impl::ParallelGraphInfo parallel_graph;

    const std::vector<stk::mesh::EntityId> &suggested_face_id_vector = elem_graph.get_suggested_side_ids();
    size_t num_suggested_face_ids_used = 0;

    for (size_t i=0; i<elem_proc_pairs_to_move.size(); i++)
    {
        stk::mesh::Entity elem_to_send = elem_proc_pairs_to_move[i].first;
        int destination_proc = elem_proc_pairs_to_move[i].second;
        stk::mesh::EntityId elem_global_id = bulkData.identifier(elem_to_send);

        size_t num_connected_elements = elem_graph.get_num_connected_elems(elem_to_send);
        stk::topology elem_topology = bulkData.bucket(elem_to_send).topology();
        const stk::mesh::Entity *elem_nodes = bulkData.begin_nodes(elem_to_send);
        for (size_t k=0; k<num_connected_elements; k++)
        {
            if (elem_graph.is_connected_elem_locally_owned(elem_to_send, k))
            {
                int side_id = elem_graph.get_side_id_to_connected_element(elem_to_send, k);
                stk::mesh::Entity connected_element = elem_graph.get_connected_element(elem_to_send, k);
                impl::LocalId local_id = elem_graph.get_local_element_id(connected_element);
                stk::topology side_topology = elem_topology.side_topology(side_id);
                std::vector<stk::mesh::Entity> side_nodes(side_topology.num_nodes());

                elem_topology.side_nodes(elem_nodes, side_id, side_nodes.begin());
                stk::mesh::OrdinalAndPermutation ordperm = get_ordinal_and_permutation(bulkData, elem_to_send, side_rank, side_nodes);

                std::pair<impl::LocalId, stk::mesh::EntityId> key(local_id, elem_global_id);
                stk::mesh::EntityId face_id = suggested_face_id_vector[num_suggested_face_ids_used];
                num_suggested_face_ids_used++;
                impl::parallel_info p_info(destination_proc, side_id, ordperm.second, face_id);
                if(active_part != NULL)
                {
                    p_info.m_in_part = bulkData.bucket(connected_element).member(*active_part);
                }
                parallel_graph.insert(std::make_pair(key, p_info));
            }
        }
    }

    elem_graph.set_num_side_ids_used(num_suggested_face_ids_used);

    bulkData.change_entity_owner(elem_proc_pairs_to_move);

    elem_graph.change_entity_owner(elem_proc_pairs_to_move, parallel_graph);
}

void change_entity_owner_hex_test_2_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);
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
        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(2, numLocallyOwnedElems);

        stk::mesh::ElemElemGraph elem_graph(bulkData);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        // Create a vector of the elements to be moved
        std::vector <stk::mesh::Entity> elems_to_move;

        stk::mesh::EntityId elem_2_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

        if (proc == 0)
        {
            int side_from_elem2_to_elem3 = elem_graph.get_side_from_element1_to_remote_element2(elem_2, stk::mesh::EntityId(3));
            int side_from_elem2_to_elem1 = elem_graph.get_side_from_element1_to_locally_owned_element2(elem_2, bulkData.get_entity(stk::topology::ELEM_RANK,1));

            EXPECT_EQ(5, side_from_elem2_to_elem3);
            EXPECT_EQ(4, side_from_elem2_to_elem1);

            elems_to_move.push_back(elem_2);

            int other_proc = 1;
            for (unsigned i=0; i<elems_to_move.size(); i++)
            {
                EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
                EXPECT_EQ(0, bulkData.parallel_owner_rank(elems_to_move[i]));
            	elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
            }
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move);

        elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

        if (proc == 1)
        {
            EXPECT_TRUE(bulkData.is_valid(elem_2));
            EXPECT_EQ(1, bulkData.parallel_owner_rank(elem_2));

            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_2));

            stk::mesh::Entity elem = elem_graph.get_connected_element(elem_2, 1);
            ASSERT_TRUE(elem_graph.is_connected_elem_locally_owned(elem_2, 1));
            EXPECT_EQ(3u, bulkData.identifier(elem));

            stk::mesh::EntityId connected_elem_global_id = elem_graph.get_entity_id_of_remote_element(elem_2, 0);
            ASSERT_FALSE(elem_graph.is_connected_elem_locally_owned(elem_2, 0));
            EXPECT_EQ(1u, connected_elem_global_id);

            stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            int side_from_elem2_to_elem3 = elem_graph.get_side_from_element1_to_locally_owned_element2(elem_2, elem_3);
            EXPECT_EQ(5, side_from_elem2_to_elem3);

            int side_from_elem2_to_elem1 = elem_graph.get_side_from_element1_to_remote_element2(elem_2, stk::mesh::EntityId(1));
            EXPECT_EQ(4, side_from_elem2_to_elem1);

            impl::parallel_info &elem_2_to_1_p_info = elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1));
            int other_side_ord = elem_2_to_1_p_info.m_other_side_ord;
            EXPECT_EQ(5, other_side_ord);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_3, stk::mesh::EntityId(2)), std::logic_error);
        }
        if (proc == 0)
        {
            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            impl::parallel_info &p_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));
            int other_side_ord = p_info.m_other_side_ord;
            EXPECT_EQ(4, other_side_ord);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(3)), std::logic_error);
       }
    }
}

TEST(ElementGraph, test_change_entity_owner_2_procs_hex_mesh_with_aura)
{
    bool aura_on = true;
    change_entity_owner_hex_test_2_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_2_procs_hex_mesh_without_aura)
{
    bool aura_on = false;
    change_entity_owner_hex_test_2_procs(aura_on);
}

void change_entity_owner_hex_test_4_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);
    std::vector<double> wall_times;
    wall_times.reserve(10);
    std::vector<std::string> msgs;
    msgs.reserve(10);

    std::vector<size_t> mem_usage;

    wall_times.push_back(stk::wall_time());
    msgs.push_back("program-start");
    mem_usage.push_back(stk::get_memory_usage_now());

    if(stk::parallel_machine_size(comm) == 4)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after mesh-read");
        mem_usage.push_back(stk::get_memory_usage_now());

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
        EXPECT_EQ(1, numLocallyOwnedElems);

        stk::mesh::ElemElemGraph elem_graph(bulkData);

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after fill-graph");
        mem_usage.push_back(stk::get_memory_usage_now());

        // Create a vector of the elements to be moved
        std::vector <stk::mesh::Entity> elems_to_move;

        stk::mesh::EntityId elem_global_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_global_id);
        if (proc == 1)
        {
            elems_to_move.push_back(elem_to_move);

            int other_proc = 2;
            for (unsigned i=0; i<elems_to_move.size(); i++)
            {
                EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
                EXPECT_EQ(1, bulkData.parallel_owner_rank(elems_to_move[i]));
                elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
            }
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move);

        elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_global_id);

        if (proc == 2)
        {
            EXPECT_TRUE(bulkData.is_valid(elem_to_move));
            EXPECT_EQ(2, bulkData.parallel_owner_rank(elem_to_move));

            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_to_move));

            stk::mesh::Entity elem = elem_graph.get_connected_element(elem_to_move, 1);
            ASSERT_TRUE(elem_graph.is_connected_elem_locally_owned(elem_to_move, 1));
            EXPECT_EQ(3u, bulkData.identifier(elem));

            stk::mesh::EntityId connected_elem_global_id = elem_graph.get_entity_id_of_remote_element(elem_to_move, 0);
            ASSERT_FALSE(elem_graph.is_connected_elem_locally_owned(elem_to_move, 0));
            EXPECT_EQ(1u, connected_elem_global_id);

            stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_3, stk::mesh::EntityId(2)), std::logic_error);
        }
        else if (proc == 0)
        {
            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            impl::parallel_info &elem_1_to_2_p_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));
            EXPECT_EQ(2, elem_1_to_2_p_info.m_other_proc);
        }
    }
}

TEST(ElementGraph, test_change_entity_owner_4_procs_hex_mesh_with_aura)
{
    bool aura_on = true;
    change_entity_owner_hex_test_4_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_4_procs_hex_mesh_without_aura)
{
    bool aura_on = false;
    change_entity_owner_hex_test_4_procs(aura_on);
}

void setup_hex_shell_hex_mesh(stk::mesh::BulkData& bulkData)
{
//
//                proc 0               proc 1           proc 2
//
//               block_1          |   block_2  |      block_3
//
//          3---------------7        7            7-------------11
//          /|             /|       /|           /|             /|
//         / |            / |      / |          / |            / |
//        /  |           /  |     /  |         /  |           /  |
//       4--------------8   |    8   |        8--------------12  |
//       |   |          |   |    |   |        |   |          |   |
//       |   |   1      |   |    | 2 |        |   |   3      |   |
//       |   |          |   |    |   |        |   |          |   |
//       |   2----------|---6    |   6        |   6----------|---10
//       |  /           |  /     |  /         |  /           |  /
//       | /            | /      | /          | /            | /
//       |/             |/       |/           |/             |/
//       1--------------5        5            5--------------9

    stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    unsigned spatial_dimension = 3;
    meta.initialize(spatial_dimension, stk::mesh::entity_rank_names());

    stk::mesh::Field<double>& field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
    stk::mesh::put_field(field, meta.universal_part());

    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::HEX_8);
    meta.commit();

    bulkData.modification_begin();

    stk::mesh::EntityIdVector elem1_nodes {1, 2, 3, 4, 5, 6, 7, 8};
    stk::mesh::EntityIdVector elem2_nodes {5, 6, 7, 8};
    stk::mesh::EntityIdVector elem3_nodes {5, 6, 7, 8, 9, 10, 11, 12};

    stk::mesh::EntityId elemId = 1;
    if (bulkData.parallel_rank() == 0) {
        stk::mesh::declare_element(bulkData, block_1, elemId, elem1_nodes);
        stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
        stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
        stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
        bulkData.add_node_sharing(node5, 1);
        bulkData.add_node_sharing(node6, 1);
        bulkData.add_node_sharing(node7, 1);
        bulkData.add_node_sharing(node8, 1);
        bulkData.add_node_sharing(node5, 2);
        bulkData.add_node_sharing(node6, 2);
        bulkData.add_node_sharing(node7, 2);
        bulkData.add_node_sharing(node8, 2);
    }
    else if (bulkData.parallel_rank() == 1) {
        elemId = 2;
        stk::mesh::declare_element(bulkData, block_2, elemId, elem2_nodes);
        stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
        stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
        stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
        bulkData.add_node_sharing(node5, 0);
        bulkData.add_node_sharing(node6, 0);
        bulkData.add_node_sharing(node7, 0);
        bulkData.add_node_sharing(node8, 0);
        bulkData.add_node_sharing(node5, 2);
        bulkData.add_node_sharing(node6, 2);
        bulkData.add_node_sharing(node7, 2);
        bulkData.add_node_sharing(node8, 2);
    }
    else if (bulkData.parallel_rank() == 2) {
        elemId = 3;
        stk::mesh::declare_element(bulkData, block_3, elemId, elem3_nodes);
        stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
        stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
        stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
        bulkData.add_node_sharing(node5, 0);
        bulkData.add_node_sharing(node6, 0);
        bulkData.add_node_sharing(node7, 0);
        bulkData.add_node_sharing(node8, 0);
        bulkData.add_node_sharing(node5, 1);
        bulkData.add_node_sharing(node6, 1);
        bulkData.add_node_sharing(node7, 1);
        bulkData.add_node_sharing(node8, 1);
    }

    bulkData.modification_end();
}

void change_entity_owner_hex_shell_hex_test_3_procs(bool aura_on)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);
    if(stk::parallel_machine_size(comm) == 3)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
        if (!aura_on)
        {
            aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
        }
        stk::mesh::BulkData bulkData(meta, comm, aura_option);

        setup_hex_shell_hex_mesh(bulkData);

        stk::mesh::ElemElemGraph elem_graph(bulkData);

        const Entity hex1   = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        const Entity hex3   = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
        const Entity shell2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);

        if (proc == 0) {
            // Connectivity for Hex Element 1
            EXPECT_EQ(1u, elem_graph.get_num_connected_elems(hex1));
            EXPECT_EQ(5,  elem_graph.get_side_id_to_connected_element(hex1, 0));
            EXPECT_EQ(2u, elem_graph.get_entity_id_of_remote_element(hex1, 0));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(hex1, 0));
        }
        else if (proc == 1) {
            // Connectivity for Shell Element 3
            EXPECT_EQ(2u, elem_graph.get_num_connected_elems(shell2));
            EXPECT_EQ(0,  elem_graph.get_side_id_to_connected_element(shell2, 0));
            EXPECT_EQ(1,  elem_graph.get_side_id_to_connected_element(shell2, 1));
            EXPECT_EQ(3u, elem_graph.get_entity_id_of_remote_element(shell2, 0));
            EXPECT_EQ(1u, elem_graph.get_entity_id_of_remote_element(shell2, 1));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(shell2, 0));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(shell2, 1));
        }
        else if (proc == 2) {
            // Connectivity for Hex Element 2
            EXPECT_EQ(1u, elem_graph.get_num_connected_elems(hex3));
            EXPECT_EQ(4,  elem_graph.get_side_id_to_connected_element(hex3, 0));
            EXPECT_EQ(2u, elem_graph.get_entity_id_of_remote_element(hex3, 0));
            EXPECT_FALSE(elem_graph.is_connected_elem_locally_owned(hex3, 0));
        }

        stk::mesh::EntityId elem_to_move_global_id = 2;
        std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
        stk::mesh::Entity elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_to_move_global_id);

        if (proc == 1)
        {
            int destination_proc = 2;
            elem_proc_pairs_to_move.push_back(std::make_pair(elem_to_move, destination_proc));
        }

        change_entity_owner(bulkData, elem_graph, elem_proc_pairs_to_move);

        std::vector<unsigned> counts;
        stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
        int numLocallyOwnedElemsInMesh = counts[stk::topology::ELEM_RANK];

        size_t size_of_elem_graph = elem_graph.size();

        if (proc == 0)
        {
            EXPECT_EQ(1, numLocallyOwnedElemsInMesh);
            EXPECT_EQ(1u, size_of_elem_graph);

            stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
            impl::parallel_info& elem1_to_elem2_info = elem_graph.get_parallel_edge_info(elem_1, stk::mesh::EntityId(2));
            EXPECT_EQ(2, elem1_to_elem2_info.m_other_proc);
        }
        if (proc == 1)
        {
            EXPECT_EQ(0, numLocallyOwnedElemsInMesh);
            EXPECT_EQ(0u, size_of_elem_graph);

            stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1)), std::logic_error);

            ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(3)), std::logic_error);
        }
        if (proc == 2)
        {
            EXPECT_EQ(2, numLocallyOwnedElemsInMesh);
            EXPECT_EQ(2u, size_of_elem_graph);

            stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
            impl::parallel_info& elem2_to_elem1_info = elem_graph.get_parallel_edge_info(elem_2, stk::mesh::EntityId(1));
            EXPECT_EQ(0, elem2_to_elem1_info.m_other_proc);
        }
    }
}

TEST(ElementGraph, test_change_entity_owner_3_procs_hex_shell_hex_mesh_with_aura)
{
    bool aura_on = true;
    change_entity_owner_hex_shell_hex_test_3_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_3_procs_hex_shell_hex_mesh_without_aura)
{
    bool aura_on = false;
    change_entity_owner_hex_shell_hex_test_3_procs(aura_on);
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

        stk::mesh::EntityVector local_id_to_element_entity(num_locally_owned_elems, Entity());
        std::vector<stk::topology> element_topologies(num_locally_owned_elems);
        impl::set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        ElemElemGraphTester elemElemGraph(bulkData);

        ASSERT_EQ(num_locally_owned_elems, elemElemGraph.get_graph_size());

        SidesForElementGraph & via_sides = elemElemGraph.get_via_sides();

        wall_times.push_back(stk::wall_time());
        msgs.push_back("after ElemElemGraphTester constructor");
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
        bulkData.my_internal_modification_end_for_skin_mesh(stk::topology::FACE_RANK, stk::mesh::impl::MeshModification::MOD_END_SORT, element_selector, NULL);

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

TEST(ElementSide, get_or_create_element_side_with_permutation)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) == 1)
    {
        unsigned spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::mesh::PartVector face_parts {&new_faces_part};
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

        EXPECT_FALSE(impl::is_id_already_in_use_locally(bulkData, side_rank, side_global_id));
        EXPECT_FALSE(impl::does_side_exist_with_different_permutation(bulkData, element1, side_ordinal, side_permutation));
        EXPECT_FALSE(impl::does_element_side_exist(bulkData, element1, side_ordinal));

        impl::connect_side_to_element(bulkData, element1, side_global_id, side_ordinal, side_permutation, face_parts);

        bulkData.modification_end();

        stk::mesh::Entity side1 = get_element_side(bulkData, element1, side_ordinal);
        EXPECT_TRUE(bulkData.is_valid(side1));

        EXPECT_TRUE(impl::is_id_already_in_use_locally(bulkData, side_rank, side_global_id));

        stk::mesh::Permutation side_permutation1 = static_cast<stk::mesh::Permutation>(0);
        EXPECT_TRUE(impl::does_side_exist_with_different_permutation(bulkData, element1, side_ordinal, side_permutation1));

        size_t num_sides = stk::mesh::count_selected_entities(new_faces_part, bulkData.buckets(side_rank));
        EXPECT_EQ(1u, num_sides);
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
        int other_proc = 1;
        int other_side_ord = 2;
        LocalId local_element = 1;
        LocalId other_element = 3;
        int permutation = 0;
        stk::mesh::EntityId chosen_face_id = 1;

        parallel_graph_info.insert(std::make_pair(std::make_pair(local_element, other_element), parallel_info(other_proc, other_side_ord, permutation,
                chosen_face_id)));

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
        size_t numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];

        stk::mesh::EntityVector local_id_to_element_entity(numLocallyOwnedElems, Entity());
        std::vector<stk::topology> element_topologies(numLocallyOwnedElems);
        impl::set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

        ElemElemGraphTester elemElemGraph(bulkData);

        ASSERT_EQ(numLocallyOwnedElems, elemElemGraph.get_graph_size());

        ElementGraph & elem_graph = elemElemGraph.get_element_graph();

        ParallelGraphInfo & parallel_graph_info = elemElemGraph.get_parallel_graph_info();

        if(stk::parallel_machine_rank(comm)==0)
        {
            LocalId local_element = 1;
            LocalId other_element = 3;
            int other_proc = 1;
            int other_side_ord = 4; // 4 left, 5 right
            int permutation = 4;

            test_parallel_graph_info(elem_graph, parallel_graph_info, local_element, other_element, other_proc, other_side_ord, permutation);
        }
        else
        {
            LocalId local_element = 0;
            LocalId other_element = 2;
            int other_proc = 0;
            int other_side_ord = 5; // 4 left, 5 right
            int permutation = 4;

            test_parallel_graph_info(elem_graph, parallel_graph_info, local_element, other_element, other_proc, other_side_ord, permutation);
        }
    }
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

            stk::mesh::EntityVector local_id_to_element_entity(num_locally_owned_elems, Entity());
            std::vector<stk::topology> element_topologies(num_locally_owned_elems);
            impl::set_local_ids_and_fill_element_entities_and_topologies(bulkData, local_id_to_element_entity, element_topologies);

            ElemElemGraphTester elemElemGraph(bulkData);

            SidesForElementGraph & via_sides = elemElemGraph.get_via_sides();

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
            bulkData.my_internal_modification_end_for_skin_mesh(stk::topology::FACE_RANK, stk::mesh::impl::MeshModification::MOD_END_SORT, element_selector, NULL);

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

TEST(ElementGraph, make_items_inactive)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) <= 2)
    {
        unsigned spatialDim = 3;

        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::mesh::PartVector boundary_mesh_parts { &faces_part };
        stk::io::put_io_part_attribute(faces_part);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names

        ASSERT_TRUE(active.primary_entity_rank() == stk::topology::INVALID_RANK);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        ElementDeathUtils::put_mesh_into_part(bulkData, active);

        stk::mesh::ElemElemGraph graph(bulkData);

        size_t num_gold_edges =  6/bulkData.parallel_size();
        ASSERT_EQ(num_gold_edges, graph.num_edges());

        stk::mesh::EntityVector deactivated_elems;

        stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
        if (bulkData.is_valid(elem1) && bulkData.bucket(elem1).owned() )
        {
            deactivated_elems.push_back(elem1);
        }
        stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
        if (bulkData.is_valid(elem3) && bulkData.bucket(elem3).owned())
        {
            deactivated_elems.push_back(elem3);
        }

        active.set_primary_entity_rank(stk::topology::ELEM_RANK);

        bulkData.modification_begin();

        for(size_t i=0; i<deactivated_elems.size(); ++i)
        {
            bulkData.change_entity_parts(deactivated_elems[i], stk::mesh::PartVector(), stk::mesh::PartVector(1, &active));
        }

        bulkData.modification_end();

        active.set_primary_entity_rank(stk::topology::INVALID_RANK);

        for(size_t i=0; i<deactivated_elems.size(); ++i)
        {
            EXPECT_FALSE(bulkData.bucket(deactivated_elems[i]).member(active));
        }

        const stk::mesh::BucketVector& all_node_buckets = bulkData.buckets(stk::topology::NODE_RANK);
        for(size_t i=0; i<all_node_buckets.size(); ++i)
        {
            const stk::mesh::Bucket& bucket = *all_node_buckets[i];
            for(size_t node_index=0; node_index<bucket.size(); ++node_index)
            {
                stk::mesh::EntityId id = bulkData.identifier(bucket[node_index]);
                if (id >=1 && id <= 4)
                {
                    EXPECT_FALSE(bucket.member(active));
                }
                else
                {
                    EXPECT_TRUE(bucket.member(active)) << "for node id " << id << std::endl;
                }
            }
        }
   }
}

// FIXME: Fails with local face ID already in use
TEST(ElementGraph, test_element_death)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) <= 2)
    {
        std::string dimension = unitTestUtils::getOption("--zdim", "none");

        //IO error when this is <4.  Shared face being attached to the wrong element
        int xdim = 4;
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
            stk::mesh::PartVector boundary_mesh_parts { &faces_part };
            stk::io::put_io_part_attribute(faces_part);
            BulkDataElementGraphTester bulkData(meta, comm);

            stk::mesh::Part& active = meta.declare_part("active", stk::topology::ELEMENT_RANK);
            stk::unit_test_util::fill_mesh_using_stk_io(filename, bulkData, comm);
            stk::unit_test_util::write_mesh_using_stk_io("orig.exo", bulkData, bulkData.parallel());

            double start_graph = stk::wall_time();

            ASSERT_TRUE(meta.get_part("block_1") != NULL);

            stk::mesh::Part& block_1 = *meta.get_part("block_1");

            ElementDeathUtils::put_mesh_into_part(bulkData, active);

            std::ostringstream os;
            os << "Proc id: " << bulkData.parallel_rank() << std::endl;

            stk::mesh::ElemElemGraph elementGraph(bulkData);

            double elapsed_graph_time = stk::wall_time() - start_graph;
            os << "Time to create graph: " << elapsed_graph_time << std::endl;

            stk::mesh::EntityRank side_rank = meta.side_rank();

            int num_time_steps = xdim * ydim * zdim;
            double elapsed_death_time = 0;

            for(int i = 0; i < num_time_steps; ++i)
            {
                stk::mesh::EntityVector killedElements = get_killed_elements(bulkData, i, active);
                move_killled_elements_to_part(bulkData, killedElements, block_1, active);
                double start_time = stk::wall_time();
                perform_element_death(bulkData, elementGraph, killedElements, active, boundary_mesh_parts);
                elapsed_death_time += (stk::wall_time() - start_time);
            }

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
                os << "Total time: " << elapsed_death_time << std::endl;
                os << "Total # of alive elements: " << num_active << std::endl;
                os << "Total # of faces: " << num_faces << std::endl;
            }

//                std::ostringstream os;
//                os << bulkData.parallel_rank() << std::endl;
//                {
//                    const stk::mesh::BucketVector &face_buckets = bulkData.buckets(stk::topology::FACE_RANK);
//                    for(size_t i=0;i<face_buckets.size();++i)
//                    {
//                        const stk::mesh::Bucket &bucket = *face_buckets[i];
//                        if(bucket.owned())
//                        {
//                            for(size_t j=0;j<bucket.size();++j)
//                            {
//                                os << "Face " << bulkData.identifier(bucket[j]) << " exists.\n";
//                            }
//                        }
//                    }
//                }
//
//                {
//                    const stk::mesh::BucketVector &buckets = bulkData.buckets(stk::topology::ELEM_RANK);
//                    for(size_t i=0;i<buckets.size();++i)
//                    {
//                        const stk::mesh::Bucket &bucket = *buckets[i];
//                        if(!bucket.member(inactive) && bucket.owned())
//                        {
//                            for(size_t j=0;j<bucket.size();++j)
//                            {
//                                os << "Element " << bulkData.identifier(bucket[j]) << " exists.\n";
//                            }
//                        }
//                    }
//                }
            std::cerr << os.str();
            stk::unit_test_util::write_mesh_using_stk_io("out.exo", bulkData, bulkData.parallel());
        }
    }
}

//EndDocExample1

TEST( ElementGraph, HexHexHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0-----------15.0
    //          /|             /|             /|             /|
    //         / |            / |            / |            / |
    //        /  |           /  |           /  |           /  |
    //      4.0------------8.0-----------12.0-----------16.0  |
    //       |   |          |   |          |   |          |   |
    //       |   |   1.0    |   |   2.0    |   |   3.0    |   |
    //       |   |          |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.0---------|-14.0
    //       |  /           |  /           |  /           |  /
    //       | /            | /            | /            | /
    //       |/             |/             |/             |/
    //      1.0------------5.0------------9.0-----------13.0
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 1u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1,  2,  3,  4,  5,  6,  7,  8 },
        { 5,  6,  7,  8,  9, 10, 11, 12 },
        { 9, 10, 11, 12, 13, 14, 15, 16 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 0, 0 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,    elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,    elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(5,    elemElemGraph.get_side_id_to_connected_element(hex2, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_EQ(hex3, elemElemGraph.get_connected_element(hex2, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Hex Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex3));
    EXPECT_EQ(4,    elemElemGraph.get_side_id_to_connected_element(hex3, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(hex3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
}

TEST( ElementGraph, HexShellSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |   |
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(shell2, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,             elemElemGraph.get_num_connected_elems(shell2));
    EXPECT_EQ(1,              elemElemGraph.get_side_id_to_connected_element(shell2, 0));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
}

TEST( ElementGraph, AdjacentHexShellSerial )
{
    //  ID.proc
    //
    //         12.0-----------11.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      9.0-----------10.0  |
    //       |   |          |   |
    //       |   |   2.0    |4.0|
    //       |   |          |   |
    //       |  8.0---------|--7.0
    //       |  /|          |  /|
    //       | / |          | / |
    //       |/  |          |/  |
    //      5.0------------6.0  |
    //       |   |          |   |
    //       |   |   1.0    |3.0|
    //       |   |          |   |
    //       |  4.0---------|--3.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------2.0
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 2, 3,  7,  6 },
        { 6, 7, 11, 10 },
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(1,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
    EXPECT_EQ(hex2,   elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(1,      elemElemGraph.get_side_id_to_connected_element(hex2, 1));
    EXPECT_EQ(hex1,   elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex2, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

    // Connectivity for Shell Element 4
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell4));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
}

TEST( ElementGraph, DISABLED_HexShellShellSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |3.0|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2, 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
    EXPECT_EQ(shell2, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell2, 0));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
}

TEST( ElementGraph, HexShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Hex Element 2
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell3, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
}

TEST( ElementGraph, DISABLED_HexShellShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex1, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 1));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex2, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell3, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell4, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell4, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell4, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
}

TEST( ElementGraph, DISABLED_HexShellReversedShellHexSerial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.0    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Two stacked shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int p_size = stk::parallel_machine_size(pm);

    if(p_size > 1)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 8, 7, 6 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
    EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex1, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
    EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 1));
    EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex2, 0));
    EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex2, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell3, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
    EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 1));
    EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell4, 0));
    EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell4, 1));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
}

void setup_node_sharing(stk::mesh::BulkData &mesh, const std::vector< std::vector<unsigned> > & shared_nodeIDs_and_procs )
{
    const unsigned p_rank = mesh.parallel_rank();

    for (size_t nodeIdx = 0, end = shared_nodeIDs_and_procs.size(); nodeIdx < end; ++nodeIdx) {
        if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
            stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
            int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
            mesh.add_node_sharing(node, sharingProc);
        }
    }
}

TEST( ElementGraph, Hex0Hex0Hex1Serial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0-----------15.1
    //          /|             /|             /|             /|
    //         / |            / |            / |            / |
    //        /  |           /  |           /  |           /  |
    //      4.0------------8.0-----------12.0-----------16.1  |
    //       |   |          |   |          |   |          |   |
    //       |   |   1.0    |   |   2.0    |   |   3.1    |   |
    //       |   |          |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.0---------|-14.1
    //       |  /           |  /           |  /           |  /
    //       | /            | /            | /            | /
    //       |/             |/             |/             |/
    //      1.0------------5.0------------9.0-----------13.1
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1,  2,  3,  4,  5,  6,  7,  8 },
        { 5,  6,  7,  8,  9, 10, 11, 12 },
        { 9, 10, 11, 12, 13, 14, 15, 16 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 0, 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0,  9, 1 },  // proc 0
        { 0, 10, 1 },
        { 0, 11, 1 },
        { 0, 12, 1 },
        { 1,  9, 0 },  // proc 1
        { 1, 10, 0 },
        { 1, 11, 0 },
        { 1, 12, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,    elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(hex2, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

        // Connectivity for Hex Element 2
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,    elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(5,    elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(hex2, 0));
        EXPECT_EQ(3u,   elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 3
        EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex3));
        EXPECT_EQ(4,    elemElemGraph.get_side_id_to_connected_element(hex3, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(hex3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
    }
}

TEST( ElementGraph, Hex0Hex1Hex0Serial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0-----------15.0
    //          /|             /|             /|             /|
    //         / |            / |            / |            / |
    //        /  |           /  |           /  |           /  |
    //      4.0------------8.0-----------12.0-----------16.0  |
    //       |   |          |   |          |   |          |   |
    //       |   |   1.0    |   |   2.1    |   |   3.0    |   |
    //       |   |          |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.0---------|-14.0
    //       |  /           |  /           |  /           |  /
    //       | /            | /            | /            | /
    //       |/             |/             |/             |/
    //      1.0------------5.0------------9.0-----------13.0
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1,  2,  3,  4,  5,  6,  7,  8 },
        { 5,  6,  7,  8,  9, 10, 11, 12 },
        { 9, 10, 11, 12, 13, 14, 15, 16 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0,  5, 1 },  // proc 0
        { 0,  6, 1 },
        { 0,  7, 1 },
        { 0,  8, 1 },
        { 0,  9, 1 },
        { 0, 10, 1 },
        { 0, 11, 1 },
        { 0, 12, 1 },
        { 1,  5, 0 },  // proc 1
        { 1,  6, 0 },
        { 1,  7, 0 },
        { 1,  8, 0 },
        { 1,  9, 0 },
        { 1, 10, 0 },
        { 1, 11, 0 },
        { 1, 12, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

        // Connectivity for Hex Element 3
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex3));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex3, 0));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(hex3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
}

TEST( ElementGraph, Hex0Hex1Hex2Serial )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1-----------15.2
    //          /|             /|             /|             /|
    //         / |            / |            / |            / |
    //        /  |           /  |           /  |           /  |
    //      4.0------------8.0-----------12.1-----------16.2  |
    //       |   |          |   |          |   |          |   |
    //       |   |   1.0    |   |   2.1    |   |   3.2    |   |
    //       |   |          |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.1---------|-14.2
    //       |  /           |  /           |  /           |  /
    //       | /            | /            | /            | /
    //       |/             |/             |/             |/
    //      1.0------------5.0------------9.1-----------13.2
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 3u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1,  2,  3,  4,  5,  6,  7,  8 },
        { 5,  6,  7,  8,  9, 10, 11, 12 },
        { 9, 10, 11, 12, 13, 14, 15, 16 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1, 2 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0,  5, 1 },  // proc 0
        { 0,  6, 1 },
        { 0,  7, 1 },
        { 0,  8, 1 },
        { 1,  5, 0 },  // proc 1
        { 1,  6, 0 },
        { 1,  7, 0 },
        { 1,  8, 0 },
        { 1,  9, 2 },
        { 1, 10, 2 },
        { 1, 11, 2 },
        { 1, 12, 2 },
        { 2,  9, 1 },  // proc 2
        { 2, 10, 1 },
        { 2, 11, 1 },
        { 2, 12, 1 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
    else if (p_rank == 2) {
        // Connectivity for Hex Element 3
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex3));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex3, 0));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(hex3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
    }
}

TEST( ElementGraph, Hex0Shell1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.1|
    //       |   |          |   |
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2 };
    stk::mesh::EntityId shellElemOwningProc[] = { 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 2
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell2, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
    }
}

TEST( ElementGraph, AdjacentHex0Shell1Parallel )
{
    //  ID.proc
    //
    //         12.0-----------11.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      9.0-----------10.0  |
    //       |   |          |   |
    //       |   |   2.0    |4.1|
    //       |   |          |   |
    //       |  8.0---------|--7.0
    //       |  /|          |  /|
    //       | / |          | / |
    //       |/  |          |/  |
    //      5.0------------6.0  |
    //       |   |          |   |
    //       |   |   1.0    |3.1|
    //       |   |          |   |
    //       |  4.0---------|--3.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------2.0
    //                      ^
    //                      |
    //                       ---- Two adjacent shell elements

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 0 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 2, 3,  7,  6 },
        { 6, 7, 11, 10 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 1, 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0,  2, 1 },  // proc 0
        { 0,  3, 1 },
        { 0,  7, 1 },
        { 0,  6, 1 },
        { 0, 11, 1 },
        { 0, 10, 1 },
        { 1,  2, 0 },  // proc 1
        { 1,  3, 0 },
        { 1,  7, 0 },
        { 1,  6, 0 },
        { 1, 11, 0 },
        { 1, 10, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0u) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,    elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(hex1, 1));
        EXPECT_EQ(hex2, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_EQ(3u,   elemElemGraph.get_entity_id_of_remote_element(hex1, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Hex Element 2
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,    elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(hex2, 0));
        EXPECT_EQ(4u,   elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
    else if (p_rank == 1u) {
        // Connectivity for Shell Element 3
        EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(1u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

        // Connectivity for Shell Element 4
        EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell4));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    }
}

TEST( ElementGraph, DISABLED_Hex0Shell0Shell1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.0|
    //       |   |          |3.1|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2, 3 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(shell2, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_EQ(3u,     elemElemGraph.get_entity_id_of_remote_element(hex1, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 2
        EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell2, 0));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell2, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 3
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    }
}

TEST( ElementGraph, DISABLED_Hex0Shell1Shell1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0
    //          /|             /|
    //         / |            / |
    //        /  |           /  |
    //      4.0------------8.0  |
    //       |   |          |   |
    //       |   |   1.0    |2.1|
    //       |   |          |3.1|
    //       |  2.0---------|--6.0
    //       |  /           |  /
    //       | /            | /
    //       |/             |/
    //      1.0------------5.0
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 2, 3 };
    stk::mesh::EntityId shellElemOwningProc[] = { 1, 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 1));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex1, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 2
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell2, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

        // Connectivity for Shell Element 3
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    }
}

TEST( ElementGraph, Hex0Shell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    }
}

TEST( ElementGraph, Hex0Shell1Hex2Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.2
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.2  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.1|   2.2    |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.2
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.2
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 3u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3 };
    stk::mesh::EntityId shellElemOwningProc[] = { 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 0, 5, 2 },
        { 0, 6, 2 },
        { 0, 7, 2 },
        { 0, 8, 2 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 },
        { 1, 5, 2 },
        { 1, 6, 2 },
        { 1, 7, 2 },
        { 1, 8, 2 },
        { 2, 5, 0 },  // proc 2
        { 2, 6, 0 },
        { 2, 7, 0 },
        { 2, 8, 0 },
        { 2, 5, 1 },
        { 2, 6, 1 },
        { 2, 7, 1 },
        { 2, 8, 1 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 3
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(0,  elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(shell3, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    }
    else if (p_rank == 2) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    }
}

TEST( ElementGraph, Hex0Shell1Hex0Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.0
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.0  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.1|   2.0    |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.0
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.0
    //                      ^
    //                      |
    //                       ---- Single shell element

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 0 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3 };
    stk::mesh::EntityId shellElemOwningProc[] = { 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

        // Connectivity for Hex Element 2
        EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 0));
        EXPECT_EQ(1u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    }
}

TEST( ElementGraph, DISABLED_Hex0Shell0Shell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
        EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex1, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

        // Connectivity for Shell Element 4
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell4, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell4, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell4, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(3u,   elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_EQ(4u,   elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
}

TEST( ElementGraph, DISABLED_Hex0Shell0Shell1Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |4.1|          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 1 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
        EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_EQ(4u,     elemElemGraph.get_entity_id_of_remote_element(hex1, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Shell Element 4
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 1));
        EXPECT_EQ(hex2, elemElemGraph.get_connected_element(shell4, 0));
        EXPECT_EQ(1u,   elemElemGraph.get_entity_id_of_remote_element(shell4, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

        // Connectivity for Hex Element 2
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(4,      elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex2, 0));
        EXPECT_EQ(3u,     elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
}

TEST( ElementGraph, DISABLED_Hex0Shell0ReversedShell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.0    |3.0|   2.1    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells, opposite orientation

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 8, 7, 6 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(5,      elemElemGraph.get_side_id_to_connected_element(hex1, 1));
        EXPECT_EQ(shell3, elemElemGraph.get_connected_element(hex1, 0));
        EXPECT_EQ(shell4, elemElemGraph.get_connected_element(hex1, 1));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Shell Element 3
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell3, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

        // Connectivity for Shell Element 4
        EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
        EXPECT_EQ(0,    elemElemGraph.get_side_id_to_connected_element(shell4, 0));
        EXPECT_EQ(1,    elemElemGraph.get_side_id_to_connected_element(shell4, 1));
        EXPECT_EQ(hex1, elemElemGraph.get_connected_element(shell4, 0));
        EXPECT_EQ(2u,   elemElemGraph.get_entity_id_of_remote_element(shell4, 1));
        EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 2
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_EQ(4u, elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
}

TEST( ElementGraph, DISABLED_Hex1Shell0Shell0Hex1Parallel )
{
    //  ID.proc
    //
    //          3.0------------7.1-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.1-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1.1    |3.0|   2.1    |   |
    //       |   |          |4.0|          |   |
    //       |  2.0---------|--6.1---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.1------------9.1
    //                      ^
    //                      |
    //                       ---- Two stacked shells

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 2u)
    {
        return;
    }

    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
        { 1, 2, 3, 4, 5,  6,  7,  8 },
        { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
    stk::mesh::EntityId hexElemOwningProc[] = { 1, 1 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
        { 5, 6, 7, 8 },
        { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
    stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

    // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
    {
        { 0, 5, 1 },  // proc 0
        { 0, 6, 1 },
        { 0, 7, 1 },
        { 0, 8, 1 },
        { 1, 5, 0 },  // proc 1
        { 1, 6, 0 },
        { 1, 7, 0 },
        { 1, 8, 0 }
    };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
        if (hexElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
        }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
        if (shellElemOwningProc[i] == p_rank) {
            stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
        }
    }
    setup_node_sharing(mesh, shared_nodeIDs_and_procs );
    mesh.modification_end();

    stk::mesh::ElemElemGraph elemElemGraph(mesh);

    const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
    const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
    const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
    const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

    if (p_rank == 0) {
        // Connectivity for Shell Element 3
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
        EXPECT_EQ(0,  elemElemGraph.get_side_id_to_connected_element(shell3, 0));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell3, 1));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(shell3, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell3, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

        // Connectivity for Shell Element 4
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell4));
        EXPECT_EQ(0,  elemElemGraph.get_side_id_to_connected_element(shell4, 0));
        EXPECT_EQ(1,  elemElemGraph.get_side_id_to_connected_element(shell4, 1));
        EXPECT_EQ(2u, elemElemGraph.get_entity_id_of_remote_element(shell4, 0));
        EXPECT_EQ(1u, elemElemGraph.get_entity_id_of_remote_element(shell4, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
    }
    else if (p_rank == 1) {
        // Connectivity for Hex Element 1
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 0));
        EXPECT_EQ(5,  elemElemGraph.get_side_id_to_connected_element(hex1, 1));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex1, 0));
        EXPECT_EQ(4u, elemElemGraph.get_entity_id_of_remote_element(hex1, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

        // Connectivity for Hex Element 2
        EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 0));
        EXPECT_EQ(4,  elemElemGraph.get_side_id_to_connected_element(hex2, 1));
        EXPECT_EQ(3u, elemElemGraph.get_entity_id_of_remote_element(hex2, 0));
        EXPECT_EQ(4u, elemElemGraph.get_entity_id_of_remote_element(hex2, 1));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
        EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
}


// element ids / proc_id:
// |-------|-------|-------|
// |       |       |       |
// |  1/0  |  4/2  |  7/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  2/0  |  5/1  |  8/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  3/0  |  6/2  |  9/2  |
// |       |       |       |
// |-------|-------|-------|


TEST(ElementGraph, TestKeyHoleSimilarProblemAInParallel)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) == 3)
    {
        const int procRank = stk::parallel_machine_rank(comm);
        unsigned spatialDim = 3;

        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::BulkData bulkData(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:3x1x3", bulkData, comm);

        stk::mesh::EntityProcVec elementProcChanges;
        if (procRank == 1) {
            elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,4),2));
            elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,6),2));
        }
        bulkData.change_entity_owner(elementProcChanges);

        stk::mesh::ElemElemGraph graph(bulkData);
        if (procRank == 0) {
            stk::mesh::Entity local_element = bulkData.get_entity(stk::topology::ELEM_RANK,2);
            ASSERT_TRUE(bulkData.bucket(local_element).owned());
            ASSERT_EQ(3u, graph.get_num_connected_elems(local_element));

            EXPECT_EQ( 3, graph.get_side_id_to_connected_element(local_element,0));
            EXPECT_EQ( 1, graph.get_side_id_to_connected_element(local_element,1));
            EXPECT_EQ( 5, graph.get_side_id_to_connected_element(local_element,2));

            EXPECT_TRUE(graph.is_connected_elem_locally_owned(local_element, 0));
            EXPECT_TRUE(graph.is_connected_elem_locally_owned(local_element, 1));
            EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 2));

            EXPECT_EQ( 1u, bulkData.identifier(graph.get_connected_element(local_element, 0)));
            EXPECT_EQ( 3u, bulkData.identifier(graph.get_connected_element(local_element, 1)));
            EXPECT_EQ( 5u, graph.get_entity_id_of_remote_element(local_element, 2));

            EXPECT_EQ( 1, graph.get_owning_proc_id_of_remote_element(local_element, 5));
        }
        if (procRank == 1) {
            stk::mesh::Entity local_element = bulkData.get_entity(stk::topology::ELEM_RANK,5);
            ASSERT_TRUE(bulkData.bucket(local_element).owned());
            size_t numConnectedElems = graph.get_num_connected_elems(local_element);
            ASSERT_EQ(4u, numConnectedElems);

            EXPECT_EQ( 4, graph.get_side_id_to_connected_element(local_element,0));
            EXPECT_EQ( 3, graph.get_side_id_to_connected_element(local_element,1));
            EXPECT_EQ( 1, graph.get_side_id_to_connected_element(local_element,2));
            EXPECT_EQ( 5, graph.get_side_id_to_connected_element(local_element,3));

            EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 0));
            EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 1));
            EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 2));
            EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 3));

            EXPECT_EQ( 2u, graph.get_entity_id_of_remote_element(local_element, 0));
            EXPECT_EQ( 4u, graph.get_entity_id_of_remote_element(local_element, 1));
            EXPECT_EQ( 6u, graph.get_entity_id_of_remote_element(local_element, 2));
            EXPECT_EQ( 8u, graph.get_entity_id_of_remote_element(local_element, 3));

            EXPECT_EQ( 0, graph.get_owning_proc_id_of_remote_element(local_element, 2));
            EXPECT_EQ( 2, graph.get_owning_proc_id_of_remote_element(local_element, 8));
            EXPECT_EQ( 2, graph.get_owning_proc_id_of_remote_element(local_element, 4));
            EXPECT_EQ( 2, graph.get_owning_proc_id_of_remote_element(local_element, 6));
        }
    }
}

// element ids / proc_id:
// |-------|-------|-------|
// |       |       |       |
// |  1/0  |  4/2  |  7/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  2/0  |  n/a  |  8/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  3/0  |  6/2  |  9/2  |
// |       |       |       |
// |-------|-------|-------|
// The element in the middle has been deleted

TEST(ElementGraph, TestKeyHoleSimilarProblemBInParallel)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) == 3)
    {
        const int procRank = stk::parallel_machine_rank(comm);
        unsigned spatialDim = 3;

        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::BulkData bulkData(meta, comm);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:3x1x3", bulkData, comm);

        stk::mesh::EntityProcVec elementProcChanges;
        if (procRank == 1) {
            elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,4),2));
            elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,6),2));
        }
        bulkData.change_entity_owner(elementProcChanges);
        bulkData.modification_begin();
        if (procRank == 1) {
            stk::mesh::Entity local_element5 = bulkData.get_entity(stk::topology::ELEM_RANK,5);
            bulkData.destroy_entity(local_element5);
        }
        bulkData.modification_end();

        stk::mesh::ElemElemGraph graph(bulkData);
        if (procRank == 0) {
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,1)));
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,2)));
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,3)));
        }
        if (procRank == 2) {
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,4)));
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,6)));
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,7)));
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,8)));
            EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,9)));
        }
    }
}



bool is_valid_graph_element(const impl::ElementGraph &elem_graph, stk::mesh::impl::LocalId elem_id)
{
    stk::mesh::impl::LocalId max_elem_id = static_cast<stk::mesh::impl::LocalId>(elem_graph.size());
    return (elem_id >= 0 && elem_id < max_elem_id);
}

int check_connectivity(const impl::ElementGraph &elem_graph, const impl::SidesForElementGraph &via_sides,
                       stk::mesh::impl::LocalId element_id1, stk::mesh::impl::LocalId element_id2)
{
    int side=-1;
    if (is_valid_graph_element(elem_graph, element_id1) && is_valid_graph_element(elem_graph, element_id2)) {
        side = get_side_from_element1_to_element2(elem_graph, via_sides, element_id1, element_id2);
    }
    return side;
}

int get_side_from_element1_to_element2(const impl::ElementGraph &elem_graph,
                                       const impl::SidesForElementGraph &via_sides,
                                       stk::mesh::impl::LocalId element1_local_id,
                                       stk::mesh::impl::LocalId other_element_id)
{
    int side = -1;
    const std::vector<stk::mesh::impl::LocalId>& conn_elements = elem_graph[element1_local_id];

    std::vector<stk::mesh::impl::LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), other_element_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = via_sides[element1_local_id][index];
    }
    return side;
}



void pack_move_elements_and_via_sides_to_comm(stk::CommSparse &comm, stk::mesh::BulkData &bulkData, ElementGraph &elem_graph, SidesForElementGraph &via_sides, std::vector<LocalId> &elems_to_move)
{
    const unsigned num_elems_to_move = elems_to_move.size();
    std::vector<int> sharing_procs;
    for (unsigned i = 0; i < num_elems_to_move; ++i)
    {
        sharing_procs.clear();
        stk::mesh::EntityKey elem_key = bulkData.entity_key(bulkData.get_entity(stk::topology::ELEM_RANK, elems_to_move[i]));
        bulkData.comm_shared_procs(elem_key, sharing_procs);
        LocalId elem_id = elems_to_move[i];
        for (unsigned j = 0; j < elem_graph[elem_id].size(); ++j)
        {
            for(size_t proc_index=0; proc_index<sharing_procs.size(); ++proc_index)
            {
                if (elem_graph[elem_id][j] >= 0)
                {
                    stk::mesh::Entity connected_elem = bulkData.get_entity(stk::topology::ELEM_RANK, elem_graph[elem_id][j]);
                    stk::mesh::EntityId locally_connected_elem_gid = -bulkData.identifier(connected_elem);
                    int locally_connected_side_id = via_sides[elem_id][j];
                    comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityId>(locally_connected_elem_gid);
                    comm.send_buffer(sharing_procs[proc_index]).pack<int>(locally_connected_side_id);
                }
                else
                {
                    stk::mesh::EntityId remote_connected_elem_gid = -elem_graph[elem_id][j];
                    int remote_connected_side_id = via_sides[elem_id][j];
                    comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityId>(remote_connected_elem_gid);
                    comm.send_buffer(sharing_procs[proc_index]).pack<int>(remote_connected_side_id);
                }
            }
        }
    }
}

void move_elems_to_remote_proc(stk::mesh::BulkData &bulkData, ElementGraph &elem_graph, SidesForElementGraph &via_sides, std::vector<LocalId> &elems_to_move, int from_proc, int to_proc)
{
    stk::CommSparse comm(bulkData.parallel());

    pack_move_elements_and_via_sides_to_comm(comm, bulkData, elem_graph, via_sides, elems_to_move);
    comm.allocate_buffers();

    pack_move_elements_and_via_sides_to_comm(comm, bulkData, elem_graph, via_sides, elems_to_move);
    comm.communicate();

    for(int proc_id=0; proc_id<bulkData.parallel_size(); ++proc_id)
    {
        if (proc_id != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityId element_id;
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(element_id);
                int side_index = 0;
                comm.recv_buffer(proc_id).unpack<int>(side_index);

            }
        }
    }
}
//EndDocExample1

}
