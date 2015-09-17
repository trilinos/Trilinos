#ifndef elementgraphtesterhpp
#define elementgraphtesterhpp

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
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "UnitTestElementDeathUtils.hpp"
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_tests/stk_mesh/SetupKeyholeMesh.hpp>

#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/heterogeneous_mesh.hpp>
#include <stk_mesh/fixtures/degenerate_mesh.hpp>

class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
    ElemElemGraphTester(stk::mesh::BulkData& bulkData)
      : ElemElemGraph(bulkData, bulkData.mesh_meta_data().universal_part()) {};

    ElemElemGraphTester(stk::mesh::BulkData& bulkData, stk::mesh::Part &part)
      : ElemElemGraph(bulkData, part) {};

    virtual ~ElemElemGraphTester() {}

    void fill_graph() { ElemElemGraph::fill_graph(); }

    void fill_parallel_graph(stk::mesh::impl::ElemSideToProcAndFaceId& elem_side_comm) { ElemElemGraph::fill_parallel_graph(elem_side_comm); }

    stk::mesh::impl::ElementGraph & get_element_graph() { return m_elem_graph; }
    stk::mesh::impl::SidesForElementGraph & get_via_sides() { return m_via_sides; }
    stk::mesh::impl::ParallelGraphInfo & get_parallel_graph_info() { return m_parallel_graph_info; }
    stk::mesh::BulkData & get_bulk_data() { return m_bulk_data; }
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
    int check_connectivity(stk::mesh::EntityId elem1_id, stk::mesh::EntityId elem2_id)
    {
        int side = -1;
        stk::mesh::Entity elem1 = m_bulk_data.get_entity(stk::topology::ELEM_RANK, elem1_id);
        stk::mesh::Entity elem2 = m_bulk_data.get_entity(stk::topology::ELEM_RANK, elem2_id);
        bool isElem1Local = m_bulk_data.is_valid(elem1) && m_bulk_data.bucket(elem1).owned();
        bool isElem2Local = m_bulk_data.is_valid(elem2) && m_bulk_data.bucket(elem2).owned();

        ThrowRequire(isElem1Local);

        if(isElem2Local)
        {
            side = check_local_connectivity(elem1, elem2);
        }
        else
        {
            side = check_remote_connectivity(elem1, elem2_id);
        }

        return side;
    }

    std::vector<stk::mesh::EntityId> get_copy_of_all_ids() const
    {
        return m_suggested_side_ids;
    }
};

#endif
