
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, ostream, cerr, etc
#include <stk_io/IossBridge.hpp>        // for put_io_part_attribute
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>  // for process_killed_elements, etc
#include <stk_mesh/baseImpl/elementGraph/ProcessKilledElements.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <vector>                       // for vector
#include "mpi.h"                        // for ompi_communicator_t, etc
#include "stk_mesh/base/Entity.hpp"     // for operator<<, Entity
#include "stk_mesh/base/Selector.hpp"   // for Selector, etc
#include "stk_mesh/base/Types.hpp"      // for EntityVector, etc
#include "stk_unit_test_utils/ElemGraphTestUtils.hpp"
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk { namespace mesh { class Part; } }

namespace
{
using stk::unit_test_util::build_mesh;
using stk::mesh::MeshBuilder;

class BulkDataTester : public stk::mesh::BulkData
{
public:
  BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
    stk::mesh::BulkData(std::shared_ptr<stk::mesh::MetaData>(&mesh_meta_data, [](auto pointerWeWontDelete){}), comm)
  {
  }

  virtual ~BulkDataTester()
  {
  }

  void set_sorting_by_face()
  {
#ifdef SIERRA_MIGRATION
    m_shouldSortFacesByNodeIds = true;
#endif
  }
};


void kill_element(stk::mesh::Entity element, stk::mesh::BulkData& bulkData, stk::mesh::Part& active, stk::mesh::Part& skin)
{
  stk::mesh::ElemElemGraph graph(bulkData);
  stk::mesh::EntityVector deactivated_elems;
  if(bulkData.is_valid(element) && bulkData.parallel_owner_rank(element) == bulkData.parallel_rank())
  {
    deactivated_elems.push_back(element);
  }

  stk::mesh::PartVector boundary_mesh_parts={&active, &skin};
  ElemGraphTestUtils::deactivate_elements(deactivated_elems, bulkData,  active);

  stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
  stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData, graph, active, remoteActiveSelector);

  stk::mesh::process_killed_elements(bulkData, deactivated_elems, active, remoteActiveSelector, boundary_mesh_parts, &boundary_mesh_parts);
}

stk::mesh::EntityVector get_entities(stk::mesh::BulkData& bulkData, const stk::mesh::ConstPartVector& parts)
{
  stk::mesh::EntityVector entities;
  stk::mesh::Selector sel = stk::mesh::selectIntersection(parts);
  stk::mesh::get_selected_entities(sel, bulkData.buckets(stk::topology::FACE_RANK), entities);
  return entities;
}

void compare_faces(const stk::mesh::BulkData& bulkData, const std::vector<size_t> &num_gold_skinned_faces, const stk::mesh::EntityVector& skinned_faces, const stk::mesh::EntityVector &active_faces)
{
  EXPECT_EQ(num_gold_skinned_faces[bulkData.parallel_rank()], skinned_faces.size());
  EXPECT_EQ(num_gold_skinned_faces[bulkData.parallel_rank()], active_faces.size());

  for(size_t i=0;i<skinned_faces.size();++i)
  {
    EXPECT_EQ(bulkData.identifier(skinned_faces[i]), bulkData.identifier(active_faces[i]));
  }
}

void compare_skin(const std::vector<size_t>& num_gold_skinned_faces, stk::mesh::BulkData& bulkData, const stk::mesh::Part& skin, const stk::mesh::Part& active)
{
  stk::mesh::EntityVector skinned_faces = get_entities(bulkData, {&skin, &active});
  stk::mesh::EntityVector active_faces  = get_entities(bulkData, {&active} );
  compare_faces(bulkData, num_gold_skinned_faces, skinned_faces, active_faces);
}

TEST(ElementDeath, compare_death_and_skin_mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) == 1)
  {
    unsigned spatialDim = 3;

    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& skin  = meta.declare_part_with_topology("skin", stk::topology::QUAD_4);
    stk::io::put_io_part_attribute(skin);
    BulkDataTester bulkData(meta, comm);
    bulkData.set_sorting_by_face();

    stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names
    stk::io::fill_mesh("generated:1x1x4", bulkData);
    stk::unit_test_util::put_mesh_into_part(bulkData, active);

    ElemGraphTestUtils::skin_boundary(bulkData, active, {&skin, &active});

    std::vector<size_t> num_gold_skinned_faces = { 18 };
    compare_skin(num_gold_skinned_faces, bulkData, skin, active);

    stk::mesh::Entity element1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
    kill_element(element1, bulkData, active, skin);

    ElemGraphTestUtils::skin_part(bulkData, active, {&skin, &active});

    num_gold_skinned_faces[0] = 14;
    compare_skin(num_gold_skinned_faces, bulkData, skin, active);
  }
}


} // end namespace


