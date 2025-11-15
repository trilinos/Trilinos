#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_tools/FixNodeSharingViaSearch.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include <vector>
#include <string>
#include <algorithm>

using stk::unit_test_util::build_mesh;

TEST(FixNodeSharingViaSearch, twoHex_2proc)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }
  const int myProc = stk::parallel_machine_rank(comm);

  constexpr int dim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(dim, stk::parallel_machine_world());
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part& hexPart = meta.declare_part_with_topology("hexPart", stk::topology::HEX_8);
  stk::mesh::Field<double>& coordField = meta.declare_field<double>(stk::topology::NODE_RANK, "coords");
  const double* initValue = nullptr;
  stk::mesh::put_field_on_mesh(coordField, meta.universal_part(), dim, initValue);
  meta.set_coordinate_field(&coordField);

  stk::mesh::EntityId elemId = myProc==0 ? 1 : 2;
  stk::mesh::EntityIdVector nodeIDs = (myProc==0) ?
    stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8} : stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12};

  bulkPtr->modification_begin();

  stk::mesh::Entity elem = stk::mesh::declare_element(*bulkPtr, hexPart, elemId, nodeIDs);

  std::vector<double> coords = (myProc==0) ?
   std::vector<double>{0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0,
                       0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0, 0.0,1.0,1.0}
   :
   std::vector<double>{0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0, 0.0,1.0,1.0,
                       0.0,0.0,2.0, 1.0,0.0,2.0, 1.0,1.0,2.0, 0.0,1.0,2.0};

  unsigned counter = 0;
  stk::mesh::ConnectedEntities nodes = bulkPtr->get_connected_entities(elem, stk::topology::NODE_RANK);
  auto coordData = coordField.data<stk::mesh::ReadWrite>();
  for(unsigned n=0; n<nodes.size(); ++n) {
    auto nodeCoords = coordData.entity_values(nodes[n]);
    for(stk::mesh::ComponentIdx d=0_comp; d<dim; ++d) {
      nodeCoords(d) = coords[counter++];
    }
  }

  stk::tools::fix_node_sharing_via_search(*bulkPtr);

  EXPECT_NO_THROW(bulkPtr->modification_end());

  const unsigned expectedNumSharedNodes = 4;
  EXPECT_EQ(expectedNumSharedNodes, stk::mesh::count_entities(*bulkPtr, stk::topology::NODE_RANK, meta.globally_shared_part()));
}

TEST(FixNodeSharingViaSearch, oneD_userTolerance)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }
  const int myProc = stk::parallel_machine_rank(comm);

  constexpr int dim = 1;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(dim, stk::parallel_machine_world());
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part& linePart = meta.declare_part_with_topology("linePart", stk::topology::LINE_2_1D);
  stk::mesh::Field<double>& coordField = meta.declare_field<double>(stk::topology::NODE_RANK, "coords");
  const double* initValue = nullptr;
  stk::mesh::put_field_on_mesh(coordField, meta.universal_part(), dim, initValue);
  meta.set_coordinate_field(&coordField);

  stk::mesh::EntityId elemId = myProc==0 ? 1 : 2;
  stk::mesh::EntityIdVector nodeIDs = (myProc==0) ?
    stk::mesh::EntityIdVector{1,2} : stk::mesh::EntityIdVector{2,3};

  bulkPtr->modification_begin();

  stk::mesh::Entity elem = stk::mesh::declare_element(*bulkPtr, linePart, elemId, nodeIDs);

  std::vector<double> coords = (myProc==0) ?
   std::vector<double>{0.0, 1.0}
   :
   std::vector<double>{1.0001, 2.0};

  unsigned counter = 0;
  stk::mesh::ConnectedEntities nodes = bulkPtr->get_connected_entities(elem, stk::topology::NODE_RANK);
  auto coordFieldData = coordField.data<stk::mesh::ReadWrite>();
  for(unsigned n=0; n<nodes.size(); ++n) {
    auto nodeCoords = coordFieldData.entity_values(nodes[n]);
    for(stk::mesh::ComponentIdx d=0_comp; d<dim; ++d) {
      nodeCoords(d) = coords[counter++];
    }
  }

  constexpr double tolerance = 1.e-2;
  stk::tools::fix_node_sharing_via_search(*bulkPtr, tolerance);

  EXPECT_NO_THROW(bulkPtr->modification_end());

  const unsigned expectedNumSharedNodes = 1;
  EXPECT_EQ(expectedNumSharedNodes, stk::mesh::count_entities(*bulkPtr, stk::topology::NODE_RANK, meta.globally_shared_part()));
}

