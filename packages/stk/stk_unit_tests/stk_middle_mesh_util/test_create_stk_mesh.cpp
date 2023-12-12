
#include <stk_middle_mesh_util/create_stk_mesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include "stk_middle_mesh/utils.hpp"
#include "gtest/gtest.h"

#include <iostream>

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

void read_stk_mesh(const std::string& fname, stk::mesh::BulkData& bulkData)
{
  stk::io::StkMeshIoBroker reader(bulkData.parallel());
  reader.set_bulk_data(bulkData);
  reader.add_mesh_database(fname, stk::io::READ_MESH);
  reader.create_input_mesh();
  reader.add_all_mesh_fields_as_input_fields();
  reader.populate_mesh();
  reader.populate_field_data();
}

void verify_correct_sharing(const stk::mesh::BulkData& stkMesh,
                            const mesh::Mesh& surfaceMesh,
                            const stk_interface::MeshPart::MeshFieldPtr& stkElsFieldPtr)
{
  std::vector<int> sharingProcs;
  constexpr unsigned maxNumEdgeNodes = 3;
  std::vector<stk::mesh::Entity> edgeNodes(maxNumEdgeNodes);

  const std::vector<mesh::MeshEntityPtr>& surfaceElems = surfaceMesh.get_elements();
  for(const mesh::MeshEntityPtr& elem : surfaceElems) {
    if (elem) {
      const stk::mesh::SideSetEntry& ssetEntry = (*stkElsFieldPtr)(elem, 0, 0);
      stk::mesh::Entity stkEl = ssetEntry.element;

      const bool stkElemIsFace = ssetEntry.side != stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
      if (stkElemIsFace) {
        stkEl = stk::mesh::get_side_entity_for_elem_side_pair(stkMesh, stkEl, ssetEntry.side);
      }

      stk::topology stkTopo = stkMesh.bucket(stkEl).topology();

      if (stkElemIsFace) {
        EXPECT_FALSE(stkTopo.is_shell());
      }
      else {
        EXPECT_TRUE(stkTopo.is_shell());
      }

      const stk::mesh::Entity* nodes = stkMesh.begin_nodes(stkEl);

      for(int dn=0; dn<elem->count_down(); ++dn) {
        mesh::MeshEntityPtr edgeEnt = elem->get_down(dn);
        EXPECT_TRUE((edgeEnt && edgeEnt->get_type() == mesh::MeshEntityType::Edge));

        mesh::EntityOrientation orientation = elem->get_down_orientation(dn);

        stkTopo.edge_nodes(nodes, dn, edgeNodes.data());

        int numSharedVerts = 0;
        for(int n=0; n<edgeEnt->count_down(); ++n) {
          int vertIdx = orientation==mesh::EntityOrientation::Standard ? n : (edgeEnt->count_down() - n - 1);
          mesh::MeshEntityPtr vert = edgeEnt->get_down(vertIdx);
          EXPECT_TRUE((vert && vert->get_type() == mesh::MeshEntityType::Vertex));

          stkMesh.comm_shared_procs(edgeNodes[n], sharingProcs);
          EXPECT_EQ(static_cast<unsigned>(vert->count_remote_shared_entities()), sharingProcs.size())<<"dn="<<dn<<", n="<<n;

          if (vert->count_remote_shared_entities() > 0) {
            ++numSharedVerts;
          }
        }

        if (numSharedVerts == edgeEnt->count_down()) {
          EXPECT_TRUE((edgeEnt->count_remote_shared_entities() > 0));
        }
      }
    }
  }
}
}

TEST(StkMeshCreator, create_mesh_from_sideset_part)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 2) { GTEST_SKIP(); }

  stk_interface::StkMeshCreator stkMeshCreator("generated:1x2x2|sideset:x", "NONE", MPI_COMM_WORLD);
  const std::string surfacePartName = "surface_1";
  stk::mesh::Part* surfacePart = stkMeshCreator.get_meta_data().get_part(surfacePartName);
  EXPECT_TRUE((surfacePart != nullptr));
  stk::mesh::Selector ownedSurface = stkMeshCreator.get_meta_data().locally_owned_part() & *surfacePart;
  const unsigned stkNumFaces = stk::mesh::count_entities(stkMeshCreator.get_bulk_data(), stk::topology::FACE_RANK, ownedSurface);
  stk_interface::MeshPart meshAndField = stkMeshCreator.create_mesh_from_part(surfacePartName);
  std::shared_ptr<mesh::Mesh> surfaceMesh = meshAndField.mesh;
  EXPECT_EQ(stkNumFaces, surfaceMesh->get_elements().size());

  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1) {
    verify_correct_sharing(stkMeshCreator.get_bulk_data(), *surfaceMesh, meshAndField.stkEls);
  }

  mesh::check_topology(surfaceMesh);
}

TEST(StkMeshCreator, create_mesh_from_shell_part)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 2) { GTEST_SKIP(); }

  stk_interface::StkMeshCreator stkMeshCreator("generated:1x2x2|shell:x", "NONE", MPI_COMM_WORLD);
  const std::string surfacePartName = "block_2";
  stk::mesh::Part* surfacePart = stkMeshCreator.get_meta_data().get_part(surfacePartName);
  EXPECT_TRUE((surfacePart != nullptr));
  stk::mesh::Selector ownedSurface = stkMeshCreator.get_meta_data().locally_owned_part() & *surfacePart;
  const unsigned stkNumFaces = stk::mesh::count_entities(stkMeshCreator.get_bulk_data(), stk::topology::ELEM_RANK, ownedSurface);
  stk_interface::MeshPart meshAndField = stkMeshCreator.create_mesh_from_part(surfacePartName);
  std::shared_ptr<mesh::Mesh> surfaceMesh = meshAndField.mesh;
  EXPECT_EQ(stkNumFaces, surfaceMesh->get_elements().size());

  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1) {
    verify_correct_sharing(stkMeshCreator.get_bulk_data(), *surfaceMesh, meshAndField.stkEls);
  }

  mesh::check_topology(surfaceMesh);
}

TEST(StkMeshCreator, create_mesh_from_bulk_data)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 2) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  read_stk_mesh("generated:1x2x2|sideset:x", *bulkDataPtr);

  stk_interface::StkMeshCreator stkMeshCreator(bulkDataPtr);

  const std::string surfacePartName = "surface_1";
  stk::mesh::Part* surfacePart = stkMeshCreator.get_meta_data().get_part(surfacePartName);
  EXPECT_TRUE((surfacePart != nullptr));
  stk::mesh::Selector ownedSurface = stkMeshCreator.get_meta_data().locally_owned_part() & *surfacePart;
  const unsigned stkNumFaces = stk::mesh::count_entities(stkMeshCreator.get_bulk_data(), stk::topology::FACE_RANK, ownedSurface);
  stk_interface::MeshPart meshAndField = stkMeshCreator.create_mesh_from_part(surfacePartName);
  std::shared_ptr<mesh::Mesh> surfaceMesh = meshAndField.mesh;
  EXPECT_EQ(stkNumFaces, surfaceMesh->get_elements().size());

  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1) {
    verify_correct_sharing(stkMeshCreator.get_bulk_data(), *surfaceMesh, meshAndField.stkEls);
  }

  mesh::check_topology(surfaceMesh);
}

}
}
}

