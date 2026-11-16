#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>

namespace {

void read_mesh_and_check_elems(MPI_Comm comm, const std::string fileName,
                               unsigned expectedNumElems,
                               stk::topology expectedTopology)
{
  std::shared_ptr<stk::mesh::BulkData> readBulkPtr = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh(fileName, *readBulkPtr);
  const stk::mesh::MetaData& readMeta = readBulkPtr->mesh_meta_data();

  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(*readBulkPtr, stk::topology::ELEM_RANK, readMeta.universal_part(), elems);
  EXPECT_EQ(expectedNumElems, elems.size());
  for(stk::mesh::Entity elem : elems) {
    stk::topology elemTopo = readBulkPtr->bucket(elem).topology();
    EXPECT_EQ(expectedTopology, elemTopo) << "elem "<<readBulkPtr->entity_key(elem) << " has topology " << elemTopo << " but expected " << expectedTopology;
  }
}

void write_surface_mesh(const std::string& fileName, stk::mesh::BulkData& bulk,
                        unsigned numSurfaces)
{
  stk::io::StkMeshIoBroker ioBroker;
  ioBroker.set_bulk_data(bulk);

  ASSERT_TRUE(numSurfaces >= 1);

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::Part* surface1 = meta.get_part("surface_1");
  stk::mesh::Selector skinSelector = *surface1;

  for(unsigned surfIdx=1; surfIdx<numSurfaces; ++surfIdx) {
    std::string surfName = "surface_"+std::to_string(surfIdx+1);
    stk::mesh::Part* surface = meta.get_part(surfName);
    ASSERT_TRUE(surface != nullptr);
    skinSelector |= *surface;
  }

  size_t fileIdx = ioBroker.create_output_mesh(fileName, stk::io::WRITE_RESULTS);

  ioBroker.set_skin_mesh_selector(fileIdx, skinSelector);

  ioBroker.write_output_mesh(fileIdx);
}

TEST(OutputExposedSurface, DISABLED_Quad4_2D)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();

  std::string meshDesc = "textmesh: 0,1,QUAD_4_2D,1,2,5,4, block_1\n"
                                   "0,2,QUAD_4_2D,2,3,6,5, block_2\n"
                         "|dimension: 2"
                         "|coordinates: 0,0, 1,0, 2,0, 0,1, 1,1, 2,1"
                         "|sideset:name=surface_1; data=1,1, 1,2, 1,3, 1,4"
                         "|sideset:name=surface_2; data=2,1, 2,2, 2,3, 2,4";

  stk::io::fill_mesh(meshDesc, *bulkPtr);

  std::string fileName("quadSurf.exo");
  write_surface_mesh(fileName, *bulkPtr, 2);

  read_mesh_and_check_elems(comm, fileName, 6, stk::topology::BEAM_2);
  unlink(fileName.c_str());
}

TEST(OutputExposedSurface, Hex8_3D)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();

  std::string meshDesc = "textmesh: 0,1,HEX_8,1,2,3,4,5,6,7,8, block_1\n"
                                   "0,2,HEX_8,5,6,7,8,9,10,11,12, block_2\n"
                         "|coordinates: 0,0,0, 1,0,0, 1,1,0, 0,1,0,"
                                       "0,0,1, 1,0,1, 1,1,1, 0,1,1,"
                                       "0,0,2, 1,0,2, 1,1,2, 0,1,2"
                         "|sideset:name=surface_1; data=1,1, 1,2, 1,3, 1,4, 1,5, 1,6"
                         "|sideset:name=surface_2; data=2,1, 2,2, 2,3, 2,4, 2,5, 2,6";

  stk::io::fill_mesh(meshDesc, *bulkPtr);

  std::string fileName("hexSurf.exo");
  write_surface_mesh(fileName, *bulkPtr, 2);

  read_mesh_and_check_elems(comm, fileName, 10, stk::topology::SHELL_QUAD_4);
  unlink(fileName.c_str());
}

TEST(OutputExposedSurface, Tet4_3D)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();

  std::string meshDesc = "textmesh: 0,1,TET_4,1,2,3,4, block_1\n"
                                   "0,2,TET_4,2,5,3,4, block_2\n"
                         "|coordinates: 0,0,0, 1,0,0, 0,1,0, 1,1,0, 0.5,0.5,1"
                         "|sideset:name=surface_1; data=1,1, 1,2, 1,3, 1,4"
                         "|sideset:name=surface_2; data=2,1, 2,2, 2,3, 2,4";

  stk::io::fill_mesh(meshDesc, *bulkPtr);

  std::string fileName("tetSurf.exo");
  write_surface_mesh(fileName, *bulkPtr, 2);

  read_mesh_and_check_elems(comm, fileName, 6, stk::topology::SHELL_TRI_3);
  unlink(fileName.c_str());
}

TEST(OutputExposedSurface, Tet10)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();

  std::string meshDesc = "textmesh: 0,1,TET_10,1,2,3,4,5,6,7,8,9,10,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, "
         "0.5,0.0,0.0, 0.5,0.5,0.0, 0.0,0.5,0.0, 0.0,0.0,0.5, 0.5,0.0,0.5, 0.0,0.5,0.5"
                         "|sideset:name=surface_1; data=1,1, 1,2, 1,3, 1,4";

  stk::io::fill_mesh(meshDesc, *bulkPtr);

  std::string fileName("tet10Surf.exo");
  write_surface_mesh(fileName, *bulkPtr, 1);

  read_mesh_and_check_elems(comm, fileName, 4, stk::topology::SHELL_TRI_6);
}

TEST(OutputExposedSurface, Hex20)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();

  std::string meshDesc = "textmesh: 0,1,HEX_20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,block_1"
         "|coordinates: 0.0,0.0,0.0, 1.0,0.0,0.0, 1.0,1.0,0.0, 0.0,1.0,0.0, "
                       "0.0,0.0,1.0, 1.0,0.0,1.0, 1.0,1.0,1.0, 0.0,1.0,1.0, "
                       "0.5,0.0,0.0, 1.0,0.5,0.0, 0.5,1.0,0.0, 0.0,0.5,0.0, "
                       "0.0,0.0,0.5, 1.0,0.0,0.5, 1.0,1.0,0.5, 0.0,1.0,0.5, "
                       "0.5,0.0,1.0, 1.0,0.5,1.0, 0.5,1.0,1.0, 0.0,0.5,1.0  "
         "|sideset:name=surface_1; data=1,1, 1,2, 1,3, 1,4, 1,5, 1,6";

  stk::io::fill_mesh(meshDesc, *bulkPtr);

  std::string fileName("hex20Surf.exo");
  write_surface_mesh(fileName, *bulkPtr, 1);

  read_mesh_and_check_elems(comm, fileName, 6, stk::topology::SHELL_QUAD_8);
  unlink(fileName.c_str());
}

TEST(OutputExposedSurface, Hex8_with_shells)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();

  std::string meshDesc = "generated:1x1x1|shell:xX|sideset:xXyYzZ";

  stk::io::fill_mesh(meshDesc, *bulkPtr);

  std::string fileName("hex8ShellSurf.exo");
  write_surface_mesh(fileName, *bulkPtr, 6);

  read_mesh_and_check_elems(comm, fileName, 6, stk::topology::SHELL_QUAD_4);
  unlink(fileName.c_str());
}

}

