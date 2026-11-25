// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MeshBlockSkinSurface.hpp>
#include <Akri_ScopedTestFile.hpp>
#include <Akri_Unit_MultiBlockMesh.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <string>
#include <vector>

namespace krino
{

namespace
{

// The two hexes of two_hex_blocks, with the quadrilateral between them not planar, in a file written by
// processor 0 and removed at the end.
class TwoHexBlocksFile
{
public:
  static constexpr double sharedNodeShift = 0.2;

  TwoHexBlocksFile()
  {
    write_mesh_file_on_processor_0(myFile.name(),
        [](const stk::ParallelMachine comm) { return two_hex_blocks(sharedNodeShift, comm); });
  }

  const std::string & filename() const { return myFile.name(); }

private:
  ScopedTestFile myFile{"MeshBlockSkinSurfaceTwoHexBlocks.e"};
};

std::unique_ptr<MeshBlockSkinSurface> build_prepared_block_skin(std::shared_ptr<BlockSkinMeshFile> file,
    const std::vector<std::string> & blockNames,
    const int sign)
{
  auto skin = std::make_unique<MeshBlockSkinSurface>(3, file, blockNames, sign);
  skin->prepare_to_compute(0., skinQueryBox, 0.);
  return skin;
}

size_t global_count_of_owned_facets(const MeshFacetedSurface & surface)
{
  const size_t local = surface.num_owned_facets();
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

std::vector<stk::math::Vector3d> points_in_and_around_two_hexes()
{
  std::vector<stk::math::Vector3d> points;
  for (double x = -0.45; x < 2.5; x += 0.3)
    for (double y = -0.45; y < 1.5; y += 0.3)
      for (double z = -0.45; z < 1.5; z += 0.3)
        points.push_back(stk::math::Vector3d(x, y, z));
  return points;
}

} // namespace

// Reading the file gives the skins of the mesh it was written from: of each block and of both, with any
// sign, the same number of facets and vertices and the same signed distances.
TEST(MeshBlockSkinSurface, twoHexBlocksFile_skinsAreThoseOfTheMeshItWasWrittenFrom)
{
  const TwoHexBlocksFile file;
  const MultiBlockMesh mesh = two_hex_blocks(TwoHexBlocksFile::sharedNodeShift);
  for (const std::vector<unsigned> & blockIndices : std::vector<std::vector<unsigned>>{{0}, {1}, {0, 1}})
    for (const int sign : {1, -1})
    {
      std::vector<std::string> blockNames;
      for (auto b : blockIndices)
        blockNames.push_back(mesh.myBlocks[b]->name());
      const auto fileSkin = build_prepared_block_skin(get_block_skin_mesh_file(MPI_COMM_WORLD, file.filename()), blockNames, sign);
      const auto memorySkin = mesh.build_skin(blockIndices, sign);

      EXPECT_EQ(global_count_of_owned_facets(*memorySkin), global_count_of_owned_facets(*fileSkin));
      EXPECT_EQ(memorySkin->global_num_vertices(), fileSkin->global_num_vertices());
      for (auto & x : points_in_and_around_two_hexes())
      {
        EXPECT_DOUBLE_EQ(memorySkin->point_signed_distance(x), fileSkin->point_signed_distance(x)) << "at " << x;
      }
    }
}

// The bounding box is that of the blocks' nodes, before the facets are built and after the file's mesh is
// let go.
TEST(MeshBlockSkinSurface, twoHexBlocksFile_boundingBoxIsThatOfTheBlocksNodes)
{
  const TwoHexBlocksFile file;
  MeshBlockSkinSurface skin(3, get_block_skin_mesh_file(MPI_COMM_WORLD, file.filename()), {"block_2"}, 1);
  const BoundingBox expected(stk::math::Vector3d(1, 0, 0), stk::math::Vector3d(2, 1, 1));
  for (const bool isPrepared : {false, true})
  {
    if (isPrepared)
      skin.prepare_to_compute(0., skinQueryBox, 0.);
    BoundingBox bbox;
    skin.insert_into(bbox);
    EXPECT_EQ(expected.get_min(), bbox.get_min()) << "prepared " << isPrepared;
    EXPECT_EQ(expected.get_max(), bbox.get_max()) << "prepared " << isPrepared;
  }
}

// The surfaces of one file share its mesh, read once, which is let go when the last of them has its facets.
TEST(MeshBlockSkinSurface, surfacesOfOneFile_shareItsMeshAndLetItGoOnceBuilt)
{
  const TwoHexBlocksFile file;
  auto sharedFile = get_block_skin_mesh_file(MPI_COMM_WORLD, file.filename());
  EXPECT_EQ(sharedFile, get_block_skin_mesh_file(MPI_COMM_WORLD, file.filename()));
  MeshBlockSkinSurface skin1(3, sharedFile, {"block_1"}, 1);
  MeshBlockSkinSurface skin2(3, sharedFile, {"block_2"}, 1);
  const std::weak_ptr<BlockSkinMeshFile> heldFile = sharedFile;
  sharedFile.reset();

  skin1.prepare_to_compute(0., skinQueryBox, 0.);
  EXPECT_FALSE(heldFile.expired());
  skin2.prepare_to_compute(0., skinQueryBox, 0.);
  EXPECT_TRUE(heldFile.expired());
}

// A name that is not a block of the file is an error that lists the file's blocks.
TEST(MeshBlockSkinSurface, nameThatIsNotABlock_throwsListingTheBlocks)
{
  const TwoHexBlocksFile file;
  MeshBlockSkinSurface skin(3, get_block_skin_mesh_file(MPI_COMM_WORLD, file.filename()), {"block_7"}, 1);
  try
  {
    skin.prepare_to_compute(0., skinQueryBox, 0.);
    ADD_FAILURE() << "expected an error for block_7";
  }
  catch (const std::exception & error)
  {
    const std::string message = error.what();
    EXPECT_NE(std::string::npos, message.find("no block named block_7")) << message;
    EXPECT_NE(std::string::npos, message.find("block_1 block_2")) << message;
  }
}

} // namespace krino
