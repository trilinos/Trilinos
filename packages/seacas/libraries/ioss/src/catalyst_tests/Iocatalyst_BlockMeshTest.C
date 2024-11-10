// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "gtest/gtest.h"
#include <catalyst_tests/Iocatalyst_BlockMesh.h>

class BlockMeshTest : public ::testing::Test
{
protected:
  Iocatalyst::BlockMesh            bmOne;
  Iocatalyst::BlockMesh            bmTwo;
  Iocatalyst::BlockMesh            bmThree;
  Iocatalyst::BlockMesh            bmFour;
  Iocatalyst::BlockMesh::Partition part;
  Iocatalyst::BlockMesh::Extent    numBlocks;
  Iocatalyst::BlockMesh::Extent    origin;

  BlockMeshTest()
  {
    origin.i = 0;
    origin.j = 0;
    origin.k = 0;
  }

  void checkLocalNumBlocks(const Iocatalyst::BlockMesh &bm, int i, int j, int k)
  {
    EXPECT_EQ(bm.getPartitionExtents().i, i);
    EXPECT_EQ(bm.getPartitionExtents().j, j);
    EXPECT_EQ(bm.getPartitionExtents().k, k);
  }

  void checkLocalBlockStart(const Iocatalyst::BlockMesh &bm, int i, int j, int k)
  {
    EXPECT_EQ(bm.getPartitionStart().i, i);
    EXPECT_EQ(bm.getPartitionStart().j, j);
    EXPECT_EQ(bm.getPartitionStart().k, k);
  }
};

TEST_F(BlockMeshTest, GetID)
{
  EXPECT_EQ(bmOne.getID(), 1);
  EXPECT_EQ(bmTwo.getID(), 2);
  EXPECT_EQ(bmThree.getID(), 3);
  EXPECT_EQ(bmFour.getID(), 4);
  Iocatalyst::BlockMesh bm = bmOne;
  EXPECT_EQ(bm.getID(), 1);
}

TEST_F(BlockMeshTest, Defaults)
{

  EXPECT_EQ(bmOne.getPartition().id, 0);
  EXPECT_EQ(bmOne.getPartition().size, 1);

  EXPECT_EQ(bmOne.getExtents().i, 1);
  EXPECT_EQ(bmOne.getExtents().j, 1);
  EXPECT_EQ(bmOne.getExtents().k, 1);

  EXPECT_EQ(bmOne.getPartitionExtents().i, 1);
  EXPECT_EQ(bmOne.getPartitionExtents().j, 1);
  EXPECT_EQ(bmOne.getPartitionExtents().k, 1);

  EXPECT_EQ(bmOne.getPartitionStart().i, 0);
  EXPECT_EQ(bmOne.getPartitionStart().j, 0);
  EXPECT_EQ(bmOne.getPartitionStart().k, 0);

  EXPECT_FALSE(bmOne.isPartitionEmpty());
}

TEST_F(BlockMeshTest, InitPartitionSizeOne)
{

  numBlocks.i = 12;
  numBlocks.j = 45;
  numBlocks.k = 176;
  part.id     = 0;
  part.size   = 1;

  bmOne.init(part, numBlocks, origin);
  EXPECT_EQ(bmOne.getPartition().id, 0);
  EXPECT_EQ(bmOne.getPartition().size, 1);
  EXPECT_EQ(bmOne.getExtents().i, 12);
  EXPECT_EQ(bmOne.getExtents().j, 45);
  EXPECT_EQ(bmOne.getExtents().k, 176);
  checkLocalNumBlocks(bmOne, 12, 45, 176);
  checkLocalBlockStart(bmOne, 0, 0, 0);
}

TEST_F(BlockMeshTest, InitPartitionSizeTwoSmallestGrid)
{

  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 1;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmOne, 1, 1, 1);
  checkLocalBlockStart(bmOne, 0, 0, 0);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  EXPECT_TRUE(bmTwo.isPartitionEmpty());
}

TEST_F(BlockMeshTest, InitPartitionSizeTwoEvenExtents)
{

  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 12;
  numBlocks.j = 6;
  numBlocks.k = 24;
  bmOne.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmOne, 12, 6, 12);
  checkLocalBlockStart(bmOne, 0, 0, 0);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmTwo, 12, 6, 12);
  checkLocalBlockStart(bmTwo, 0, 0, 12);
}

TEST_F(BlockMeshTest, InitPartitionSizeTwoOddExtents)
{

  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 12;
  numBlocks.j = 6;
  numBlocks.k = 27;
  bmOne.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmOne, 12, 6, 13);
  checkLocalBlockStart(bmOne, 0, 0, 0);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmTwo, 12, 6, 14);
  checkLocalBlockStart(bmTwo, 0, 0, 13);
}

TEST_F(BlockMeshTest, InitPartitionSizeTwoZandYnumBlocksOne)
{
  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 13;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmOne, 6, 1, 1);
  checkLocalBlockStart(bmOne, 0, 0, 0);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmTwo, 7, 1, 1);
  checkLocalBlockStart(bmTwo, 6, 0, 0);
}

TEST_F(BlockMeshTest, InitPartitionSizeThree)
{
  part.id     = 0;
  part.size   = 3;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  EXPECT_TRUE(bmOne.isPartitionEmpty());
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmTwo, 2, 1, 1);
  checkLocalBlockStart(bmTwo, 0, 0, 0);
  part.id = 2;
  bmThree.init(part, numBlocks, origin);
  checkLocalNumBlocks(bmThree, 2, 1, 1);
  checkLocalBlockStart(bmThree, 0, 1, 0);
}

TEST_F(BlockMeshTest, GetPartitionPointIDsSizeOneOneBlock)
{
  std::vector<int> points = {1, 2, 3, 4, 5, 6, 7, 8};
  EXPECT_EQ(bmOne.getPartitionPointIDs(), points);
}

TEST_F(BlockMeshTest, GetPartitionPointIDsSizeOneTwoBlocksInX)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> points = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  EXPECT_EQ(bmOne.getPartitionPointIDs(), points);
}

TEST_F(BlockMeshTest, GetPartitionPointIDsSizeTwoOneBlock)
{

  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 1;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> points = {1, 2, 3, 4, 5, 6, 7, 8};
  EXPECT_EQ(bmOne.getPartitionPointIDs(), points);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  points.clear();
  EXPECT_EQ(bmTwo.getPartitionPointIDs(), points);
}

TEST_F(BlockMeshTest, GetPartitionPointIDsSizeTwoTwoBlocksInY)
{
  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 1;
  numBlocks.j = 2;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> points = {1, 2, 3, 4, 7, 8, 9, 10};
  EXPECT_EQ(bmOne.getPartitionPointIDs(), points);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  std::vector<int> pointsTwo = {3, 4, 5, 6, 9, 10, 11, 12};
  EXPECT_EQ(bmTwo.getPartitionPointIDs(), pointsTwo);
}

TEST_F(BlockMeshTest, GetPartitionPointIDsSizeFourEightBlocks)
{
  part.id     = 0;
  part.size   = 4;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> points = {1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15};
  EXPECT_EQ(bmOne.getPartitionPointIDs(), points);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  std::vector<int> pointsTwo = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18};
  EXPECT_EQ(bmTwo.getPartitionPointIDs(), pointsTwo);
  part.id = 2;
  bmThree.init(part, numBlocks, origin);
  std::vector<int> pointsThree = {10, 11, 12, 13, 14, 15, 19, 20, 21, 22, 23, 24};
  EXPECT_EQ(bmThree.getPartitionPointIDs(), pointsThree);
  part.id = 3;
  bmFour.init(part, numBlocks, origin);
  std::vector<int> pointsFour = {13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27};
  EXPECT_EQ(bmFour.getPartitionPointIDs(), pointsFour);
}

TEST_F(BlockMeshTest, GetPointCoordsForGlobalIDSizeOneOneBlock)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 1;
  numBlocks.j = 1;
  numBlocks.k = 1;
  origin      = {0, 0, 0};
  bmOne.init(part, numBlocks, origin);
  Iocatalyst::BlockMesh::Point p = bmOne.getPointCoordsForPointID(8);
  EXPECT_DOUBLE_EQ(p.x, 1.0);
  EXPECT_DOUBLE_EQ(p.y, 1.0);
  EXPECT_DOUBLE_EQ(p.z, 1.0);
}

TEST_F(BlockMeshTest, GetPointCoordsForGlobalIDSizeOneEightBlocks)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  origin      = {2, 2, 2};
  bmOne.init(part, numBlocks, origin);
  Iocatalyst::BlockMesh::Point p = bmOne.getPointCoordsForPointID(27);
  EXPECT_DOUBLE_EQ(p.x, 4.0);
  EXPECT_DOUBLE_EQ(p.y, 4.0);
  EXPECT_DOUBLE_EQ(p.z, 4.0);
  p = bmOne.getPointCoordsForPointID(1);
  EXPECT_DOUBLE_EQ(p.x, 2.0);
  EXPECT_DOUBLE_EQ(p.y, 2.0);
  EXPECT_DOUBLE_EQ(p.z, 2.0);
  p = bmOne.getPointCoordsForPointID(10);
  EXPECT_DOUBLE_EQ(p.x, 2.0);
  EXPECT_DOUBLE_EQ(p.y, 2.0);
  EXPECT_DOUBLE_EQ(p.z, 3.0);
}

TEST_F(BlockMeshTest, GetPartitionBlockIDsSizeOneOneBlock)
{
  std::vector<int> points = {1};
  EXPECT_EQ(bmOne.getPartitionBlockIDs(), points);
}

TEST_F(BlockMeshTest, GetPartitionBlockIDsSizeOneTwoBlocksInX)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> ids = {1, 2};
  EXPECT_EQ(bmOne.getPartitionBlockIDs(), ids);
}

TEST_F(BlockMeshTest, GetPartitionBlockIDsSizeTwoOneBlock)
{

  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 1;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> ids = {1};
  EXPECT_EQ(bmOne.getPartitionBlockIDs(), ids);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  ids.clear();
  EXPECT_EQ(bmTwo.getPartitionBlockIDs(), ids);
}

TEST_F(BlockMeshTest, GetPartitionBlockIDsSizeTwoTwoBlocksInY)
{
  part.id     = 0;
  part.size   = 2;
  numBlocks.i = 1;
  numBlocks.j = 2;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> ids = {1};
  EXPECT_EQ(bmOne.getPartitionBlockIDs(), ids);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  std::vector<int> idsTwo = {2};
  EXPECT_EQ(bmTwo.getPartitionBlockIDs(), idsTwo);
}

TEST_F(BlockMeshTest, GetPartitionBlockIDsSizeFourEightBlocks)
{
  part.id     = 0;
  part.size   = 4;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  std::vector<int> ids = {1, 2};
  EXPECT_EQ(bmOne.getPartitionBlockIDs(), ids);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  std::vector<int> idsTwo = {3, 4};
  EXPECT_EQ(bmTwo.getPartitionBlockIDs(), idsTwo);
  part.id = 2;
  bmThree.init(part, numBlocks, origin);
  std::vector<int> idsThree = {5, 6};
  EXPECT_EQ(bmThree.getPartitionBlockIDs(), idsThree);
  part.id = 3;
  bmFour.init(part, numBlocks, origin);
  std::vector<int> idsFour = {7, 8};
  EXPECT_EQ(bmFour.getPartitionBlockIDs(), idsFour);
}

TEST_F(BlockMeshTest, GetBlockConnectivityPointIDsSizeOneOneBlock)
{
  Iocatalyst::BlockMesh::BlockConn ids = {1, 2, 4, 3, 5, 6, 8, 7};
  EXPECT_EQ(bmOne.getBlockConnectivityPointIDs(1), ids);
}

TEST_F(BlockMeshTest, GetBlockConnectivityPointIDsSizeOneTwoBlocksInX)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 1;
  numBlocks.k = 1;
  bmOne.init(part, numBlocks, origin);
  Iocatalyst::BlockMesh::BlockConn ids = {2, 3, 6, 5, 8, 9, 12, 11};
  EXPECT_EQ(bmOne.getBlockConnectivityPointIDs(2), ids);
}

TEST_F(BlockMeshTest, GetBlockConnectivityPointIDsSizeOneEightBlocks)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  Iocatalyst::BlockMesh::BlockConn ids = {14, 15, 18, 17, 23, 24, 27, 26};
  EXPECT_EQ(bmOne.getBlockConnectivityPointIDs(8), ids);
  ids = {1, 2, 5, 4, 10, 11, 14, 13};
  EXPECT_EQ(bmOne.getBlockConnectivityPointIDs(1), ids);
  ids = {2, 3, 6, 5, 11, 12, 15, 14};
  EXPECT_EQ(bmOne.getBlockConnectivityPointIDs(2), ids);
}

TEST_F(BlockMeshTest, GetNumPointsSizeFourEightBlocks)
{
  part.id     = 0;
  part.size   = 4;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  EXPECT_EQ(bmOne.getNumPoints(), 27);
  EXPECT_EQ(bmOne.getNumPartitionPoints(), 12);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  EXPECT_EQ(bmTwo.getNumPoints(), 27);
  EXPECT_EQ(bmTwo.getNumPartitionPoints(), 12);
  part.id = 2;
  bmThree.init(part, numBlocks, origin);
  EXPECT_EQ(bmThree.getNumPoints(), 27);
  EXPECT_EQ(bmThree.getNumPartitionPoints(), 12);
  part.id = 3;
  bmFour.init(part, numBlocks, origin);
  EXPECT_EQ(bmFour.getNumPoints(), 27);
  EXPECT_EQ(bmFour.getNumPartitionPoints(), 12);
}

TEST_F(BlockMeshTest, GetNumBlocksSizeFourEightBlocks)
{
  part.id     = 0;
  part.size   = 4;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  EXPECT_EQ(bmOne.getNumBlocks(), 8);
  EXPECT_EQ(bmOne.getNumPartitionBlocks(), 2);
  part.id = 1;
  bmTwo.init(part, numBlocks, origin);
  EXPECT_EQ(bmTwo.getNumBlocks(), 8);
  EXPECT_EQ(bmTwo.getNumPartitionBlocks(), 2);
  part.id = 2;
  bmThree.init(part, numBlocks, origin);
  EXPECT_EQ(bmThree.getNumBlocks(), 8);
  EXPECT_EQ(bmThree.getNumPartitionBlocks(), 2);
  part.id = 3;
  bmFour.init(part, numBlocks, origin);
  EXPECT_EQ(bmFour.getNumBlocks(), 8);
  EXPECT_EQ(bmFour.getNumPartitionBlocks(), 2);
}

TEST_F(BlockMeshTest, GetGlobalIDForBlockID)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  Iocatalyst::BlockMesh::Extent coords = {0, 0, 0};
  Iocatalyst::BlockMesh::ID     id = bmOne.getIDfromCoords(coords, bmOne.getGlobalBlockExtents());
  EXPECT_EQ(bmOne.getGlobalIDForBlockID(1), id);
  coords = {1, 1, 1};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalBlockExtents());
  EXPECT_EQ(bmOne.getGlobalIDForBlockID(8), id);

  origin.i = 151;
  origin.j = 56;
  origin.k = 667;
  bmOne.init(part, numBlocks, origin);
  coords = {origin.i, origin.j, origin.k};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalBlockExtents());
  EXPECT_EQ(bmOne.getGlobalIDForBlockID(1), id);
  coords = {origin.i + 1, origin.j, origin.k};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalBlockExtents());
  EXPECT_EQ(bmOne.getGlobalIDForBlockID(2), id);
  coords = {origin.i + 1, origin.j + 1, origin.k + 1};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalBlockExtents());
  EXPECT_EQ(bmOne.getGlobalIDForBlockID(8), id);
}

TEST_F(BlockMeshTest, GetGlobalIDForPointID)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);
  Iocatalyst::BlockMesh::Extent coords = {0, 0, 0};
  Iocatalyst::BlockMesh::ID     id = bmOne.getIDfromCoords(coords, bmOne.getGlobalPointExtents());
  EXPECT_EQ(bmOne.getGlobalIDForPointID(1), id);
  coords = {2, 2, 2};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalPointExtents());
  EXPECT_EQ(bmOne.getGlobalIDForPointID(27), id);

  origin.i = 34;
  origin.j = 898;
  origin.k = 454;
  bmOne.init(part, numBlocks, origin);
  coords = {origin.i, origin.j, origin.k};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalPointExtents());
  EXPECT_EQ(bmOne.getGlobalIDForPointID(1), id);
  coords = {origin.i + 1, origin.j, origin.k};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalPointExtents());
  EXPECT_EQ(bmOne.getGlobalIDForPointID(2), id);
  coords = {origin.i + 2, origin.j + 2, origin.k + 2};
  id     = bmOne.getIDfromCoords(coords, bmOne.getGlobalPointExtents());
  EXPECT_EQ(bmOne.getGlobalIDForPointID(27), id);
}

TEST_F(BlockMeshTest, AddTransientField)
{
  part.id     = 0;
  part.size   = 1;
  numBlocks.i = 2;
  numBlocks.j = 2;
  numBlocks.k = 2;
  bmOne.init(part, numBlocks, origin);

  bmOne.addTransientCellField("foo_x", 2);
  bmOne.addTransientPointField("bar_x", 3);

  auto cell_fields = bmOne.getTransientCellFieldMap();
  EXPECT_EQ((*cell_fields)["foo_x"], 2);
  auto point_fields = bmOne.getTransientPointFieldMap();
  EXPECT_EQ((*point_fields)["bar_x"], 3);
}
