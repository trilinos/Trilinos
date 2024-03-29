// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>

TEST_F(Iocatalyst_DatabaseIOTest, WriteThreeElementBlocksWith24Cells)
{
  Iocatalyst::BlockMesh bmOne;
  setBlockMeshSize(2, 2, 2);
  addBlockMesh(bmOne);

  Iocatalyst::BlockMesh bmTwo;
  setOrigin(2, 0, 0);
  setBlockMeshSize(2, 2, 2);
  addBlockMesh(bmTwo);

  Iocatalyst::BlockMesh bmThree;
  setOrigin(4, 0, 0);
  setBlockMeshSize(2, 2, 2);
  addBlockMesh(bmThree);

  runUnstructuredTest("test_eb_3_cells_24");
}

TEST_F(Iocatalyst_DatabaseIOTest, WriteOneElementBlockWith8Cells)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);
  addBlockMesh(bm);
  runUnstructuredTest("test_eb_1_cells_8");
}

TEST_F(Iocatalyst_DatabaseIOTest, WriteOneElementBlockWith300Cells)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(10, 10, 3);
  addBlockMesh(bm);
  runUnstructuredTest("test_eb_1_cells_300");
}

TEST_F(Iocatalyst_DatabaseIOTest, WriteThreeElementBlocksWith835Cells)
{
  Iocatalyst::BlockMesh bmOne;
  setBlockMeshSize(5, 15, 3);
  addBlockMesh(bmOne);

  Iocatalyst::BlockMesh bmTwo;
  setOrigin(5, 0, 0);
  setBlockMeshSize(8, 14, 4);
  addBlockMesh(bmTwo);

  Iocatalyst::BlockMesh bmThree;
  setOrigin(13, 0, 0);
  setBlockMeshSize(3, 6, 9);
  addBlockMesh(bmThree);

  runUnstructuredTest("test_eb_3_cells_835");
}