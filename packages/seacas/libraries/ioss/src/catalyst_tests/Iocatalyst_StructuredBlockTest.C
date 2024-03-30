// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>

TEST_F(Iocatalyst_DatabaseIOTest, WriteOneStructuredBlockWith8Cells)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);
  addBlockMesh(bm);
  runStructuredTest("test_sb_1_cells_8");
}

TEST_F(Iocatalyst_DatabaseIOTest, WriteOneStructuredBlockWith200Cells)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(10, 10, 2);
  addBlockMesh(bm);
  runStructuredTest("test_sb_1_cells_200");
}

TEST_F(Iocatalyst_DatabaseIOTest, WriteThreeStructuredBlocksWith835Cells)
{
  Iocatalyst::BlockMesh bmOne;
  setBlockMeshSize(5, 15, 3);
  addBlockMesh(bmOne);

  Iocatalyst::BlockMesh bmTwo;
  setBlockMeshSize(8, 14, 4);
  setOrigin(5, 0, 0);
  addBlockMesh(bmTwo);

  Iocatalyst::BlockMesh bmThree;
  setBlockMeshSize(3, 6, 9);
  setOrigin(13, 0, 0);
  addBlockMesh(bmThree);
  runStructuredTest("test_sb_3_cells_835");
}