// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>

TEST_F(Iocatalyst_DatabaseIOTest, GetNumLocalPointsInMeshSet)
{
  EXPECT_EQ(bmSet.getNumLocalPointsInMeshSet(), 0);
  Iocatalyst::BlockMesh bmOne;
  setBlockMeshSize(2, 2, 2);
  addBlockMesh(bmOne);
  EXPECT_EQ(bmSet.getNumLocalPointsInMeshSet(), 27);

  Iocatalyst::BlockMesh bmTwo;
  setBlockMeshSize(1, 1, 1);
  setOrigin(2, 0, 0);
  addBlockMesh(bmTwo);
  EXPECT_EQ(bmSet.getNumLocalPointsInMeshSet(), 31);

  Iocatalyst::BlockMesh bmThree;
  setBlockMeshSize(2, 2, 2);
  setOrigin(3, 0, 0);
  addBlockMesh(bmThree);
  EXPECT_EQ(bmSet.getNumLocalPointsInMeshSet(), 54);

  Iocatalyst::BlockMesh bmFour;
  setBlockMeshSize(8, 3, 3);
  setOrigin(5, 0, 0);
  addBlockMesh(bmFour);
  EXPECT_EQ(bmSet.getNumLocalPointsInMeshSet(), 189);
}

TEST_F(Iocatalyst_DatabaseIOTest, AddTransientFieldToBlockMesh)
{
  Iocatalyst::BlockMesh bmOne;
  setBlockMeshSize(2, 2, 2);

  bmOne.addTransientCellField("foo_x", 2);
  bmOne.addTransientPointField("bar_x", 3);

  addBlockMesh(bmOne);

  std::string exodusFileName =
      "AddTransientFieldToBlockMesh" + CATALYST_TEST_FILE_NP + std::to_string(part.size) + EXODUS_FILE_EXTENSION;
  Iocatalyst::BlockMeshSet::IOSSparams iop(exodusFileName, EXODUS_DATABASE_TYPE);
  bmSet.writeIOSSFile(iop);
}