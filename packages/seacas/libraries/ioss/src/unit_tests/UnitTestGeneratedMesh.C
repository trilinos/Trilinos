// Copyright(C) 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "generated/Iogn_GeneratedMesh.h"
#include "gtest/gtest.h"
#include <stdexcept>

TEST(GeneratedMesh, NoSideSets)
{
  Iogn::GeneratedMesh gm("10x12x14");
  EXPECT_THROW(gm.sideset_touching_blocks(0), std::runtime_error);
}

TEST(GeneratedMesh, TwoSideSets)
{
  Iogn::GeneratedMesh gm("10x12x14|sideset:XY");
  EXPECT_EQ(gm.sideset_count(), 2);
  EXPECT_EQ(gm.sideset_touching_blocks(1)[0], "block_1");
  EXPECT_EQ(gm.sideset_touching_blocks(2)[0], "block_1");
}

TEST(GeneratedMesh, OneShellOneSideSetDifferentSides)
{
  Iogn::GeneratedMesh gm("10x12x14|shell:z|sideset:Z");
  EXPECT_EQ(gm.sideset_count(), 1);
  EXPECT_EQ(gm.block_count(), 2);
  EXPECT_EQ(gm.sideset_touching_blocks(1)[0], "block_1");
}

TEST(GeneratedMesh, OneShellOneSideSetSameSide)
{
  Iogn::GeneratedMesh gm("10x12x14|shell:z|sideset:z");
  EXPECT_EQ(gm.sideset_count(), 1);
  EXPECT_EQ(gm.block_count(), 2);
  EXPECT_EQ(gm.sideset_touching_blocks(1)[0], "block_2");
}

TEST(GeneratedMesh, FourShellsOppositeSides)
{
  Iogn::GeneratedMesh gm("10x12x14|shell:yyYY|sideset:Y");
  EXPECT_EQ(gm.sideset_count(), 1);
  EXPECT_EQ(gm.block_count(), 5);
  EXPECT_EQ(gm.sideset_touching_blocks(1)[0], "block_5");
}

TEST(GeneratedMesh, ShellIdsatLocationPXEmpty)
{
  Iogn::GeneratedMesh gm("10x12x14");
  auto                shell_ids = gm.shell_ids_at_location(Iogn::GeneratedMesh::ShellLocation::PX);
  EXPECT_TRUE(shell_ids.empty());
}

TEST(GeneratedMesh, ShellIdsatLocationMX)
{
  Iogn::GeneratedMesh gm("10x12x14|shell:xxx");
  auto                shell_ids = gm.shell_ids_at_location(Iogn::GeneratedMesh::ShellLocation::MX);
  EXPECT_EQ(shell_ids.size(), 3);
  EXPECT_EQ(shell_ids.back(), 2);
  EXPECT_EQ(shell_ids.front(), 4);
  EXPECT_EQ(shell_ids[1], 3);
}

TEST(GeneratedMesh, ShellIdsatLocationPX)
{
  Iogn::GeneratedMesh gm("10x12x14|shell:xxxXX");
  auto                shell_ids = gm.shell_ids_at_location(Iogn::GeneratedMesh::ShellLocation::PX);
  EXPECT_EQ(shell_ids.size(), 2);
  EXPECT_EQ(shell_ids.back(), 6);
  EXPECT_EQ(shell_ids.front(), 5);
}

TEST(GeneratedMesh, ShellIdsatLocationMZ)
{
  Iogn::GeneratedMesh gm("10x12x14|shell:xxxXXyyyYYzZZ");
  auto                shell_ids = gm.shell_ids_at_location(Iogn::GeneratedMesh::ShellLocation::MZ);
  EXPECT_EQ(shell_ids.size(), 1);
  EXPECT_EQ(shell_ids.front(), 12);
}
