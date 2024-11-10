// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_ElementBlock.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_StructuredBlock.h>
#include <catalyst/Iocatalyst_DatabaseIO.h>
#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>
#include <cstdint>

TEST_F(Iocatalyst_DatabaseIOTest, ReadConduitCanExo)
{
  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ""));

  auto db = getCatalystDatabaseFromConduitFiles("Iocatalyst_can_ex2_MPI_1", iossProp);
  ASSERT_TRUE(db != nullptr);
  Ioss::Region reg(db);

  auto eb = reg.get_element_block("block_1");
  EXPECT_TRUE(eb->field_exists("eqps"));

  auto nb = reg.get_node_block("nodeblock_1");
  EXPECT_TRUE(nb->field_exists("accl"));
  EXPECT_TRUE(nb->field_exists("vel"));
  EXPECT_TRUE(nb->field_exists("displ"));
}

TEST_F(Iocatalyst_DatabaseIOTest, ReadConduitSPARCOneCGNS)
{
  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ""));

  auto db = getCatalystDatabaseFromConduitFiles("Iocatalyst_sparc1_cgns_MPI_1", iossProp);
  ASSERT_TRUE(db != nullptr);
  Ioss::Region reg(db);

  auto sb = reg.get_structured_block("blk-1");
  EXPECT_TRUE(sb != nullptr);
  EXPECT_TRUE(sb->field_exists("Density1"));
  EXPECT_TRUE(sb->field_exists("Temperature1"));
  EXPECT_TRUE(sb->field_exists("Velocity"));

  auto nb = sb->get_node_block();
  EXPECT_TRUE(nb.field_exists("mesh_model_coordinates"));
  EXPECT_TRUE(nb.field_exists("mesh_model_coordinates_x"));
  EXPECT_TRUE(nb.field_exists("mesh_model_coordinates_y"));
  EXPECT_TRUE(nb.field_exists("mesh_model_coordinates_z"));
}

TEST_F(Iocatalyst_DatabaseIOTest, SetReaderTimeStepWithIOSSProp)
{
  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("CATALYST_READER_TIME_STEP", 12));

  auto db = getCatalystDatabaseFromConduitFiles("Iocatalyst_can_ex2_MPI_1", iossProp);
  ASSERT_TRUE(db != nullptr);

  Iocatalyst::DatabaseIO *catdb = static_cast<Iocatalyst::DatabaseIO *>(db);

  Ioss::Region reg(db);

  auto mint = reg.get_min_time();
  EXPECT_EQ(mint.first, 1);
  EXPECT_DOUBLE_EQ(mint.second, 0.0011999331181868911);

  auto maxt = reg.get_max_time();
  EXPECT_EQ(maxt.first, 1);
  EXPECT_DOUBLE_EQ(maxt.second, 0.0011999331181868911);
}

TEST_F(Iocatalyst_DatabaseIOTest, SetReaderTimeStepWithIOSSEnvVar)
{
  setenv("CATALYST_READER_TIME_STEP", "24", 1);

  auto db = getCatalystDatabaseFromConduitFiles("Iocatalyst_can_ex2_MPI_1");
  ASSERT_TRUE(db != nullptr);

  Ioss::Region reg(db);

  auto mint = reg.get_min_time();
  EXPECT_EQ(mint.first, 1);
  EXPECT_DOUBLE_EQ(mint.second, 0.0024000538978725672);

  auto maxt = reg.get_max_time();
  EXPECT_EQ(maxt.first, 1);
  EXPECT_DOUBLE_EQ(maxt.second, 0.0024000538978725672);
}

TEST_F(Iocatalyst_DatabaseIOTest, SetRankNumRanksSerialParallel)
{
  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("my_processor", 0));
  iossProp.add(Ioss::Property("processor_count", 1));

  auto db = getCatalystDatabaseFromConduitFiles("Iocatalyst_can_ex2_MPI_1", iossProp);
  ASSERT_TRUE(db != nullptr);
}

TEST_F(Iocatalyst_DatabaseIOTest, CellIdsAndCellNodeIds)
{
  setenv("CATALYST_READER_TIME_STEP", "1", 1);

  auto db = getCatalystDatabaseFromConduitFiles("Iocatalyst_sparc1_cgns_MPI_1");
  ASSERT_TRUE(db != nullptr);

  Ioss::Region reg(db);

  auto sb = reg.get_structured_block("blk-1");

  EXPECT_TRUE(sb->field_exists("cell_ids"));
  EXPECT_TRUE(sb->field_exists("cell_node_ids"));

  auto cids = sb->get_fieldref("cell_ids");
  EXPECT_TRUE(cids.get_type() == Ioss::Field::INTEGER);
  std::vector<int32_t> cidBuff(cids.raw_count());
  EXPECT_EQ(sb->get_field_data("cell_ids", Data(cidBuff), sizeof(int32_t) * cidBuff.size()),
            cids.raw_count());
  EXPECT_EQ(cidBuff[0], 1);
  EXPECT_EQ(cidBuff[cids.raw_count() - 1], 256);


  auto cnids = sb->get_fieldref("cell_node_ids");
  EXPECT_TRUE(cnids.get_type() == Ioss::Field::INTEGER);
  std::vector<int32_t> cnidsBuff(cnids.raw_count());
  EXPECT_EQ(
      sb->get_field_data("cell_node_ids", Data(cnidsBuff), sizeof(int32_t) * cnidsBuff.size()),
      cnids.raw_count());
  EXPECT_EQ(cnidsBuff[0], 1);
  EXPECT_EQ(cnidsBuff[cnids.raw_count() - 1], 594);
}