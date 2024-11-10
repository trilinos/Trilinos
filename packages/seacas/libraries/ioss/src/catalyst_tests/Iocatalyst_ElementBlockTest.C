// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_Compare.h>
#include <Ioss_CopyDatabase.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_IOFactory.h>
#include <Ioss_MeshCopyOptions.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_StructuredBlock.h>
#include <catalyst/Iocatalyst_DatabaseIO.h>
#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>
#include <catalyst/Iocatalyst_CatalystManager.h>



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

TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_ENABLE_FIELD_RECOGNITION_ON)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  bm.addTransientCellField("foo_x", 2);
  bm.addTransientCellField("foo_y", 3);
  bm.addTransientCellField("foo_z", 4);


  addBlockMesh(bm);

  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("ENABLE_FIELD_RECOGNITION", "ON"));
  Ioss::DatabaseIO *exo_d = writeAndGetExodusDatabaseOnRead("test_eb_1_enable_field_recog", iossProp);

  Iocatalyst::BlockMeshSet::IOSSparams iop("cat", EXODUS_DATABASE_TYPE, iossProp);

  Ioss::DatabaseIO *cat_d = bmSet.getCatalystDatabase(iop);

  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize on read"; }
  Ioss::Region cat_reg(cat_d);
  
  Ioss::Region exo_reg(exo_d);
  
  auto cat_elemBlock = cat_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));
  auto exo_elemBlock = exo_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));

  checkEntityContainerZeroCopyFields(cat_reg.get_element_blocks());

  bool exo_foo_exists = exo_elemBlock->field_exists("foo");
  bool cat_foo_exists = cat_elemBlock->field_exists("foo");
  EXPECT_TRUE(exo_foo_exists);
  EXPECT_TRUE(cat_foo_exists);
  if(exo_foo_exists && cat_foo_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("foo") == cat_elemBlock->get_field("foo"));
  
  //Check field data for equality
  auto cat_field = cat_elemBlock->get_fieldref("foo");
  std::vector<std::byte> dcBuffer(cat_field.get_size());
  cat_elemBlock->get_field_data("foo", Data(dcBuffer), dcBuffer.size());

  exo_reg.begin_state(1);
  auto exo_field = exo_elemBlock->get_fieldref("foo");
  std::vector<std::byte> deBuffer(exo_field.get_size());
  exo_elemBlock->get_field_data("foo", Data(deBuffer), deBuffer.size());
  EXPECT_EQ(dcBuffer, deBuffer);

  //Check foo_x doesn't exist
  bool exo_foo_x_exists = exo_elemBlock->field_exists("foo_x");
  bool cat_foo_x_exists = cat_elemBlock->field_exists("foo_x");
  EXPECT_FALSE(exo_foo_x_exists);
  EXPECT_FALSE(cat_foo_x_exists);
}

//Sanity Test
TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_ENABLE_FIELD_RECOGNITION_OFF)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  bm.addTransientCellField("foo_x", 2);
  bm.addTransientCellField("foo_y", 3);
  bm.addTransientCellField("foo_z", 4);


  addBlockMesh(bm);

  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("ENABLE_FIELD_RECOGNITION", "OFF"));
  Ioss::DatabaseIO *exo_d = writeAndGetExodusDatabaseOnRead("test_eb_1_enable_field_recog", iossProp);
  
  Iocatalyst::BlockMeshSet::IOSSparams iop("cat", EXODUS_DATABASE_TYPE, iossProp);

  Ioss::DatabaseIO *cat_d = bmSet.getCatalystDatabase(iop);

  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize on read"; }
  Ioss::Region cat_reg(cat_d);
  
  Ioss::Region exo_reg(exo_d);
  
  auto cat_elemBlock = cat_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));
  auto exo_elemBlock = exo_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));

  checkEntityContainerZeroCopyFields(cat_reg.get_element_blocks());

  bool exo_foo_x_exists = exo_elemBlock->field_exists("foo_x");
  bool cat_foo_x_exists = cat_elemBlock->field_exists("foo_x");
  EXPECT_TRUE(exo_foo_x_exists);
  EXPECT_TRUE(cat_foo_x_exists);
  if(exo_foo_x_exists && cat_foo_x_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("foo_x") == cat_elemBlock->get_field("foo_x"));
  
}


TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_IGNORE_REALN_FIELDS_ON)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  bm.addTransientCellField("foo_1", 2);
  bm.addTransientCellField("foo_2", 3);
  bm.addTransientCellField("foo_3", 4);


  addBlockMesh(bm);

  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("IGNORE_REALN_FIELDS", "ON"));
  Ioss::DatabaseIO *exo_d = writeAndGetExodusDatabaseOnRead("test_eb_1_ignore_realn_fields", iossProp);

  Iocatalyst::BlockMeshSet::IOSSparams iop("cat", EXODUS_DATABASE_TYPE, iossProp);

  Ioss::DatabaseIO *cat_d = bmSet.getCatalystDatabase(iop);

  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize on read"; }
  Ioss::Region cat_reg(cat_d);
  
  Ioss::Region exo_reg(exo_d);
  
  auto cat_elemBlock = cat_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));
  auto exo_elemBlock = exo_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));

  checkEntityContainerZeroCopyFields(cat_reg.get_element_blocks());

  bool exo_foo_1_exists = exo_elemBlock->field_exists("foo_1");
  bool cat_foo_1_exists = cat_elemBlock->field_exists("foo_1");
  EXPECT_TRUE(exo_foo_1_exists);
  EXPECT_TRUE(cat_foo_1_exists);
  if(exo_foo_1_exists && cat_foo_1_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("foo_1") == cat_elemBlock->get_field("foo_1"));
  
}

TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_IGNORE_REALN_FIELDS_OFF)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  bm.addTransientCellField("foo_1", 3);
  bm.addTransientCellField("foo_2", 2);
  bm.addTransientCellField("foo_3", 4);


  addBlockMesh(bm);

  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("IGNORE_REALN_FIELDS", "OFF"));
  Ioss::DatabaseIO *exo_d = writeAndGetExodusDatabaseOnRead("test_eb_1_ignore_realn_fields_off", iossProp);

  Iocatalyst::BlockMeshSet::IOSSparams iop("cat", EXODUS_DATABASE_TYPE, iossProp);

  Ioss::DatabaseIO *cat_d = bmSet.getCatalystDatabase(iop);

  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize on read"; }
  Ioss::Region cat_reg(cat_d);
  
  Ioss::Region exo_reg(exo_d);
  
  auto cat_elemBlock = cat_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));
  auto exo_elemBlock = exo_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));

  checkEntityContainerZeroCopyFields(cat_reg.get_element_blocks());

  bool exo_foo_exists = exo_elemBlock->field_exists("foo");
  bool cat_foo_exists = cat_elemBlock->field_exists("foo");
  EXPECT_TRUE(exo_foo_exists);
  EXPECT_TRUE(cat_foo_exists);
  if(exo_foo_exists && cat_foo_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("foo") == cat_elemBlock->get_field("foo"));
  
  //Check field data for equality
  auto cat_field = cat_elemBlock->get_fieldref("foo");
  std::vector<std::byte> dcBuffer(cat_field.get_size());
  cat_elemBlock->get_field_data("foo", Data(dcBuffer), dcBuffer.size());

  exo_reg.begin_state(1);
  auto exo_field = exo_elemBlock->get_fieldref("foo");
  std::vector<std::byte> deBuffer(exo_field.get_size());
  exo_elemBlock->get_field_data("foo", Data(deBuffer), deBuffer.size());
  EXPECT_EQ(dcBuffer, deBuffer);
  
}

TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_FIELD_SUFFIX_SEPARATOR)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  bm.addTransientCellField("foo_x", 2);
  bm.addTransientCellField("foo_y", 3);
  bm.addTransientCellField("foo_z", 4);
  bm.addTransientCellField("bar:x", 5);
  bm.addTransientCellField("bar:y", 6);
  bm.addTransientCellField("bar:z", 7);


  addBlockMesh(bm);

  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ":"));
  
  Ioss::DatabaseIO *exo_d = writeAndGetExodusDatabaseOnRead("test_eb_1_field_suf_sep", iossProp);

  Iocatalyst::BlockMeshSet::IOSSparams iop("cat", EXODUS_DATABASE_TYPE, iossProp);
  Ioss::DatabaseIO *cat_d = bmSet.getCatalystDatabase(iop);


  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize on read"; }
  Ioss::Region cat_reg(cat_d);
  
  Ioss::Region exo_reg(exo_d);
  
  auto cat_elemBlock = cat_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));
  auto exo_elemBlock = exo_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));

  checkEntityContainerZeroCopyFields(cat_reg.get_element_blocks());

  bool exo_foo_x_exists = exo_elemBlock->field_exists("foo_x");
  bool cat_foo_x_exists = cat_elemBlock->field_exists("foo_x");
  EXPECT_TRUE(exo_foo_x_exists);
  EXPECT_TRUE(cat_foo_x_exists);
  if(exo_foo_x_exists && cat_foo_x_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("foo_x") == cat_elemBlock->get_field("foo_x"));
  
  bool exo_bar_exists = exo_elemBlock->field_exists("bar");
  bool cat_bar_exists = cat_elemBlock->field_exists("bar");
  EXPECT_TRUE(exo_bar_exists);
  EXPECT_TRUE(cat_bar_exists);
  if(exo_bar_exists && cat_bar_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("bar") == cat_elemBlock->get_field("bar"));
  
  //Check bar field data for equality
  auto cat_field = cat_elemBlock->get_fieldref("bar");
  std::vector<std::byte> dcBuffer(cat_field.get_size());
  cat_elemBlock->get_field_data("bar", Data(dcBuffer), dcBuffer.size());

  exo_reg.begin_state(1);
  auto exo_field = exo_elemBlock->get_fieldref("bar");
  std::vector<std::byte> deBuffer(exo_field.get_size());
  exo_elemBlock->get_field_data("bar", Data(deBuffer), deBuffer.size());
  EXPECT_EQ(dcBuffer, deBuffer);

}

TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_FIELD_STRIP_TRAILING_UNDERSCORE)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  bm.addTransientCellField("foo_x", 2);
  bm.addTransientCellField("foo_y", 3);
  bm.addTransientCellField("foo_z", 4);

  addBlockMesh(bm);

  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("FIELD_STRIP_TRAILING_UNDERSCORE", "ON"));
  iossProp.add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ""));

  Ioss::DatabaseIO *exo_d = writeAndGetExodusDatabaseOnRead("test_eb_1_field_strip_tr_unders", iossProp);

  Iocatalyst::BlockMeshSet::IOSSparams iop("cat", EXODUS_DATABASE_TYPE, iossProp);

  Ioss::DatabaseIO *cat_d = bmSet.getCatalystDatabase(iop);

  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize on read"; }
  Ioss::Region cat_reg(cat_d);
  
  Ioss::Region exo_reg(exo_d);
  
  auto cat_elemBlock = cat_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));
  auto exo_elemBlock = exo_reg.get_element_block(bmSet.getUnstructuredBlockName(bm.getID()));

  checkEntityContainerZeroCopyFields(cat_reg.get_element_blocks());

  bool exo_foo_exists = exo_elemBlock->field_exists("foo");
  bool cat_foo_exists = cat_elemBlock->field_exists("foo");
  EXPECT_TRUE(exo_foo_exists);
  EXPECT_TRUE(cat_foo_exists);
  if(exo_foo_exists && cat_foo_exists) 
    EXPECT_TRUE(exo_elemBlock->get_field("foo") == cat_elemBlock->get_field("foo"));
  
  //Check field data for equality
  auto cat_field = cat_elemBlock->get_fieldref("foo");
  std::vector<std::byte> dcBuffer(cat_field.get_size());
  cat_elemBlock->get_field_data("foo", Data(dcBuffer), dcBuffer.size());

  exo_reg.begin_state(1);
  auto exo_field = exo_elemBlock->get_fieldref("foo");
  std::vector<std::byte> deBuffer(exo_field.get_size());
  exo_elemBlock->get_field_data("foo", Data(deBuffer), deBuffer.size());
  EXPECT_EQ(dcBuffer, deBuffer);

}

//Read from file. Can. Or available exodus file in test suite.
TEST_F(Iocatalyst_DatabaseIOTest, Exodus_Prop_SURFACE_SPLIT_TYPE)
{
  Iocatalyst::BlockMesh bm;
  setBlockMeshSize(2, 2, 2);

  addBlockMesh(bm);

  //If write and read have same split type, we can handle. Else no.
  Ioss::PropertyManager iossProp;
  iossProp.add(Ioss::Property("SURFACE_SPLIT_TYPE", "BLOCK"));

  std::string exoFile = "can.ex2";
  
  conduit_cpp::Node c_node = getConduitFromExodusFile(exoFile, iossProp);

  Ioss::DatabaseIO *cat_d = getCatalystDatabaseFromConduit(c_node, iossProp);
  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize"; }
 
  Ioss::Region cat_reg(cat_d);
  
  Ioss::SideSetContainer cat_sideSets = cat_reg.get_sidesets();
  checkEntityContainerZeroCopyFields(cat_sideSets);
  
  EXPECT_TRUE(cat_sideSets.empty())<<"Cat sidesets not empty when different SURFACE_SPLIT_TYPE";

  Ioss::PropertyManager iossProp_s;
  iossProp_s.add(Ioss::Property("SURFACE_SPLIT_TYPE", "TOPOLOGY"));
  cat_d = getCatalystDatabaseFromConduit(c_node, iossProp_s);

  if(cat_d == nullptr){ EXPECT_TRUE(false) << "Catalyst db unable to initialize"; }

  Ioss::Region cat_reg_same(cat_d);
  
  cat_sideSets = cat_reg_same.get_sidesets();

  EXPECT_TRUE(!cat_sideSets.empty())<<"Cat sidesets empty when identical SURFACE_SPLIT_TYPE";

}
