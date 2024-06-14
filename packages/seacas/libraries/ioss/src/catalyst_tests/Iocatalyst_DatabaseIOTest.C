// Copyright(C) 1999-2021 National Technology & Engineering Solutions
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
#include <catalyst/Iocatalyst_CatalystManager.h>
#include <catalyst/Iocatalyst_DatabaseIO.h>
#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>
#include <stdlib.h>

Iocatalyst_DatabaseIOTest::Iocatalyst_DatabaseIOTest()
{
  part.id         = putils.parallel_rank();
  part.size       = putils.parallel_size();
  blockMeshSize.i = 2;
  blockMeshSize.j = 2;
  blockMeshSize.k = 2;
  origin.i        = 0;
  origin.j        = 0;
  origin.k        = 0;
}

bool Iocatalyst_DatabaseIOTest::regionsAreEqual(const std::string &fileName,
                                                const std::string &catFileName,
                                                const std::string &iossDatabaseType)
{
  Ioss::PropertyManager dbProps;
  auto                  inputFileName         = fileName;
  auto                  inputCatalystFileName = catFileName;
  Ioss::ParallelUtils   pu;
  int                   numRanks = pu.parallel_size();
  int                   rank     = pu.parallel_rank();
  if (iossDatabaseType == EXODUS_DATABASE_TYPE && numRanks > 1) {
    inputFileName += "." + std::to_string(numRanks) + "." + std::to_string(rank);
    inputCatalystFileName += "." + std::to_string(numRanks) + "." + std::to_string(rank);
  }
  Ioss::DatabaseIO *dbi =
      Ioss::IOFactory::create(iossDatabaseType, inputFileName, Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_self(), dbProps);
  if (dbi == nullptr || !dbi->ok(true)) {
    return false;
  }

  Ioss::PropertyManager dbCatProps;
  Ioss::DatabaseIO     *dbiCat =
      Ioss::IOFactory::create(iossDatabaseType, inputCatalystFileName, Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_self(), dbCatProps);
  if (dbiCat == nullptr || !dbiCat->ok(true)) {
    return false;
  }

  Ioss::Region          ir(dbi);
  Ioss::Region          rCat(dbiCat);
  Ioss::MeshCopyOptions options;
  options.data_storage_type = 1;
  return Ioss::Compare::compare_database(ir, rCat, options);
}

void Iocatalyst_DatabaseIOTest::runStructuredTest(const std::string &testName)
{
  std::string cgnsFileName =
      testName + CATALYST_TEST_FILE_NP + std::to_string(part.size) + CGNS_FILE_EXTENSION;
  std::string catalystFileName = CATALYST_TEST_FILE_PREFIX + testName + CATALYST_TEST_FILE_NP +
                                 std::to_string(part.size) + CGNS_FILE_EXTENSION;
  Iocatalyst::BlockMeshSet::IOSSparams iop(cgnsFileName, CGNS_DATABASE_TYPE);
  bmSet.writeIOSSFile(iop);
  iop.fileName = catalystFileName;
  bmSet.writeCatalystIOSSFile(iop);
  checkZeroCopyFields(iop);
  EXPECT_TRUE(regionsAreEqual(cgnsFileName, catalystFileName, CGNS_DATABASE_TYPE));
}

void Iocatalyst_DatabaseIOTest::runUnstructuredTest(const std::string &testName)
{
  std::string exodusFileName =
      testName + CATALYST_TEST_FILE_NP + std::to_string(part.size) + EXODUS_FILE_EXTENSION;
  std::string catalystFileName = CATALYST_TEST_FILE_PREFIX + testName + CATALYST_TEST_FILE_NP +
                                 std::to_string(part.size) + EXODUS_FILE_EXTENSION;
  Iocatalyst::BlockMeshSet::IOSSparams iop(exodusFileName, EXODUS_DATABASE_TYPE);
  bmSet.writeIOSSFile(iop);
  iop.fileName = catalystFileName;
  bmSet.writeCatalystIOSSFile(iop);
  checkZeroCopyFields(iop);
  EXPECT_TRUE(regionsAreEqual(exodusFileName, catalystFileName, EXODUS_DATABASE_TYPE));
}

Ioss::DatabaseIO *
Iocatalyst_DatabaseIOTest::writeAndGetExodusDatabaseOnRead(const std::string    &testName,
                                                           Ioss::PropertyManager dbProps)
{
  std::string exodusFileName =
      testName + CATALYST_TEST_FILE_NP + std::to_string(part.size) + EXODUS_FILE_EXTENSION;
  Iocatalyst::BlockMeshSet::IOSSparams iop(exodusFileName, EXODUS_DATABASE_TYPE);
  bmSet.writeIOSSFile(iop);
  Ioss::DatabaseIO *exo_db =
      getDatabaseOnReadFromFileName(exodusFileName, EXODUS_DATABASE_TYPE, dbProps);
  if (exo_db == nullptr) {
    EXPECT_TRUE(false) << "Exodus db unable to initialize on read";
  }
  return exo_db;
}

Ioss::DatabaseIO *
Iocatalyst_DatabaseIOTest::getExodusDatabaseFromFile(std::string          &filename,
                                                     Ioss::PropertyManager dbProps)
{
  Ioss::PropertyManager edbProps(dbProps);

  std::string inputFileName = filename;

  Ioss::DatabaseIO *edbi =
      Ioss::IOFactory::create(EXODUS_DATABASE_TYPE, inputFileName, Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_self(), edbProps);
  if (edbi == nullptr || !edbi->ok(true)) {
    return nullptr;
  }

  return edbi;
}

Ioss::DatabaseIO *
Iocatalyst_DatabaseIOTest::getCatalystDatabaseFromConduitFiles(const std::string    &dirName,
                                                               Ioss::PropertyManager dbProps)
{
  setenv("CATALYST_DATA_DUMP_DIRECTORY", dirName.c_str(), 1);
  Ioss::DatabaseIO *cdbi =
      Ioss::IOFactory::create(CATALYST_DATABASE_TYPE, "catalyst.bin", Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_self(), dbProps);
  if (cdbi == nullptr || !cdbi->ok(true)) {
    return nullptr;
  }

  return cdbi;
}

conduit_cpp::Node Iocatalyst_DatabaseIOTest::getConduitFromExodusFile(std::string &filename,
                                                                      Ioss::PropertyManager dbProps)
{
  Iocatalyst::CatalystManager::getInstance().reset();

  Ioss::PropertyManager edbProps;
  edbProps.add(Ioss::Property("SURFACE_SPLIT_TYPE", "TOPOLOGY"));
  Ioss::DatabaseIO *edbi = getExodusDatabaseFromFile(filename, edbProps);

  // Create Cat Db on write
  Ioss::PropertyManager cdbwProps(edbi->get_property_manager());
  Ioss::DatabaseIO     *cdb_on_write =
      Ioss::IOFactory::create(CATALYST_DATABASE_TYPE, CATALYST_DUMMY_DATABASE, Ioss::WRITE_RESULTS,
                              Ioss::ParallelUtils::comm_world(), cdbwProps);
  if (cdb_on_write == nullptr || !cdb_on_write->ok(true)) {
    return conduit_cpp::Node();
  }

  Ioss::Region          cor(edbi);
  Ioss::Region          cir(cdb_on_write);
  Ioss::MeshCopyOptions options;
  options.data_storage_type = 1;
  Ioss::copy_database(cor, cir, options);

  auto c_node = reinterpret_cast<conduit_node *>(
      ((Iocatalyst::DatabaseIO *)cdb_on_write)->get_catalyst_conduit_node());
  conduit_cpp::Node conduitNode;
  auto              cpp_node = conduit_cpp::cpp_node(c_node);
  conduitNode.set(cpp_node);
  return conduitNode;
}

Ioss::DatabaseIO *
Iocatalyst_DatabaseIOTest::getCatalystDatabaseFromConduit(conduit_cpp::Node    &conduitNode,
                                                          Ioss::PropertyManager dbProps)
{

  Ioss::PropertyManager cdbrProps = Ioss::PropertyManager(dbProps);
  cdbrProps.add(Ioss::Property("CATALYST_CONDUIT_NODE", conduit_cpp::c_node(&conduitNode)));

  // Give to Cat Db on read
  Ioss::DatabaseIO *cdb_on_read =
      Ioss::IOFactory::create(CATALYST_DATABASE_TYPE, CATALYST_DUMMY_DATABASE, Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_world(), cdbrProps);
  if (cdb_on_read == nullptr || !cdb_on_read->ok(true)) {
    return nullptr;
  }

  return cdb_on_read;
}

Ioss::DatabaseIO *Iocatalyst_DatabaseIOTest::getDatabaseOnReadFromFileName(
    const std::string &fileName, const std::string &iossDatabaseType, Ioss::PropertyManager dbProps)
{
  Ioss::PropertyManager dbaseProps = Ioss::PropertyManager(dbProps);
  auto                inputFileName = fileName;
  Ioss::ParallelUtils pu;
  int                 numRanks = pu.parallel_size();
  int                 rank     = pu.parallel_rank();
  if (iossDatabaseType == EXODUS_DATABASE_TYPE && numRanks > 1) {
    inputFileName += "." + std::to_string(numRanks) + "." + std::to_string(rank);
  }
  Ioss::DatabaseIO *dbi =
      Ioss::IOFactory::create(iossDatabaseType, inputFileName, Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_self(), dbaseProps);
  if (dbi == nullptr || !dbi->ok(true)) {
    return nullptr;
  }
  return dbi;
}

void Iocatalyst_DatabaseIOTest::checkZeroCopyFields(Iocatalyst::BlockMeshSet::IOSSparams &iop)
{
  Ioss::PropertyManager cdbProps;
  cdbProps.add(Ioss::Property("CATALYST_CONDUIT_NODE", iop.getCatalystConduitNode()));

  Ioss::DatabaseIO *cdbi =
      Ioss::IOFactory::create(Iocatalyst::BlockMeshSet::CATALYST_DATABASE_TYPE,
                              Iocatalyst::BlockMeshSet::CATALYST_DUMMY_DATABASE, Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_world(), cdbProps);
  if (cdbi == nullptr || !cdbi->ok(true)) {
    return;
  }

  Ioss::Region                            cir(cdbi);
  Iocatalyst::DatabaseIO::RegionContainer rc;
  rc.push_back(&cir);
  checkEntityContainerZeroCopyFields(rc);
  checkEntityContainerZeroCopyFields(cir.get_node_blocks());
  checkEntityContainerZeroCopyFields(cir.get_element_blocks());
  checkEntityContainerZeroCopyFields(cir.get_structured_blocks());
}

void Iocatalyst_DatabaseIOTest::setBlockMeshSize(unsigned int i, unsigned int j, unsigned int k)
{
  blockMeshSize.i = i;
  blockMeshSize.j = j;
  blockMeshSize.k = k;
}

void Iocatalyst_DatabaseIOTest::setOrigin(unsigned int i, unsigned int j, unsigned int k)
{
  origin.i = i;
  origin.j = j;
  origin.k = k;
}

void Iocatalyst_DatabaseIOTest::addBlockMesh(Iocatalyst::BlockMesh &blockMesh)
{
  blockMesh.init(part, blockMeshSize, origin);
  bmSet.addBlockMesh(blockMesh);
}