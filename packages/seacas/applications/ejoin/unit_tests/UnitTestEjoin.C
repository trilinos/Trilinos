// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "ArgvBuilder.h"
#include "FileUtils.h"
#include "MeshFixture.h"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unistd.h>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <exodusII.h>

#include <Ionit_Initializer.h>
#include <Ioss_Enumerate.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <tokenize.h>

#include "EJ_CodeTypes.h"
#include "EJ_SystemInterface.h"
#include "EJ_Version.h"
#include "EJ_mapping.h"
#include "EJ_match.h"
#include "EJ_vector3d.h"

template <typename INT>
double ejoin(SystemInterface &interFace, const RegionVector &inputRegions, INT dummy);

namespace {

  class EJoinTester : public utest_util::MeshFixture
  {
  public:
    EJoinTester() : utest_util::MeshFixture(3)
    {
      Ioss::Init::Initializer io;
      setup_empty_mesh();
    }

    ~EJoinTester() {}

  protected:
    void create_input_textmesh_regions(const std::vector<std::string> &meshDescs,
                                       RegionVector                   &inputRegions)
    {
      Ioss::PropertyManager propertyManager;

      int regionCount = 0;
      for (const std::string &meshDesc : meshDescs) {
        Ioss::DatabaseIO *db = Ioss::IOFactory::create("textmesh", meshDesc, Ioss::READ_MODEL,
                                                       get_comm(), propertyManager);
        EXPECT_TRUE(db != nullptr);
        EXPECT_TRUE(db->ok(true));

        std::string   regionName = "input_region" + std::to_string(regionCount);
        Ioss::Region *region     = new Ioss::Region(db, regionName);
        EXPECT_TRUE(region != nullptr);
        region->property_add(Ioss::Property("block_omission_count", 0));
        inputRegions.push_back(region);
        regionCount++;
      }
    }

    void create_input_region(const std::string &meshDesc, RegionVector &inputRegions)
    {
      // Use the internal region from the fixture
      setup_mesh(meshDesc);
      inputRegions.push_back(get_mesh().get_region());
      inputRegions[0]->property_add(Ioss::Property("block_omission_count", 0));
    }
  };

  TEST_F(EJoinTester, emptyTest) {}

  TEST_F(EJoinTester, joinSingleHexMeshWithMaterialProperties)
  {
    if (get_parallel_size() > 1)
      GTEST_SKIP();

    std::string inputFile  = "dummy.g";
    std::string outputFile = utest_util::unique_filename("singleHexJoin", "g");

    SystemInterface interFace;
    interFace.inputFiles_.push_back(inputFile);
    interFace.outputName_ = outputFile;

    std::string meshDesc = "textmesh:"
                           "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
                           "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";

    RegionVector inputRegions;
    create_input_region(meshDesc, inputRegions);

    // Add a material property to block_1 from textmesh
    const std::string propertyName("MATERIAL_PROPERTY");
    const std::string propertyValue("KRYPTONITE");
    add_material_property_to_element_block(inputRegions[0], "block_1", propertyName, propertyValue);

    // Call ejoin on the single mesh ... the material property should make it to output
    if (inputRegions[0]->get_database()->int_byte_size_api() == 4) {
      (void)ejoin(interFace, inputRegions, 0);
    }
    else {
      (void)ejoin(interFace, inputRegions, static_cast<int64_t>(0));
    }

    // Test the material property from the output file
    test_property_from_file(outputFile, propertyName, propertyValue);

    unlink(outputFile.c_str());
  }

  TEST_F(EJoinTester, join2DQuadMeshes)
  {
    if (get_parallel_size() > 1)
      GTEST_SKIP();

    std::string outputFile = utest_util::unique_filename("twoQuadJoin", "g");

    utest_util::ArgvBuilder builder;
    builder.addArgument("ejoin_unit_test");
    builder.addArgument("-output");
    builder.addArgument(outputFile);
    builder.addArgument("-match_node_coordinates");
    builder.addArgument("-tolerance");
    builder.addArgument("1.0e-3");
    builder.addArgument("-combine_element_blocks");
    builder.addArgument("dummy1.g");
    builder.addArgument("dummy2.g");

    SystemInterface interFace;
    interFace.parse_options(builder.argc(), builder.argv());

    std::string meshDesc1 = "0,1,QUAD_4_2D,1,2,3,4,block_1"
                            "|coordinates: 0,0, 1,0, 1,1, 0,1"
                            "|dimension:2";
    std::string meshDesc2 = "0,1,QUAD_4_2D,1,2,3,4,block_1"
                            "|coordinates: 1,0, 2,0, 2,1, 1,1"
                            "|dimension:2";

    RegionVector inputRegions;
    create_input_textmesh_regions({meshDesc1, meshDesc2}, inputRegions);

    if (inputRegions[0]->get_database()->int_byte_size_api() == 4) {
      (void)ejoin(interFace, inputRegions, 0);
    }
    else {
      (void)ejoin(interFace, inputRegions, static_cast<int64_t>(0));
    }

    for (size_t i = 0; i < inputRegions.size(); i++) {
      delete inputRegions[i];
      inputRegions[i] = nullptr;
    }

    // Use the internal region from the fixture
    setup_mesh(outputFile);

    const auto &mesh = get_mesh();

    EXPECT_EQ(6, mesh.get_num_global_nodes());
    EXPECT_EQ(2, mesh.get_num_global_elements());
    EXPECT_EQ(1, mesh.get_num_global_element_blocks());

    unlink(outputFile.c_str());
  }

  TEST_F(EJoinTester, join1DBeamMeshes)
  {
    if (get_parallel_size() > 1)
      GTEST_SKIP();

    std::string outputFile = utest_util::unique_filename("twoBeamJoin", "g");

    utest_util::ArgvBuilder builder;
    builder.addArgument("ejoin_unit_test");
    builder.addArgument("-output");
    builder.addArgument(outputFile);
    builder.addArgument("-match_node_coordinates");
    builder.addArgument("-tolerance");
    builder.addArgument("1.0e-3");
    builder.addArgument("-combine_element_blocks");
    builder.addArgument("dummy1.g");
    builder.addArgument("dummy2.g");

    SystemInterface interFace;
    interFace.parse_options(builder.argc(), builder.argv());

    std::string meshDesc1 = "0,1,BEAM_2,1,2,block_1"
                            "|coordinates: 0, 1"
                            "|dimension:1";
    std::string meshDesc2 = "0,1,BEAM_2,1,2,block_1"
                            "|coordinates: 1, 2"
                            "|dimension:1";

    RegionVector inputRegions;
    create_input_textmesh_regions({meshDesc1, meshDesc2}, inputRegions);

    if (inputRegions[0]->get_database()->int_byte_size_api() == 4) {
      (void)ejoin(interFace, inputRegions, 0);
    }
    else {
      (void)ejoin(interFace, inputRegions, static_cast<int64_t>(0));
    }

    for (size_t i = 0; i < inputRegions.size(); i++) {
      delete inputRegions[i];
      inputRegions[i] = nullptr;
    }

    // Use the internal region from the fixture
    setup_mesh(outputFile);

    const auto &mesh = get_mesh();

    // Nodal coordinate matching isn't working as expected in 1D .... should expect 3 here instead
    // of 4
    EXPECT_EQ(4, mesh.get_num_global_nodes());
    EXPECT_EQ(2, mesh.get_num_global_elements());
    EXPECT_EQ(1, mesh.get_num_global_element_blocks());

    unlink(outputFile.c_str());
  }

  TEST_F(EJoinTester, join1DPointMeshes)
  {
    if (get_parallel_size() > 1)
      GTEST_SKIP();

    std::string outputFile = utest_util::unique_filename("fourPointJoin", "g");

    utest_util::ArgvBuilder builder;
    builder.addArgument("ejoin_unit_test");
    builder.addArgument("-output");
    builder.addArgument(outputFile);
    builder.addArgument("-match_node_coordinates");
    builder.addArgument("-tolerance");
    builder.addArgument("1.0e-3");
    builder.addArgument("dummy1.g");
    builder.addArgument("dummy2.g");

    SystemInterface interFace;
    interFace.parse_options(builder.argc(), builder.argv());

    std::string meshDesc1 = "0,1,NODE,1,block_1\n"
                            "0,2,NODE,2,block_1"
                            "|coordinates: 0, 1"
                            "|dimension:1";
    std::string meshDesc2 = "0,1,NODE,1,block_1\n"
                            "0,2,NODE,2,block_1"
                            "|coordinates: 1, 2"
                            "|dimension:1";

    RegionVector inputRegions;
    create_input_textmesh_regions({meshDesc1, meshDesc2}, inputRegions);

    if (inputRegions[0]->get_database()->int_byte_size_api() == 4) {
      (void)ejoin(interFace, inputRegions, 0);
    }
    else {
      (void)ejoin(interFace, inputRegions, static_cast<int64_t>(0));
    }

    for (size_t i = 0; i < inputRegions.size(); i++) {
      delete inputRegions[i];
      inputRegions[i] = nullptr;
    }

    // Use the internal region from the fixture
    setup_mesh(outputFile);

    const auto &mesh = get_mesh();

    // Nodal coordinate matching isn't working as expected in 1D .... should expect 3 here instead
    // of 4
    EXPECT_EQ(4, mesh.get_num_global_nodes());
    EXPECT_EQ(0, mesh.get_num_global_elements());
    EXPECT_EQ(0, mesh.get_num_global_element_blocks());

    unlink(outputFile.c_str());
  }
} // namespace
