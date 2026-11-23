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
#include "elb_nem_slice.h"
#include "unistd.h"

namespace {

  class NemSliceTester : public utest_util::MeshFixture
  {
  public:
    NemSliceTester() : utest_util::MeshFixture(3)
    {
      Ioss::Init::Initializer io;
      setup_empty_mesh();
    }

    ~NemSliceTester() {}

    void SetUp() override
    {
      optind = 1;
#ifdef APPLE
      optreset = 1;
#endif
    }

    void TearDown() override
    {
      optind = 1;
#ifdef APPLE
      optreset = 1;
#endif
    }
  };

  TEST_F(NemSliceTester, emptyTest) {}

  TEST_F(NemSliceTester, beamOneDimSideIds)
  {
    std::string inputFile  = "dummy.g";
    std::string outputFile = utest_util::unique_filename("beamOneDimSideIds", "g");

    std::string meshDesc =
        "textmesh:"
        "0,7,BEAM_2,4,1,block_2\n"
        "0,5,BEAM_2,1,9,block_2\n"
        "0,6,BEAM_2,9,8,block_2\n"
        "0,8,BEAM_2,8,7,block_2\n"
        "0,3,BEAM_2,7,6,block_1\n"
        "0,1,BEAM_2,6,3,block_1\n"
        "0,2,BEAM_2,3,5,block_1\n"
        "0,4,BEAM_2,5,2,block_1\n"
        "|coordinates:   -0.874923313523287, 1.5, 0.80376037271592,\n"
        "                -1.2, 1.1935484824723, 0.52582860264049,\n"
        "                0.111515502329254, -0.153258442252325, -0.500565877886665\n"
        "|dimension:   1";

    setup_mesh(meshDesc);

    Ioss::PropertyManager props;
    write_region_to_file(get_mesh().get_region(), props, outputFile);

    utest_util::ArgvBuilder builder;
    builder.addArgument("nem_slice_unit_test");
    builder.addArgument("-e");
    builder.addArgument("-m");
    builder.addArgument("mesh=2");
    builder.addArgument("-l");
    builder.addArgument("random");
    builder.addArgument("-o");
    builder.addArgument("edge2.Mesh.8.nem");
    builder.addArgument(outputFile.c_str());

    NemSlice<int> ns(builder.argc(), builder.argv());

    ns.setup();
    ns.find_graph();
    ns.find_load_balance();

    std::vector<std::vector<int>> expected = {{2, 1}, {1, 2}};

    EXPECT_EQ(ns.lb.e_cmap_sides, expected);
  }

} // namespace
