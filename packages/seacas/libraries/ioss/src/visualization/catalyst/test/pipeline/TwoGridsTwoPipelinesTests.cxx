// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest1_cgns_ex2", "[cgns to catalyst script]") {
    //runPhactoriJSONTestTwoGridTwoPipe("test9a.json", "aero_blunt_wedge_test3.cgns",
    //        "test9b.json", "block_crush_1.ex2");
    runPhactoriJSONTestTwoGridTwoPipe(
            "test9a.json", "aero_blunt_wedge_test3.cgns",
            "test9b.json", "block_crush_1.ex2"
            );
    checkTestOutputFileExists("CatalystOutput_test9/test9_inputA.0000.png");
    checkTestOutputFileExists("CatalystOutput_test9/test9_inputB.0000.png");
}

