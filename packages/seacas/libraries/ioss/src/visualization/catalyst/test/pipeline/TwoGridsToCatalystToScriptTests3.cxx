// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest3_wedge_volume_and_wall", "[volume and wall to catalyst script]") {

    //runPhactoriJSONTestTwoGrid("test7.json", "aero_blunt_wedge_test3.cgns", "aero_blunt_wedge_wall_test3.ex2");
    runPhactoriJSONTestTwoGrid("test7.json", "aero_blunt_wedge_test3.cgns", "aero_blunt_wedge_test3.cgns");
    checkTestOutputFileExists("CatalystOutput_test7/test7_inputA.0000.png");
    checkTestOutputFileExists("CatalystOutput_test7/test7_inputB.0000.png");
}
