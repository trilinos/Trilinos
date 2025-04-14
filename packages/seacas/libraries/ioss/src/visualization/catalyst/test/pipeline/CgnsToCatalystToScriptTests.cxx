// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture,
    "SimpleAeroBluntWedgeCgnsTest1", "[cgns to catalyst script]") {

    runPhactoriJSONTest("test2.json", "aero_blunt_wedge_test3.cgns");
    checkTestOutputFileExists("CatalystOutput/test2.0000.png");
}
