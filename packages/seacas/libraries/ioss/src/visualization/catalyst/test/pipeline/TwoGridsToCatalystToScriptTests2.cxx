// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest2_ex2_cgns", "[cgns to catalyst script]") {

    runPhactoriJSONTestTwoGrid("test4.json", "block_crush_1.ex2", "aero_blunt_wedge_test3.cgns");
    checkTestOutputFileExists("CatalystOutput_test4/test4_inputA.0000.png");
    checkTestOutputFileExists("CatalystOutput_test4/test4_inputB.0000.png");
}
