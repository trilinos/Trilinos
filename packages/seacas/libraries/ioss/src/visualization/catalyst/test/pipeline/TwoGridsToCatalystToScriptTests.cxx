// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest1_cgns_ex2", "[cgns to catalyst script]") {

    runPhactoriJSONTestTwoGrid("test3.json", "aero_blunt_wedge_test3.cgns", "block_crush_1.ex2");
    checkTestOutputFileExists("CatalystOutput_test3/test3_inputA.0000.png");
    checkTestOutputFileExists("CatalystOutput_test3/test3_inputB.0000.png");
}

/*
TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest2_ex2_cgns", "[cgns to catalyst script]") {

    runPhactoriJSONTestTwoGrid("test4.json", "block_crush_1.ex2", "aero_blunt_wedge_test3.cgns");
    checkTestOutputFileExists("CatalystOutput_test4/test4.0000.png");
}

TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest3_cgns_cgns", "[cgns to catalyst script]") {

    runPhactoriJSONTestTwoGrid("test5.json", "aero_blunt_wedge_test3.cgns", "aero_blunt_wedge_test3.cgns");
    checkTestOutputFileExists("CatalystOutput_test5/test5.0000.png");
}

TEST_CASE_METHOD(CatalystTestFixture,
    "TwoGridInputTest4_ex2_ex2", "[cgns to catalyst script]") {

    runPhactoriJSONTestTwoGrid("test6.json", "block_crush_1.ex2", "block_crush_1.ex2");
    checkTestOutputFileExists("CatalystOutput/test6.0000.png");
}
*/
