// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture,
    "SimpleBlockCrushTest1", "[exodus to catalyst script]") {

    runPhactoriJSONTest("test1.json", "block_crush_1.ex2");
    checkTestOutputFileExists("CatalystOutput/test1.0010.png");
}
