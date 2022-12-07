// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"

TEST_CASE_METHOD(CatalystTestFixture,
    "CatalystV2ScriptFromParaViewGuiTest1", "[paraview gui generated catalyst script]") {

    //runPhactoriJSONTestTwoGrid("test7.json", "aero_blunt_wedge_test3.cgns", "aero_blunt_wedge_wall_test3.ex2");
    runParaViewGuiScriptTest("aero_blunt_wedge_pv590_script_1.py", "aero_blunt_wedge_test3.cgns");
    checkTestOutputFileExists("test8_000004.png");
}
