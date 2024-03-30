// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "Ioss_PropertyManager.h"
#include "catch.hpp"
#include <Iovs_CatalystLogging.h>

TEST_CASE_METHOD(CatalystTestFixture, "CatalystLoggingDefault", "[catalyst logging]")
{
  std::string           logFileName = "log.csv";
  Ioss::PropertyManager props;
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_FILE_NAME", logFileName));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "Log Info"));
  props.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", -33));
  runCatalystLoggingTest(&props, "test1.json", "block_crush_1.ex2");
  checkTestOutputFileDoesNotExist(Iovs::CatalystLogging::getDefaultLogFileName().c_str());
  checkTestOutputFileExists("CatalystOutput/test1.0010.png");
  checkTestOutputFileExists(logFileName.c_str());
}
