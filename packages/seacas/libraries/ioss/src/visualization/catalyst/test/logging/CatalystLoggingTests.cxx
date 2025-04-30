// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <Iovs_CatalystLogging.h>
#include <iostream>
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(Iovs::CatalystLogging, "CatalystLoggingDefault", "[catalystLogging]")
{

  REQUIRE(isCatalystLoggingON() == false);
  REQUIRE(getLogFileName() == getDefaultLogFileName());
  REQUIRE(getLogOutputDirectoryPath() == getDefaultLogOutputDirectoryPath());
}

TEST_CASE_METHOD(Iovs::CatalystLogging, "CatalystLoggingEnabled", "[catalystLogging]")
{
  Ioss::PropertyManager props;
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_FILE_NAME", "foo.csv"));
  props.add(Ioss::Property("CATALYST_LOGGING_OUTPUT_DIRECTORY_PATH", "/projects/bar/"));
  setProperties(&props);

  REQUIRE(isCatalystLoggingON() == true);
  REQUIRE(getLogFileName() == "foo.csv");
  REQUIRE(getLogOutputDirectoryPath() == "/projects/bar/");
}

TEST_CASE_METHOD(Iovs::CatalystLogging, "CatalystLoggingFileInvalid", "[catalystLogging]")
{
  Ioss::PropertyManager props;
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_OUTPUT_DIRECTORY_PATH", "/dev/null/invalid/"));
  setProperties(&props);
  REQUIRE_NOTHROW(writeToLogFile());
}

TEST_CASE_METHOD(Iovs::CatalystLogging, "CatalystLoggingWrite", "[catalystLogging]")
{
  Ioss::PropertyManager props;
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  remove(getLogFileName().c_str());

  SECTION("Write Log Default")
  {
    setProperties(&props);

    if (isCatalystLoggingON()) {
      writeToLogFile();
    }

    CatalystTestFixture::checkTestOutputFileExists(getLogFileName().c_str());
    std::vector<std::vector<std::string>> lfc = readLogFile();

    REQUIRE(lfc.size() == 0);
  }

  SECTION("Write Log Once With All Supported Property Types")
  {
    props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
    props.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 6));
    props.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 3.7556));
    setProperties(&props);

    std::vector<std::string> logLine;
    if (isCatalystLoggingON()) {
      logLine = writeToLogFile();
    }

    CatalystTestFixture::checkTestOutputFileExists(getLogFileName().c_str());
    std::vector<std::vector<std::string>> lfc = readLogFile();

    REQUIRE(lfc.size() == 2);
    REQUIRE(lfc[0] == getLogFileHeaders());
    REQUIRE(lfc[1] == logLine);
  }

  SECTION("Write Log Twice With All Supported Property Types")
  {
    props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "bar"));
    props.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 12));
    props.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 7.45454));
    setProperties(&props);

    std::vector<std::string> logLineOne;
    if (isCatalystLoggingON()) {
      logLineOne = writeToLogFile();
    }

    Ioss::PropertyManager propsTwo;
    propsTwo.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
    propsTwo.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 90.3));
    propsTwo.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
    propsTwo.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 107));
    setProperties(&propsTwo);

    std::vector<std::string> logLineTwo;
    if (isCatalystLoggingON()) {
      logLineTwo = writeToLogFile();
    }

    CatalystTestFixture::checkTestOutputFileExists(getLogFileName().c_str());
    std::vector<std::vector<std::string>> lfc = readLogFile();

    REQUIRE(lfc.size() == 3);
    REQUIRE(lfc[0] == getLogFileHeaders());
    REQUIRE(lfc[1] == logLineOne);
    REQUIRE(lfc[2] == logLineTwo);
  }

  SECTION("Write Log With Unsupported Property Types")
  {
    props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "bar"));
    int ten = 10;
    props.add(Ioss::Property("CATALYST_LOGGING_POINTER_PROP", (void *)&ten));
    std::vector<int> iv = {1, 2, 3};
    props.add(Ioss::Property("CATALYST_LOGGING_VEC_INTEGER_PROP", iv));
    std::vector<double> dv = {1.4, 2.2, -56.3};
    props.add(Ioss::Property("CATALYST_LOGGING_VEC_DOUBLE_PROP", dv));
    setProperties(&props);

    std::vector<std::string> logLine;
    if (isCatalystLoggingON()) {
      logLine = writeToLogFile();
    }

    CatalystTestFixture::checkTestOutputFileExists(getLogFileName().c_str());
    std::vector<std::vector<std::string>> lfc = readLogFile();

    std::vector<std::string> h = {"STRING_PROP"};
    REQUIRE(h == getLogFileHeaders());
    REQUIRE(lfc.size() == 2);
    REQUIRE(lfc[0] == getLogFileHeaders());
    REQUIRE(lfc[1] == logLine);
  }

  SECTION("Write Log with quoted commas")
  {
    props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
    props.add(Ioss::Property("CATALYST_LOGGING_QUOTED_COMMA_PROP", "\"8,7,go,3.2\""));
    props.add(Ioss::Property("CATALYST_LOGGING_QUOTED_ANOTHER_COMMA_PROP",
                             "This address \", 343 Far Lane, ND 8732, RT 3\""));
    setProperties(&props);

    std::vector<std::string> logLine;
    if (isCatalystLoggingON()) {
      logLine = writeToLogFile();
    }

    CatalystTestFixture::checkTestOutputFileExists(getLogFileName().c_str());
    std::vector<std::vector<std::string>> lfc = readLogFile();

    REQUIRE(lfc.size() == 2);
    REQUIRE(lfc[0] == getLogFileHeaders());
    REQUIRE(lfc[1] == logLine);
  }

  SECTION("Write Log with quoted quotes")
  {
    props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "bar"));
    props.add(Ioss::Property("CATALYST_LOGGING_QUOTED_QUOTE_PROP", "\"\"foo\"\" bar \"1,3,2\""));
    props.add(
        Ioss::Property("CATALYST_LOGGING_ANOTHER_QUOTED_QUOTE_PROP", "This has \"\"quote\"\""));

    setProperties(&props);

    std::vector<std::string> logLine;
    if (isCatalystLoggingON()) {
      logLine = writeToLogFile();
    }

    CatalystTestFixture::checkTestOutputFileExists(getLogFileName().c_str());
    std::vector<std::vector<std::string>> lfc = readLogFile();

    REQUIRE(lfc.size() == 2);
    REQUIRE(lfc[0] == getLogFileHeaders());
    REQUIRE(lfc[1] == logLine);
  }

  remove(getLogFileName().c_str());
}
