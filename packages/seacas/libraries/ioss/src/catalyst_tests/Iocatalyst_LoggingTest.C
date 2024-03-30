// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <catalyst_tests/Iocatalyst_LoggingTest.h>

TEST_F(LoggingTest, Defaults)
{
  EXPECT_FALSE(log.isCatalystLoggingON());
  EXPECT_EQ(log.getLogFileName(), log.getDefaultLogFileName());
  EXPECT_EQ(log.getLogOutputDirectoryPath(), log.getDefaultLogOutputDirectoryPath());
}

TEST_F(LoggingTest, CatalystLoggingEnabled)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_FILE_NAME", "foo.csv"));
  props.add(Ioss::Property("CATALYST_LOGGING_OUTPUT_DIRECTORY_PATH", "/projects/bar/"));
  log.setProperties(&props);

  EXPECT_TRUE(log.isCatalystLoggingON());
  EXPECT_EQ(log.getLogFileName(), "foo.csv");
  EXPECT_EQ(log.getLogOutputDirectoryPath(), "/projects/bar/");
}

TEST_F(LoggingTest, WriteLogDefault)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  log.setProperties(&props);

  if (log.isCatalystLoggingON()) {
    log.writeToLogFile();
  }

  checkTestOutputFileExists(log.getLogFileName().c_str());
  std::vector<std::vector<std::string>> lfc = log.readLogFile();

  EXPECT_EQ(lfc.size(), 0);
}

TEST_F(LoggingTest, WriteLogOnceWithAllSupportedPropertyTypes)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
  props.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 6));
  props.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 3.7556));
  log.setProperties(&props);

  std::vector<std::string> logLine;
  if (log.isCatalystLoggingON()) {
    logLine = log.writeToLogFile();
  }

  checkTestOutputFileExists(log.getLogFileName().c_str());
  std::vector<std::vector<std::string>> lfc = log.readLogFile();

  ASSERT_EQ(lfc.size(), 2);
  EXPECT_EQ(lfc[0], log.getLogFileHeaders());
  EXPECT_EQ(lfc[1], logLine);
}

TEST_F(LoggingTest, WriteLogTwiceWithAllSupportedPropertyTypes)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "bar"));
  props.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 12));
  props.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 7.45454));
  log.setProperties(&props);

  std::vector<std::string> logLineOne;
  if (log.isCatalystLoggingON()) {
    logLineOne = log.writeToLogFile();
  }

  Ioss::PropertyManager propsTwo;
  propsTwo.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  propsTwo.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 90.3));
  propsTwo.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
  propsTwo.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 107));
  log.setProperties(&propsTwo);

  std::vector<std::string> logLineTwo;
  if (log.isCatalystLoggingON()) {
    logLineTwo = log.writeToLogFile();
  }

  checkTestOutputFileExists(log.getLogFileName().c_str());
  std::vector<std::vector<std::string>> lfc = log.readLogFile();

  ASSERT_EQ(lfc.size(), 3);
  EXPECT_EQ(lfc[0], log.getLogFileHeaders());
  EXPECT_EQ(lfc[1], logLineOne);
  EXPECT_EQ(lfc[2], logLineTwo);
}

TEST_F(LoggingTest, WriteLogWithUnsupportedPropertyTypes)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "bar"));
  int ten = 10;
  props.add(Ioss::Property("CATALYST_LOGGING_POINTER_PROP", (void *)&ten));
  std::vector<int> iv = {1, 2, 3};
  props.add(Ioss::Property("CATALYST_LOGGING_VEC_INTEGER_PROP", iv));
  std::vector<double> dv = {1.4, 2.2, -56.3};
  props.add(Ioss::Property("CATALYST_LOGGING_VEC_DOUBLE_PROP", dv));
  log.setProperties(&props);

  std::vector<std::string> logLine;
  if (log.isCatalystLoggingON()) {
    logLine = log.writeToLogFile();
  }

  checkTestOutputFileExists(log.getLogFileName().c_str());
  std::vector<std::vector<std::string>> lfc = log.readLogFile();

  std::vector<std::string> h = {"STRING_PROP"};
  EXPECT_EQ(h, log.getLogFileHeaders());
  ASSERT_EQ(lfc.size(), 2);
  EXPECT_EQ(lfc[0], log.getLogFileHeaders());
  EXPECT_EQ(lfc[1], logLine);
}

TEST_F(LoggingTest, WriteLogWithQuotedCommas)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
  props.add(Ioss::Property("CATALYST_LOGGING_QUOTED_COMMA_PROP", "\"8,7,go,3.2\""));
  props.add(Ioss::Property("CATALYST_LOGGING_QUOTED_ANOTHER_COMMA_PROP",
                           "This address \", 343 Far Lane, ND 8732, RT 3\""));
  log.setProperties(&props);

  std::vector<std::string> logLine;
  if (log.isCatalystLoggingON()) {
    logLine = log.writeToLogFile();
  }

  checkTestOutputFileExists(log.getLogFileName().c_str());
  std::vector<std::vector<std::string>> lfc = log.readLogFile();

  ASSERT_EQ(lfc.size(), 2);
  EXPECT_EQ(lfc[0], log.getLogFileHeaders());
  EXPECT_EQ(lfc[1], logLine);
}

TEST_F(LoggingTest, WriteLogWithQuotedQuotes)
{
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "bar"));
  props.add(Ioss::Property("CATALYST_LOGGING_QUOTED_QUOTE_PROP", "\"\"foo\"\" bar \"1,3,2\""));
  props.add(Ioss::Property("CATALYST_LOGGING_ANOTHER_QUOTED_QUOTE_PROP", "This has \"\"quote\"\""));
  log.setProperties(&props);

  std::vector<std::string> logLine;
  if (log.isCatalystLoggingON()) {
    logLine = log.writeToLogFile();
  }

  checkTestOutputFileExists(log.getLogFileName().c_str());
  std::vector<std::vector<std::string>> lfc = log.readLogFile();

  ASSERT_EQ(lfc.size(), 2);
  EXPECT_EQ(lfc[0], log.getLogFileHeaders());
  EXPECT_EQ(lfc[1], logLine);
}