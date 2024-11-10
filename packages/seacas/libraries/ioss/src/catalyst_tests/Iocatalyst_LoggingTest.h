// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "gtest/gtest.h"
#include <catalyst/Iocatalyst_CatalystLogging.h>

class LoggingTest : public ::testing::Test
{
protected:
  Iocatalyst::CatalystLogging log;
  Ioss::PropertyManager       props;

  ~LoggingTest() { remove(log.getLogFileName().c_str()); }

  void checkTestOutputFileExists(const char *fileName) { EXPECT_TRUE(isFileExists(fileName)); }

  bool isFileExists(const char *fileName)
  {
    FILE *fp               = fopen(fileName, "r");
    bool  outputFileExists = false;
    if (fp != NULL) {
      outputFileExists = true;
      fclose(fp);
    }
    return outputFileExists;
  }
};
