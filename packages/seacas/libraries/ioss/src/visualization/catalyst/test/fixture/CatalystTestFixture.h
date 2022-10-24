// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef CatalystTestFixture_h
#define CatalystTestFixture_h

#include "IossApplication.h"
#include "vtk_jsoncpp.h"
#include <map>

class CatalystTestFixture
{
public:
  CatalystTestFixture();
  ~CatalystTestFixture();

  void checkPhactoriStringValidParse(const std::string &phactoriSyntax,
                                     const Json::Value &parsedJSONResult);

  void checkPhactoriStringInvalidParse(const std::string &phactoriSyntax);

  void runPhactoriJSONTest(const std::string &jsonFile, const std::string &inputFile);

  void runPhactoriJSONTestTwoGrid(const std::string &jsonFile, const std::string &inputFileA,
                                  const std::string &inputFileB);

  void runParaViewGuiScriptTest(const std::string &pythonScript, const std::string &inputFile);

  void runCatalystLoggingTest(Ioss::PropertyManager *logging_properties,
                              const std::string &jsonFile, const std::string &inputFile);

  static bool isFileExists(const char *fileName);
  static void checkTestOutputFileExists(const char *fileName);
  static void checkTestOutputFileDoesNotExist(const char *fileName);

  Json::Value getDefaultPhactoriJSON();

  Json::Value getDefaultCameraJSON();
  Json::Value getDefaultImageSetJSON();
  Json::Value getDefaultImageSetWithCameraJSON();
  Json::Value getDefaultCameraParallelProjectionJSON();
  Json::Value getDefaultOperationsJSON();

private:
  IossApplication ioapp;
};

#endif
