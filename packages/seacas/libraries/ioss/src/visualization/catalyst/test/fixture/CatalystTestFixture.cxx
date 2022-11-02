// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "TestDataDirectoryPath.h"
#include "catch.hpp"
#include <Iovs_Utils.h>
#include <cstdlib>

CatalystTestFixture::CatalystTestFixture() {}

CatalystTestFixture::~CatalystTestFixture() {}

void CatalystTestFixture::runParaViewGuiScriptTest(const std::string &pythonScript,
                                                   const std::string &inputFile)
{
  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.setParaViewExportedScript(td + pythonScript);
  ioapp.addFileName(td + inputFile);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
}

void CatalystTestFixture::runPhactoriJSONTest(const std::string &jsonFile,
                                              const std::string &inputFile)
{

  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.setPhactoriInputJSON(td + jsonFile);
  ioapp.addFileName(td + inputFile);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
}

void CatalystTestFixture::runPhactoriJSONTestTwoGrid(const std::string &jsonFile,
                                                     const std::string &inputFileA,
                                                     const std::string &inputFileB)
{

  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.setPhactoriInputJSON(td + jsonFile);
  ioapp.addFileName(td + inputFileA);
  ioapp.addFileName(td + inputFileB);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
}

void CatalystTestFixture::runCatalystLoggingTest(Ioss::PropertyManager *logging_properties,
                                                 const std::string &    jsonFile,
                                                 const std::string &    inputFile)
{
  ioapp.setAdditionalProperties(logging_properties);
  runPhactoriJSONTest(jsonFile, inputFile);
}

void CatalystTestFixture::checkPhactoriStringValidParse(const std::string &phactoriSyntax,
                                                        const Json::Value &parsedJSONResult)
{

  Iovs::CatalystManagerBase::ParseResult pres;
  Iovs::Utils::getInstance().getCatalystManager().parsePhactoriString(phactoriSyntax, pres);
  REQUIRE(!pres.parseFailed);

  Json::CharReaderBuilder builder{};
  auto                    reader = std::unique_ptr<Json::CharReader>(builder.newCharReader());
  Json::Value             parseRoot{};
  std::string             errors{};

  auto parseWorked = reader->parse(pres.jsonParseResult.c_str(),
                                   pres.jsonParseResult.c_str() + pres.jsonParseResult.length(),
                                   &parseRoot, &errors);

  REQUIRE(parseWorked);
  REQUIRE(parseRoot == parsedJSONResult);
}

void CatalystTestFixture::checkPhactoriStringInvalidParse(const std::string &phactoriSyntax)
{

  Iovs::CatalystManagerBase::ParseResult pres;
  Iovs::Utils::getInstance().getCatalystManager().parsePhactoriString(phactoriSyntax, pres);
  REQUIRE(pres.parseFailed);
}

Json::Value CatalystTestFixture::getDefaultPhactoriJSON()
{
  Json::Value defPhac;
  defPhac["camera blocks"]         = Json::objectValue;
  defPhac["representation blocks"] = Json::objectValue;
  defPhac["operation blocks"]      = Json::objectValue;
  defPhac["imageset blocks"]       = Json::objectValue;
  defPhac["scatter plot blocks"]   = Json::objectValue;
  defPhac["plot over time blocks"] = Json::objectValue;
  defPhac["onoff criteria blocks"] = Json::objectValue;
  defPhac["visual marker blocks"]  = Json::objectValue;
  defPhac["experimental blocks"]   = Json::objectValue;
  return defPhac;
}

Json::Value CatalystTestFixture::getDefaultCameraJSON()
{
  Json::Value camera;
  camera["camera type"] = "camera";

  Json::Value imageSet;
  imageSet["camera"]        = "fooCamera";
  imageSet["image size"][0] = 800;
  imageSet["image size"][1] = 450;

  Json::Value dj                       = getDefaultPhactoriJSON();
  dj["camera blocks"]["fooCamera"]     = camera;
  dj["imageset blocks"]["fooImageset"] = imageSet;
  return dj;
}

Json::Value CatalystTestFixture::getDefaultImageSetJSON()
{
  Json::Value imageSet;
  imageSet["image size"][0] = 800;
  imageSet["image size"][1] = 450;

  Json::Value dj                       = getDefaultPhactoriJSON();
  dj["imageset blocks"]["fooImageset"] = imageSet;
  return dj;
}

Json::Value CatalystTestFixture::getDefaultImageSetWithCameraJSON()
{
  Json::Value dj = getDefaultImageSetJSON();
  Json::Value camera;
  camera["camera type"]                   = "camera";
  camera["look direction"][0]             = -1.0;
  camera["look direction"][1]             = -1.0;
  camera["look direction"][2]             = -1.0;
  dj["camera blocks"]["singletestcamera"] = camera;
  Json::Value &imageset                   = dj["imageset blocks"]["fooImageset"];
  imageset["camera"]                      = "singletestcamera";
  return dj;
}

Json::Value CatalystTestFixture::getDefaultOperationsJSON()
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();
  Json::Value clip;
  clip["type"] = "clip";
  Json::Value plane;
  plane[0]                        = -0.5;
  plane[1]                        = 0.0;
  plane[2]                        = 0.0;
  clip["relative point on plane"] = plane;
  Json::Value normal;
  normal[0]            = 1.0;
  normal[1]            = 0.0;
  normal[2]            = 0.0;
  clip["plane normal"] = normal;
  clip["side to keep"] = "positive";

  Json::Value &imageset              = dj["imageset blocks"]["fooImageset"];
  imageset["operation"]              = "clipNone";
  dj["operation blocks"]["clipNone"] = clip;
  return dj;
}

Json::Value CatalystTestFixture::getDefaultCameraParallelProjectionJSON()
{
  Json::Value camParallel;
  camParallel["camera type"]       = "camera";
  camParallel["projection type"]   = "parallel";
  camParallel["look direction"][0] = -5.0;
  camParallel["look direction"][1] = -1.0;
  camParallel["look direction"][2] = -1.0;

  Json::Value camPerspective;
  camPerspective["camera type"]       = "camera";
  camPerspective["projection type"]   = "perspective";
  camPerspective["look direction"][0] = -5.0;
  camPerspective["look direction"][1] = -1.0;
  camPerspective["look direction"][2] = -1.0;

  Json::Value isParallel;
  isParallel["camera"]         = "parallel_projection_cam1";
  isParallel["image size"][0]  = 800;
  isParallel["image size"][1]  = 450;
  isParallel["image basename"] = "parallel_projection_is1.";

  Json::Value isPerspective;
  isPerspective["camera"]         = "perspective_projection_cam1";
  isPerspective["image size"][0]  = 800;
  isPerspective["image size"][1]  = 450;
  isPerspective["image basename"] = "perspective_projection_is1.";

  Json::Value dj                                      = getDefaultPhactoriJSON();
  dj["camera blocks"]["parallel_projection_cam1"]     = camParallel;
  dj["camera blocks"]["perspective_projection_cam1"]  = camPerspective;
  dj["imageset blocks"]["parallel_projection_is1"]    = isParallel;
  dj["imageset blocks"]["perspective_projection_is1"] = isPerspective;
  return dj;
}

void CatalystTestFixture::checkTestOutputFileExists(const char *fileName)
{
  REQUIRE(isFileExists(fileName) == true);
}

void CatalystTestFixture::checkTestOutputFileDoesNotExist(const char *fileName)
{
  REQUIRE(isFileExists(fileName) == false);
}

bool CatalystTestFixture::isFileExists(const char *fileName)
{
  FILE *fp               = fopen(fileName, "r");
  bool  outputFileExists = false;
  if (fp != NULL) {
    outputFileExists = true;
    fclose(fp);
  }
  return outputFileExists;
  REQUIRE(outputFileExists);
}
