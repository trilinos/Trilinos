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
  typedef std::vector<std::pair<std::string, int>> VarAndCompCountVec;
  typedef std::vector<std::string>                 StringVec;

  CatalystTestFixture();
  ~CatalystTestFixture();

  void checkPhactoriStringValidParse(const std::string &phactoriSyntax,
                                     const Json::Value &parsedJSONResult);

  void checkPhactoriStringInvalidParse(const std::string &phactoriSyntax);

  void runPhactoriJSONTest(const std::string &jsonFile, const std::string &inputFile);

  void runPhactoriJSONTestTwoGrid(const std::string &jsonFile, const std::string &inputFileA,
                                  const std::string &inputFileB);

  void runPhactoriJSONTestTwoGridTwoPipe(const std::string &jsonFileA,
                                         const std::string &inputFileA,
                                         const std::string &jsonFileB,
                                         const std::string &inputFileB);

  void runParaViewGuiScriptTest(const std::string &pythonScript, const std::string &inputFile);

  void runCatalystLoggingTest(Ioss::PropertyManager *logging_properties,
                              const std::string &jsonFile, const std::string &inputFile);

  void runCatalystMultiBlockMeshTest(const std::string &inputFile);

  void checkMeshOutputVariables(const std::string &inputFile, const VarAndCompCountVec &cellVars,
                                const VarAndCompCountVec &pointVars,
                                const VarAndCompCountVec &globalVars, const std::string &blockPath);

  void checkPartitionedDataSetCollectionStructure(const std::string &inputFile,
                                                  const StringVec &partitions, int numCells,
                                                  const StringVec &searchQueries);

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
