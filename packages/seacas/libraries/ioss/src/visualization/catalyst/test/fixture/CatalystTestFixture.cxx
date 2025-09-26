// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "TestDataDirectoryPath.h"
#include "vtkAbstractArray.h"
#include "vtkCellData.h"
#include "vtkDataAssembly.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkNew.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPartitionedDataSetCollection.h"
#include "vtkPointData.h"
#include "vtkXMLPartitionedDataSetCollectionReader.h"
#include <Iovs_Utils.h>
#include <cstdlib>
#include <catch2/catch_test_macros.hpp>

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

void CatalystTestFixture::checkMeshOutputVariables(const std::string        &inputFile,
                                                   const VarAndCompCountVec &cellVars,
                                                   const VarAndCompCountVec &pointVars,
                                                   const VarAndCompCountVec &globalVars,
                                                   const std::string        &blockPath)
{
  vtkNew<vtkXMLPartitionedDataSetCollectionReader> vpdcr;

  vpdcr->SetFileName(inputFile.c_str());
  vpdcr->Update();
  vtkPartitionedDataSetCollection *vpdc =
      vtkPartitionedDataSetCollection::SafeDownCast(vpdcr->GetOutput());
  REQUIRE(vpdc->GetDataAssembly() != nullptr);

  auto assembly   = vpdc->GetDataAssembly();
  auto childNodes = assembly->GetChildNodes(assembly->GetFirstNodeByPath(blockPath.c_str()));
  bool foundBlockThatHasAllVars = false;
  for (size_t i = 0; i < childNodes.size(); i++) {
    auto dsi = assembly->GetDataSetIndices(childNodes[i]);
    for (size_t j = 0; j < dsi.size(); j++) {
      auto pds = vpdc->GetPartitionedDataSet(dsi[j]);
      for (unsigned int k = 0; k < pds->GetNumberOfPartitions(); k++) {
        vtkDataSet *ds = pds->GetPartition(k);
        if (ds == nullptr) {
          continue;
        }

        auto hasAllVars = [](vtkFieldData *fd, const VarAndCompCountVec &vars) {
          for (auto vv : vars) {
            vtkAbstractArray *ar = fd->GetAbstractArray(vv.first.c_str());
            if (ar == nullptr) {
              return false;
            }
            if (ar->GetNumberOfComponents() != vv.second) {
              return false;
            }
          }
          return true;
        };

        if (hasAllVars(ds->GetCellData(), cellVars) && hasAllVars(ds->GetPointData(), pointVars) &&
            hasAllVars(ds->GetFieldData(), globalVars)) {
          foundBlockThatHasAllVars = true;
        }
      }
    }
  }

  REQUIRE(foundBlockThatHasAllVars);
}

void CatalystTestFixture::checkPartitionedDataSetCollectionStructure(const std::string &inputFile,
                                                                     const StringVec   &partitions,
                                                                     int                numCells,
                                                                     const StringVec &searchQueries)
{

  vtkNew<vtkXMLPartitionedDataSetCollectionReader> vpdcr;

  vpdcr->SetFileName(inputFile.c_str());
  vpdcr->Update();
  vtkPartitionedDataSetCollection *vpdc =
      vtkPartitionedDataSetCollection::SafeDownCast(vpdcr->GetOutput());
  REQUIRE(vpdc->GetDataAssembly() != nullptr);
  REQUIRE(vpdc->GetDataAssembly()->GetRootNodeName() == std::string("IOSS"));
  REQUIRE(vpdc->GetNumberOfPartitionedDataSets() == partitions.size());
  int numCellsCount = 0;
  for (unsigned int i = 0; i < vpdc->GetNumberOfPartitionedDataSets(); i++) {
    REQUIRE(vpdc->HasMetaData(i));
    REQUIRE(vpdc->GetMetaData(i)->Get(vtkCompositeDataSet::NAME()) == partitions[i]);
    auto pds = vpdc->GetPartitionedDataSet(i);
    REQUIRE(pds != nullptr);
    auto num_parts = pds->GetNumberOfPartitions();
    for (unsigned int j = 0; j < num_parts; j++) {
      int  partNumCells = pds->GetPartition(j)->GetNumberOfCells();
      REQUIRE(partNumCells > 0);
      numCellsCount += partNumCells;
    }
  }
  REQUIRE(numCells == numCellsCount);
  auto ids = vpdc->GetDataAssembly()->SelectNodes(searchQueries);
  REQUIRE(ids.size() == searchQueries.size());
}

void CatalystTestFixture::runCatalystMultiBlockMeshTest(const std::string &inputFile)
{
  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.addFileName(td + inputFile);
  ioapp.setOutputCatalystMeshOneFile(true);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
}

void CatalystTestFixture::runPhactoriJSONTest(const std::string &jsonFile,
                                              const std::string &inputFile)
{

  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.addPhactoriInputJSON(td + jsonFile);
  ioapp.addFileName(td + inputFile);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
}

void CatalystTestFixture::runPhactoriJSONTestTwoGrid(const std::string &jsonFile,
                                                     const std::string &inputFileA,
                                                     const std::string &inputFileB)
{

  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.addPhactoriInputJSON(td + jsonFile);
  ioapp.addFileName(td + inputFileA);
  ioapp.addFileName(td + inputFileB);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
}

void CatalystTestFixture::runPhactoriJSONTestTwoGridTwoPipe(const std::string &jsonFileA,
                                                            const std::string &inputFileA,
                                                            const std::string &jsonFileB,
                                                            const std::string &inputFileB)
{

  ioapp.setSendMultipleGridsToTheSamePipeline(false);
  std::string td = std::string(TEST_DATA_DIRECTORY_PATH);
  ioapp.addPhactoriInputJSON(td + jsonFileA);
  ioapp.addFileName(td + inputFileA);
  ioapp.addPhactoriInputJSON(td + jsonFileB);
  ioapp.addFileName(td + inputFileB);
  ioapp.runApplication();
  REQUIRE(ioapp.getApplicationExitCode() == EXIT_SUCCESS);
  ioapp.setSendMultipleGridsToTheSamePipeline(true);
}

void CatalystTestFixture::runCatalystLoggingTest(Ioss::PropertyManager *logging_properties,
                                                 const std::string     &jsonFile,
                                                 const std::string     &inputFile)
{
  ioapp.setAdditionalProperties(logging_properties);
  runPhactoriJSONTest(jsonFile, inputFile);
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
}
