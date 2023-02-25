// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef Ioss_Application_h
#define Ioss_Application_h

#include <string>
#include <vector>

namespace Ioss {
  class Region;
  class PropertyManager;
} // namespace Ioss

class IossApplication
{
public:
  IossApplication(int argc, char **argv);

  IossApplication();

  ~IossApplication();

  void runApplication();

  int getApplicationExitCode();

  bool printIOSSRegionReportON();
  void setPrintIOSSRegionReport(bool status);

  bool outputCopyOfInputDatabaseON();
  void setOutputCopyOfInputDatabase(bool status);

  bool outputCatalystMeshOneFileON();
  void setOutputCatalystMeshOneFile(bool status);

  bool outputCatalystMeshFilePerProcON();
  void setOutputCatalystMeshFilePerProc(bool status);

  bool forceCGNSOutputON();
  void setForceCGNSOutput(bool status);

  bool forceExodusOutputON();
  void setForceExodusOutput(bool status);

  bool        useIOSSInputDBTypeON();
  std::string getIOSSInputDBType();
  void        setIOSSInputDBType(const std::string &dbType);

  bool        usePhactoriInputScriptON();
  std::string getPhactoriInputScript();
  void        setPhactoriInputScript(const std::string &scriptFilePath);

  bool usePhactoriInputJSONON();

  int         getNumberOfPhactoriInputJSONs();
  std::string getPhactoriInputJSON(int ndx);
  void        addPhactoriInputJSON(const std::string &jsonFilePath);

  bool        useParaViewExportedScriptON();
  std::string getParaViewExportedScript();
  void        setParaViewExportedScript(const std::string &exportedScriptFilePath);

  bool useCatalystStartTimeStepON();
  int  getCatalystStartTimeStep();
  void setCatalystStartTimeStep(int timeStep);

  bool useCatalystStopTimeStepON();
  int  getCatalystStopTimeStep();
  void setCatalystStopTimeStep(int timeStep);

  bool sendMultipleGridsToTheSamePipelineON();
  void setSendMultipleGridsToTheSamePipeline(bool onOffFlag);

  std::string &getFileName(int ndx);
  int          getNumberOfFileNames();
  void         addFileName(const std::string &name);

  void                   setAdditionalProperties(Ioss::PropertyManager *additionalProperties);
  Ioss::PropertyManager *getAdditionalProperties();

private:
  int           getMyRank();
  int           getNumRanks();
  bool          isRankZero();
  bool          isSerial();
  bool          decomposedMeshExists(int ndx);
  int           getNumberOfInputIOSSRegions();
  Ioss::Region *getInputIOSSRegion(int ndx);
  void          copyInputIOSSDatabaseOnRank();
  void          printMessage(const std::string &message);
  void          printErrorMessage(const std::string &message);
  void          printIOSSRegionReportsForRank();
  void          exitApplicationSuccess();
  void          exitApplicationFailure();
  void          exitApplication();

  void        initialize();
  void        SetUpDefaultProperties(Ioss::PropertyManager *outputProperties);
  void        addAdditionalProperties(Ioss::PropertyManager *outputProperties);
  void        callCatalystIOSSDatabaseOnRank();
  void        callCatalystIOSSDatabaseOnRankOneGrid();
  void        callCatalystIOSSDatabaseOnRankMultiGrid(bool sendAllGridsToOnePipeline);
  void        openInputIOSSDatabase(int ndx);
  void        openInputIOSSDatabases();
  void        processCommandLine(int argc, char **argv);
  void        initializeMPI(int argc, char **argv);
  void        initMPIRankAndSize();
  void        finalizeMPI();
  void        printUsageMessage();
  void        checkForOnlyOneCatalystOutputPath();
  void        checkForOnlyOneCatalystOutputType();
  void        getStartStopTimeSteps(int numTimeSteps, int &startTimeStep, int &stopTimeStep);
  std::string getIOSSDatabaseTypeFromFile(int ndx);
  std::string getIOSSDatabaseType(int ndx);
  std::string getCatalystDatabaseType(int ndx);
  std::string getFileSuffix(int ndx);
  std::string getParallelFileName(int ndx);
  std::string getPhactoriDefaultJSON();
  int         myRank;
  int         numRanks;
  bool        useCatalystStartTimeStep;
  int         catalystStartTimeStep;
  bool        useCatalystStopTimeStep;
  int         catalystStopTimeStep;
  bool        printIOSSReport;
  bool        copyDatabase;
  bool        writeCatalystMeshOneFile;
  bool        writeCatalystMeshFilePerProc;
  bool        usePhactoriInputScript;
  bool        usePhactoriInputJSON;
  bool        useParaViewExportedScript;
  bool        forceCGNSOutput;
  bool        forceExodusOutput;
  bool        useIOSSInputDBType;
  bool        hasCommandLineArguments;
  int         applicationExitCode;
  bool        sendMultipleGridsToTheSamePipeline;
  std::string iossInputDBType;
  std::string phactoriInputScriptFilePath;
  std::string phactoriInputJSONFilePath;
  std::string paraViewExportedScriptFilePath;
  std::vector<std::string>    fileName;
  std::vector<std::string>    phctriInptJSONFilePathList;
  std::string                 copyOutputDatabaseName     = "iossDatabaseCopy";
  std::string                 outputCatalystMeshFileName = "iossDatabaseCatalystMesh";
  std::string                 iossReportFileName         = "IossRegionReport";
  const std::string           applicationName            = "ioss2catalyst";
  std::vector<Ioss::Region *> inputIOSSRegion;
  Ioss::PropertyManager      *additionalProperties;

#if defined(__APPLE__)
  const char *CATALYST_PLUGIN_DYNAMIC_LIBRARY = "libcatalystioss.dylib";
#else
  const char *CATALYST_PLUGIN_DYNAMIC_LIBRARY = "libcatalystioss.so";
#endif
};

#endif
