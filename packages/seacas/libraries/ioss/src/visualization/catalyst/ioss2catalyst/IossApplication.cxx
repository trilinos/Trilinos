// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "IossApplication.h"
#include "CatalystPluginPaths.h"
#include "Ionit_Initializer.h"
#include "IossRegionReport.h"
#include "Ioss_CopyDatabase.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_FileInfo.h"
#include "Ioss_IOFactory.h"
#include "Ioss_MeshCopyOptions.h"
#include "Ioss_Region.h"
#include "Ioss_Utils.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <random>
#include <unistd.h>
std::random_device myrd9;        // only used once to initialise (seed) engine
std::mt19937       rng(myrd9()); // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> myuni9(0, 1000000000); // guaranteed unbiased

IossApplication::IossApplication()
{
  initialize();
  setenv("CATALYST_ADAPTER_INSTALL_DIR", CATALYST_PLUGIN_BUILD_DIR, false);
  std::string pluginLibPath =
      std::string(CATALYST_PLUGIN_BUILD_DIR) + std::string("/") + CATALYST_PLUGIN_DYNAMIC_LIBRARY;
  setenv("CATALYST_PLUGIN", pluginLibPath.c_str(), false);
}

IossApplication::IossApplication(int argc, char **argv)
{
  initializeMPI(argc, argv);
  initialize();
  processCommandLine(argc, argv);
  setenv("CATALYST_ADAPTER_INSTALL_DIR", CATALYST_PLUGIN_INSTALL_DIR, false);
}

void IossApplication::initialize()
{
  myRank                             = 0;
  numRanks                           = 1;
  printIOSSReport                    = false;
  copyDatabase                       = false;
  writeCatalystMeshOneFile           = false;
  writeCatalystMeshFilePerProc       = false;
  usePhactoriInputScript             = false;
  usePhactoriInputJSON               = false;
  useParaViewExportedScript          = false;
  phactoriInputScriptFilePath        = "";
  phactoriInputJSONFilePath          = "";
  paraViewExportedScriptFilePath     = "";
  catalystStartTimeStep              = 0;
  useCatalystStartTimeStep           = false;
  catalystStopTimeStep               = 0;
  useCatalystStopTimeStep            = false;
  forceCGNSOutput                    = false;
  forceExodusOutput                  = false;
  useIOSSInputDBType                 = false;
  hasCommandLineArguments            = false;
  applicationExitCode                = EXIT_SUCCESS;
  additionalProperties               = nullptr;
  sendMultipleGridsToTheSamePipeline = true;
  initMPIRankAndSize();
}

IossApplication::~IossApplication()
{
  int ii;
  int numRegions = inputIOSSRegion.size();
  for (ii = 0; ii < numRegions; ii++) {
    if (inputIOSSRegion[ii] != nullptr) {
      delete inputIOSSRegion[ii];
      inputIOSSRegion[ii] = nullptr;
    }
  }
}

void IossApplication::runApplication()
{
  checkForOnlyOneCatalystOutputPath();
  checkForOnlyOneCatalystOutputType();
  Ioss::Init::Initializer io;
  openInputIOSSDatabases();

  if (printIOSSRegionReportON()) {
    printIOSSRegionReportsForRank();
  }

  if (outputCopyOfInputDatabaseON()) {
    copyInputIOSSDatabaseOnRank();
  }

  callCatalystIOSSDatabaseOnRank();
  sync();

  exitApplicationSuccess();
}

int IossApplication::getApplicationExitCode() { return applicationExitCode; }

int IossApplication::getMyRank() { return myRank; }

int IossApplication::getNumRanks() { return numRanks; }

bool IossApplication::isRankZero() { return myRank == 0; }

bool IossApplication::isSerial() { return numRanks == 1; }

void IossApplication::initializeMPI(int argc, char **argv) { MPI_Init(&argc, &argv); }

void IossApplication::initMPIRankAndSize()
{
  Ioss::ParallelUtils pu{};
  myRank   = pu.parallel_rank();
  numRanks = pu.parallel_size();
}

void IossApplication::finalizeMPI() { MPI_Finalize(); }

void IossApplication::processCommandLine(int argc, char **argv)
{
  hasCommandLineArguments = true;

  int   c;
  char *cvalue = NULL;
  while ((c = getopt(argc, argv, "a:b:cd:e:hi:mnp:rs:")) != -1) {
    char *cvalue = nullptr;
    switch (c) {
    case 'a':
      cvalue = optarg;
      setCatalystStartTimeStep(atoi(cvalue));
      break;
    case 'b':
      cvalue = optarg;
      setCatalystStopTimeStep(atoi(cvalue));
      break;
    case 'c': copyDatabase = true; break;
    case 'd':
      cvalue = optarg;
      setForceExodusOutput(true);
      setIOSSInputDBType(cvalue);
      break;
    case 'e':
      cvalue = optarg;
      setForceCGNSOutput(true);
      setIOSSInputDBType(cvalue);
      break;
    case 'h':
      printUsageMessage();
      exitApplicationSuccess();
      break;
    case 'i':
      cvalue = optarg;
      addPhactoriInputJSON(cvalue);
      break;
    case 'm': setOutputCatalystMeshOneFile(true); break;
    case 'n': setOutputCatalystMeshFilePerProc(true); break;
    case 'p':
      cvalue = optarg;
      setPhactoriInputScript(cvalue);
      break;
    case 'r': printIOSSReport = true; break;
    case 's':
      cvalue = optarg;
      setParaViewExportedScript(cvalue);
      break;
    case '?':
      printErrorMessage("Unknown command line option -" + std::string(1, c) + "\n");
      printUsageMessage();
      exitApplicationFailure();
      break;
    default:
      printUsageMessage();
      exitApplicationFailure();
      break;
    }
  }

  if (optind == argc) {
    printErrorMessage("No input filename given on command line");
    printUsageMessage();
    exitApplicationFailure();
  }
  else if (optind != argc - 1) {
    printErrorMessage("Expected one argument for input filename. " +
                      std::string("Got multiple arguments."));
    printUsageMessage();
    exitApplicationFailure();
  }
  else {
    // fileName = argv[optind];
    addFileName(argv[optind]);
  }
}

void IossApplication::getStartStopTimeSteps(int numTimeSteps, int &startTimeStep, int &stopTimeStep)
{

  if (useCatalystStartTimeStepON()) {
    startTimeStep = getCatalystStartTimeStep();
    if (startTimeStep < 1 || startTimeStep > numTimeSteps) {
      printErrorMessage("Start time-step out of range");
      exitApplicationFailure();
    }
  }
  else {
    startTimeStep = 1;
  }

  if (useCatalystStopTimeStepON()) {
    stopTimeStep = getCatalystStopTimeStep();
    if (stopTimeStep < 1 || stopTimeStep > numTimeSteps) {
      printErrorMessage("Stop time-step out of range");
      exitApplicationFailure();
    }
  }
  else {
    stopTimeStep = numTimeSteps;
  }

  if (startTimeStep > stopTimeStep) {
    printErrorMessage("Start time-step > stop time-step.");
    exitApplicationFailure();
  }
}

void IossApplication::checkForOnlyOneCatalystOutputPath()
{

  int numTimesCatalystCalled = 0;
  if (usePhactoriInputScriptON()) {
    numTimesCatalystCalled++;
  }
  if (usePhactoriInputJSONON()) {
    numTimesCatalystCalled++;
  }
  if (outputCatalystMeshOneFileON()) {
    numTimesCatalystCalled++;
  }
  if (outputCatalystMeshFilePerProcON()) {
    numTimesCatalystCalled++;
  }
  if (useParaViewExportedScriptON()) {
    numTimesCatalystCalled++;
  }

  if (numTimesCatalystCalled > 1) {
    printErrorMessage("Catalyst called with more than one option.");
    printUsageMessage();
    exitApplicationFailure();
  }

  if (numTimesCatalystCalled == 1 && (printIOSSRegionReportON() || outputCopyOfInputDatabaseON())) {
    printErrorMessage("Catalyst called with report output.");
    printUsageMessage();
    exitApplicationFailure();
  }
}

void IossApplication::checkForOnlyOneCatalystOutputType()
{
  if (forceCGNSOutputON() && forceExodusOutputON()) {
    printErrorMessage("Both CGNS and Exodus Catalyst output requested.");
    printUsageMessage();
    exitApplicationFailure();
  }
}

int IossApplication::getNumberOfFileNames() { return fileName.size(); }

std::string &IossApplication::getFileName(int ndx)
{
  if (ndx >= fileName.size()) {
    printErrorMessage("getFileName called with ndx too large.");
    printUsageMessage();
    exitApplicationFailure();
  }
  return fileName[ndx];
}

void IossApplication::addFileName(const std::string &name) { fileName.push_back(name); }

bool IossApplication::printIOSSRegionReportON() { return printIOSSReport; }

void IossApplication::setPrintIOSSRegionReport(bool status) { printIOSSReport = status; }

bool IossApplication::outputCopyOfInputDatabaseON() { return copyDatabase; }

void IossApplication::setOutputCopyOfInputDatabase(bool status) { copyDatabase = status; }

bool IossApplication::outputCatalystMeshOneFileON() { return writeCatalystMeshOneFile; }

void IossApplication::setOutputCatalystMeshOneFile(bool status)
{
  writeCatalystMeshOneFile = status;
}

bool IossApplication::outputCatalystMeshFilePerProcON() { return writeCatalystMeshFilePerProc; }

void IossApplication::setOutputCatalystMeshFilePerProc(bool status)
{
  writeCatalystMeshFilePerProc = status;
}

bool IossApplication::usePhactoriInputScriptON() { return usePhactoriInputScript; }

std::string IossApplication::getPhactoriInputScript() { return phactoriInputScriptFilePath; }

void IossApplication::setPhactoriInputScript(const std::string &scriptFilePath)
{
  phactoriInputScriptFilePath = scriptFilePath;
  usePhactoriInputScript      = true;
}

bool IossApplication::usePhactoriInputJSONON() { return usePhactoriInputJSON; }

int IossApplication::getNumberOfPhactoriInputJSONs() { return phctriInptJSONFilePathList.size(); }

std::string IossApplication::getPhactoriInputJSON(int ndx)
{
  if (ndx >= phctriInptJSONFilePathList.size()) {
    printErrorMessage("phctriInptJSONFilePathList called with ndx too large.");
    printUsageMessage();
    exitApplicationFailure();
  }
  return phctriInptJSONFilePathList[ndx];
}

void IossApplication::addPhactoriInputJSON(const std::string &jsonFilePath)
{
  phctriInptJSONFilePathList.push_back(jsonFilePath);
  usePhactoriInputJSON = true;
}

bool IossApplication::useParaViewExportedScriptON() { return useParaViewExportedScript; }

std::string IossApplication::getParaViewExportedScript() { return paraViewExportedScriptFilePath; }

void IossApplication::setParaViewExportedScript(const std::string &exportedScriptFilePath)
{
  paraViewExportedScriptFilePath = exportedScriptFilePath;
  useParaViewExportedScript      = true;
}

bool IossApplication::useCatalystStartTimeStepON() { return useCatalystStartTimeStep; }

int IossApplication::getCatalystStartTimeStep() { return catalystStartTimeStep; }

void IossApplication::setCatalystStartTimeStep(int timeStep)
{
  catalystStartTimeStep    = timeStep;
  useCatalystStartTimeStep = true;
}

bool IossApplication::useCatalystStopTimeStepON() { return useCatalystStopTimeStep; }

int IossApplication::getCatalystStopTimeStep() { return catalystStopTimeStep; }

void IossApplication::setCatalystStopTimeStep(int timeStep)
{
  catalystStopTimeStep    = timeStep;
  useCatalystStopTimeStep = true;
}

bool IossApplication::forceCGNSOutputON() { return forceCGNSOutput; }

void IossApplication::setForceCGNSOutput(bool status) { forceCGNSOutput = status; }

bool IossApplication::forceExodusOutputON() { return forceExodusOutput; }

void IossApplication::setForceExodusOutput(bool status) { forceExodusOutput = status; }

bool IossApplication::useIOSSInputDBTypeON() { return useIOSSInputDBType; }

std::string IossApplication::getIOSSInputDBType() { return iossInputDBType; }

bool IossApplication::sendMultipleGridsToTheSamePipelineON()
{
  return sendMultipleGridsToTheSamePipeline;
}

void IossApplication::setSendMultipleGridsToTheSamePipeline(bool onOffFlag)
{
  sendMultipleGridsToTheSamePipeline = onOffFlag;
}

void IossApplication::setIOSSInputDBType(const std::string &dbType)
{
  iossInputDBType    = dbType;
  useIOSSInputDBType = true;
}

void IossApplication::printUsageMessage()
{
  std::string um = "\nUSAGE\n\n" + applicationName;

  um += " [OPTIONS] [OUTPUT OPTIONS] <FILE>\n\n";

  um += "DESCRIPTION\n\n";
  um += "Read input file(s) and write to ParaView Catalyst using IOSS";

  um += "\n\nOPTIONS\n\n";

  um += "-a <n> call Catalyst starting at time-step n (starts at 1).\n\n";

  um += "-b <n> call Catalyst stopping at time-step n.\n\n";

  um += "[-d | -e] <db_type> Force IOSS input database type to db_type\n";
  um += "   -d for catalyst_exodus output or -e for catalyst_cgns output\n\n";

  um += "-h print usage message and exit program\n\n";

  um += "OUTPUT OPTIONS\n\n";
  um += " [-cr | -m | -n | -i <file> | -p <file> | -s <file>]\n\n";

  um += "-c copy input file(s) to one file per processor with output \n";
  um += "   filename(s) prefix " + copyOutputDatabaseName + "\n\n";

  um += "-i <file> run Catalyst with Phactori input JSON given in <file>.";
  um += "\n\n";

  um += "-m output Catalyst mesh representation of input file(s)\n";
  um += "   each time-step to a single file for all processors with output\n";
  um += "   filename " + outputCatalystMeshFileName + "_time_<n>.vtm";
  um += "\n\n";

  um += "-n output Catalyst mesh representation of input file(s)\n";
  um += "   each time-step to a file for each processor with output\n";
  um += "   filename " + outputCatalystMeshFileName;
  um += "_proc_<p>_time_<n>.vtm\n\n";

  um += "-p <file> run Catalyst with Phactori input command syntax given\n";
  um += "   in <file>\n\n";

  um += "-r print IOSS region report for input file(s) to one file\n";
  um += "   per processor with output filename(s) prefix ";
  um += iossReportFileName + "\n\n";

  um += "-s <file> run Catalyst with a ParaView exported Python script\n";
  um += "   given in <file>";

  um += "\n\nEXAMPLES\n\n";

  um += "mpiexec -np 1 " + applicationName + " -h\n";
  um += "    Print usage message and exit program.";
  um += "\n\n";

  um += "mpiexec -np 4 " + applicationName + " -c -r <FILE>\n";
  um += "    Output copy of input mesh and IOSS region report.";
  um += "\n\n";

  um += "mpiexec -np 4 " + applicationName + " <FILE>\n";
  um += "    Run Catalyst with default Phactori JSON script to produce\n";
  um += "    eight axis aligned external camera images of the input mesh.";
  um += "\n\n";

  um += "mpiexec -np 4 " + applicationName + " -m <FILE>\n";
  um += "    Output mesh representation for Catalyst.\n\n";

  um += "mpiexec -np 4 " + applicationName + " -i file.json <FILE>\n";
  um += "    Run Catalyst with Phactori JSON input from file.json.\n\n";

  um += "mpiexec -np 4 " + applicationName + " -p file.txt <FILE>\n";
  um += "    Run Catalyst with Phactori command syntax input from file.txt.";
  um += "\n\n";

  um += "mpiexec -np 4 " + applicationName + " -s file.py <FILE>\n";
  um += "    Run Catalyst with ParaView exported Python script in file.py.";

  um += "\n\nFILE\n\n";
  um += "Exodus or CGNS input file name for a single file\n\n";
  um += "   (file_name.ex2 or file_name.cgns)\n\n";

  um += "Exodus or CGNS file name prefix for multiple files\n\n";
  um += "   (file_name.cgns for file_name.cgns.2.0, file_name.cgns.2.1)\n";
  um += "   (file_name.ex2 for file_name.ex2.2.0, file_name.ex2.2.1)\n\n";
  printMessage(um);
}

void IossApplication::exitApplicationSuccess()
{
  applicationExitCode = EXIT_SUCCESS;
  if (hasCommandLineArguments) {
    exitApplication();
  }
}

void IossApplication::exitApplicationFailure()
{
  applicationExitCode = EXIT_FAILURE;
  if (hasCommandLineArguments) {
    exitApplication();
  }
}

void IossApplication::exitApplication()
{
  finalizeMPI();
  std::exit(getApplicationExitCode());
}

void IossApplication::printMessage(const std::string &message)
{
  if (isRankZero()) {
    std::cout << message;
  }
}

void IossApplication::printErrorMessage(const std::string &message)
{
  if (isRankZero()) {
    std::cerr << "\nERROR: " << message << "\n";
  }
}

void IossApplication::printIOSSRegionReportsForRank()
{
  std::string fn = iossReportFileName + "." + std::to_string(getNumRanks()) + "." +
                   std::to_string(getMyRank()) + ".txt";
  std::ofstream ofs(fn, std::ofstream::out);
  Ioss::Region *region;

  int ii;
  int numRegions = getNumberOfInputIOSSRegions();

  for (ii = 0; ii < numRegions; ii++) {
    region = getInputIOSSRegion(ii);

    ofs << "begin region report for region with index " << std::to_string(ii) << "\n";
    ofs << ioss_region_report::region_report(*region) << "\n";
    auto state_count = region->get_property("state_count").get_int();
    for (auto state = 1; state <= state_count; ++state) {
      region->begin_state(state);
      ofs << ioss_region_report::region_report(*region) << "\n";
      region->end_state(state);
    }
    ofs << "end region report for region with index " << std::to_string(ii) << "\n";
  }

  ofs.close();
}

std::string IossApplication::getParallelFileName(int ndx)
{
  if (isSerial()) {
    return getFileName(ndx);
  }

  std::stringstream nameStream;
  const int         numZeroes = std::ceil(log10(double(getNumRanks())));

  nameStream << getFileName(ndx) << "." << getNumRanks() << "." << std::setfill('0')
             << std::setw(numZeroes) << getMyRank();

  return nameStream.str();
}

bool IossApplication::decomposedMeshExists(int ndx)
{
  int status = 0;

  if (isRankZero()) {
    std::string   parallelFilename = getParallelFileName(ndx);
    std::ifstream fstream(parallelFilename);
    if (fstream.good()) {
      status = 1;
    }
    if (fstream.is_open()) {
      fstream.close();
    }
  }
  Ioss::ParallelUtils pu{};
  pu.broadcast(status);
  return status == 1;
}

int IossApplication::getNumberOfInputIOSSRegions() { return inputIOSSRegion.size(); }

Ioss::Region *IossApplication::getInputIOSSRegion(int ndx)
{
  if (ndx >= inputIOSSRegion.size()) {
    printErrorMessage("getInputIOSSRegion called with ndx too large.");
    printUsageMessage();
  }
  return inputIOSSRegion[ndx];
}

void IossApplication::openInputIOSSDatabase(int ndx)
{
  Ioss::PropertyManager inputProperties;
  if (decomposedMeshExists(ndx)) {
    inputProperties.add(Ioss::Property("DECOMPOSITION_METHOD", "external"));
  }
  else {
    inputProperties.add(Ioss::Property("DECOMPOSITION_METHOD", "rib"));
  }

  Ioss::DatabaseIO *dbi =
      Ioss::IOFactory::create(getIOSSDatabaseType(ndx), getFileName(ndx), Ioss::READ_RESTART,
                              Ioss::ParallelUtils::comm_world(), inputProperties);
  if (dbi == nullptr || !dbi->ok(true)) {
    printErrorMessage("Unable to open input file(s) " + getFileName(ndx));
    exitApplicationFailure();
  }
  inputIOSSRegion.push_back(new Ioss::Region(dbi));
}

void IossApplication::openInputIOSSDatabases()
{
  int numInputIossDatabases = getNumberOfFileNames();
  int ii;
  for (ii = 0; ii < numInputIossDatabases; ii++) {
    openInputIOSSDatabase(ii);
  }
}

std::string IossApplication::getPhactoriDefaultJSON()
{
  char const *phactoriDefaultJSON = R"pd({
        "camera blocks":{},
        "representation blocks": {},
        "operation blocks": {},
        "imageset blocks": {},
        "scatter plot blocks": {},
        "plot over time blocks": {},
        "marker blocks": {}
    })pd";
  return phactoriDefaultJSON;
}

void IossApplication::copyInputIOSSDatabaseOnRank()
{
  std::string           fn = copyOutputDatabaseName + "." + getFileSuffix(0);
  Ioss::PropertyManager outputProperties;
  outputProperties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
  Ioss::DatabaseIO *dbo =
      Ioss::IOFactory::create(getIOSSDatabaseType(0), fn, Ioss::WRITE_RESULTS,
                              Ioss::ParallelUtils::comm_world(), outputProperties);
  if (dbo == nullptr || !dbo->ok(true)) {
    printErrorMessage("Unable to open output file(s) " + fn);
    exitApplicationFailure();
  }

  Ioss::Region *inputRegion  = getInputIOSSRegion(0);
  Ioss::Region *outputRegion = new Ioss::Region(dbo, inputRegion->name());

  auto                  state_count = inputRegion->get_property("state_count").get_int();
  double                min_time    = inputRegion->get_state_time(1);
  double                max_time    = inputRegion->get_state_time(state_count);
  Ioss::MeshCopyOptions copyOptions;
  copyOptions.data_storage_type = 1;
  copyOptions.minimum_time      = min_time;
  copyOptions.maximum_time      = max_time;
  Ioss::copy_database(*inputRegion, *outputRegion, copyOptions);

  delete outputRegion;
}

void IossApplication::callCatalystIOSSDatabaseOnRank()
{
  if (getNumberOfInputIOSSRegions() < 2) {
    callCatalystIOSSDatabaseOnRankOneGrid();
  }
  else {
    if (this->sendMultipleGridsToTheSamePipelineON()) {
      callCatalystIOSSDatabaseOnRankMultiGrid(true);
    }
    else {
      callCatalystIOSSDatabaseOnRankMultiGrid(false);
    }
  }
}

const char *gGridInputNames[5] = {"input", "inputB", "inputC", "inputD", "inputE"};

void IossApplication::callCatalystIOSSDatabaseOnRankMultiGrid(bool sendAllGridsToOnePipeline)
{
  // create all the ioss database instances
  // create the output region for each grid (one on each database)
  // for each timestep
  //  call copy_database for each timestep for each database

  // assuming we can't share the property manager for now
  int                                  numInputRegions = getNumberOfInputIOSSRegions();
  std::vector<Ioss::PropertyManager *> outputProps;
  std::vector<Ioss::DatabaseIO *>      outputDbs;
  std::vector<Ioss::Region *>          outputRegions;
  int                                  ii;

  printMessage("IossApplication::callCatalystIOSSDatabaseOnRankMultiGrid entered\n");
  std::cerr << "IossApplication::callCatalystIOSSDatabaseOnRankMultiGrid entered2\n";
  // create a Ioss::DatabaseIO instance for each grid
  char mytmpnm[256];
  int  imyRandomInt = myuni9(rng);
  snprintf(mytmpnm, 256, "catalyst_%d", imyRandomInt);

  for (ii = 0; ii < numInputRegions; ii++) {
    Ioss::PropertyManager *newProps = new (Ioss::PropertyManager);
    SetUpDefaultProperties(newProps);
    if (sendAllGridsToOnePipeline) {
      newProps->add(Ioss::Property("CATALYST_MULTI_INPUT_PIPELINE_NAME", "multipipe1"));
    }
    newProps->add(Ioss::Property("CATALYST_INPUT_NAME", gGridInputNames[ii % 5]));
    if (getNumberOfPhactoriInputJSONs() > 1) {
      printMessage("getNumberOfPhactoriInputJSONs() > 1\n");
      printMessage(getPhactoriInputJSON(ii) + "\n");
      newProps->add(Ioss::Property("PHACTORI_JSON_SCRIPT", getPhactoriInputJSON(ii)));
    }
    outputProps.push_back(newProps);
    Ioss::DatabaseIO *newDbo =
        Ioss::IOFactory::create(getCatalystDatabaseType(ii), mytmpnm, Ioss::WRITE_RESULTS,
                                Ioss::ParallelUtils::comm_world(), *newProps);
    outputDbs.push_back(newDbo);
  }
  // create an output Ioss::Region instance for each grid
  for (ii = 0; ii < numInputRegions; ii++) {
    Ioss::Region *inputRegion  = getInputIOSSRegion(ii);
    Ioss::Region *outputRegion = new Ioss::Region(outputDbs[ii], inputRegion->name());
    outputRegions.push_back(outputRegion);
  }

  // use the minimum number of steps from all regions as the overall number of steps
  auto state_count = getInputIOSSRegion(0)->get_property("state_count").get_int();
  for (ii = 1; ii < numInputRegions; ii++) {
    int oneRegionNumStates = getInputIOSSRegion(ii)->get_property("state_count").get_int();
    if (oneRegionNumStates < state_count) {
      state_count = oneRegionNumStates;
    }
  }
  int startTimeStep;
  int stopTimeStep;
  getStartStopTimeSteps(state_count, startTimeStep, stopTimeStep);

  // right now, copy_database can't be called repeatedly with different
  // timesteps because it tries to create the Node block (and other stuff)
  // each time. So we are only doing one timestep for testing at this point.
  // startTimeStep = stopTimeStep;

  // for each timestep, call copy_database (one timestep) for each output
  // database
  int currentTimeStep;
  for (currentTimeStep = startTimeStep; currentTimeStep <= stopTimeStep; currentTimeStep++) {
    for (ii = 0; ii < numInputRegions; ii++) {
      Ioss::Region *inputRegion               = getInputIOSSRegion(ii);
      Ioss::Region *outputRegion              = outputRegions[ii];
      int           thisRegionCurrentTimeStep = currentTimeStep;
      int thisRegionStateCount = getInputIOSSRegion(ii)->get_property("state_count").get_int();
      if (thisRegionCurrentTimeStep > thisRegionStateCount) {
        thisRegionCurrentTimeStep = thisRegionStateCount;
      }
      double min_time = inputRegion->get_state_time(thisRegionCurrentTimeStep);
      double max_time = inputRegion->get_state_time(thisRegionCurrentTimeStep);

      Ioss::MeshCopyOptions copyOptions;
      copyOptions.data_storage_type = 1;
      copyOptions.minimum_time      = min_time;
      copyOptions.maximum_time      = max_time;
      // Ioss::copy_database(*inputRegion, *outputRegion, copyOptions);
      bool defineFlag;
      if (currentTimeStep == startTimeStep) {
        defineFlag = true;
      }
      else {
        defineFlag = false;
      }
      printf("region, ts1, ts2, time: %d, %d, %d, %f\n", ii, currentTimeStep,
             thisRegionCurrentTimeStep, (float)max_time);
      printf("defineFlag: %d\n", (int)defineFlag);
      copyOptions.define_geometry = defineFlag;
      Ioss::copy_database(*inputRegion, *outputRegion, copyOptions);
    }
  }
  for (ii = 0; ii < numInputRegions; ii++) {
    Ioss::Region *outputRegion = outputRegions[ii];
    delete outputRegion;
  }
}

void IossApplication::SetUpDefaultProperties(Ioss::PropertyManager *outputProperties)
{
  if (usePhactoriInputScriptON()) {
    outputProperties->add(Ioss::Property("PHACTORI_INPUT_SYNTAX_SCRIPT", getPhactoriInputScript()));
  }
  else if (usePhactoriInputJSONON()) {
    printMessage("SetUpDefaultProperties, usePhactoriInputJSONON true\n");
    if (getNumberOfPhactoriInputJSONs() == 1) {
      printMessage("1 entry: " + getPhactoriInputJSON(0) + "\n");
      outputProperties->add(Ioss::Property("PHACTORI_JSON_SCRIPT", getPhactoriInputJSON(0)));
    }
  }
  else if (useParaViewExportedScriptON()) {
    outputProperties->add(Ioss::Property("CATALYST_SCRIPT", getParaViewExportedScript()));
  }
  else if (outputCatalystMeshOneFileON()) {
    outputProperties->add(
        Ioss::Property("WRITE_CATALYST_MESH_ONE_FILE_WITH_PREFIX", outputCatalystMeshFileName));
  }
  else if (outputCatalystMeshFilePerProcON()) {
    outputProperties->add(Ioss::Property("WRITE_CATALYST_MESH_FILE_PER_PROC_WITH_PREFIX",
                                         outputCatalystMeshFileName));
  }
  else {
    outputProperties->add(
        Ioss::Property("CATALYST_BLOCK_PARSE_JSON_STRING", getPhactoriDefaultJSON()));
  }

  outputProperties->add(Ioss::Property("CATALYST_BLOCK_PARSE_INPUT_DECK_NAME", applicationName));
  addAdditionalProperties(outputProperties);
}

void IossApplication::callCatalystIOSSDatabaseOnRankOneGrid()
{
  Ioss::PropertyManager outputProperties;

  SetUpDefaultProperties(&outputProperties);

  Ioss::DatabaseIO *dbo =
      Ioss::IOFactory::create(getCatalystDatabaseType(0), "catalyst", Ioss::WRITE_RESULTS,
                              Ioss::ParallelUtils::comm_world(), outputProperties);
  if (dbo == nullptr || !dbo->ok(true)) {
    printErrorMessage("Unable to open catalyst database");
    exitApplicationFailure();
  }

  Ioss::Region *inputRegion  = getInputIOSSRegion(0);
  Ioss::Region *outputRegion = new Ioss::Region(dbo, inputRegion->name());

  auto state_count = inputRegion->get_property("state_count").get_int();

  int startTimeStep;
  int stopTimeStep;
  getStartStopTimeSteps(state_count, startTimeStep, stopTimeStep);

  double min_time = inputRegion->get_state_time(startTimeStep);
  double max_time = inputRegion->get_state_time(stopTimeStep);

  Ioss::MeshCopyOptions copyOptions;
  copyOptions.data_storage_type = 1;
  copyOptions.minimum_time      = min_time;
  copyOptions.maximum_time      = max_time;
  Ioss::copy_database(*inputRegion, *outputRegion, copyOptions);

  delete outputRegion;
}

std::string IossApplication::getIOSSDatabaseType(int ndx)
{
  std::string retVal = getIOSSDatabaseTypeFromFile(ndx);
  if (useIOSSInputDBTypeON()) {
    retVal = getIOSSInputDBType();
  }
  return retVal;
}

std::string IossApplication::getIOSSDatabaseTypeFromFile(int ndx)
{
  if (ndx >= fileName.size()) {
    printErrorMessage("getIOSSDatabaseTypeFromFile called with ndx too large.");
    printUsageMessage();
    exitApplicationFailure();
  }
  Ioss::FileInfo file(fileName[ndx]);
  auto           extension = file.extension();
  if (extension == "e" || extension == "g" || extension == "gen" || extension == "exo") {
    return "exodus";
  }
  else if (extension == "cgns") {
    return "cgns";
  }
  else {
    return "exodus";
  }
}

std::string IossApplication::getFileSuffix(int ndx)
{
  std::string dbType = getIOSSDatabaseType(ndx);
  if (dbType == "exodus") {
    return "ex2";
  }
  else if (dbType == "cgns") {
    return "cgns";
  }
  else {
    return dbType;
  }
}

std::string IossApplication::getCatalystDatabaseType(int ndx)
{
  std::string retVal = "catalyst_exodus";
  if (!forceExodusOutputON()) {
    if (getIOSSDatabaseType(ndx) == "cgns") {
      retVal = "catalyst_cgns";
    }
    if (forceCGNSOutputON()) {
      retVal = "catalyst_cgns";
    }
  }
  return retVal;
}

void IossApplication::setAdditionalProperties(Ioss::PropertyManager *additionalProperties)
{
  this->additionalProperties = additionalProperties;
}

Ioss::PropertyManager *IossApplication::getAdditionalProperties() { return additionalProperties; }

void IossApplication::addAdditionalProperties(Ioss::PropertyManager *outputProperties)
{
  if (additionalProperties) {

    Ioss::NameList nlist;
    additionalProperties->describe(&nlist);
    for (auto &name : nlist) {
      outputProperties->add(additionalProperties->get(name));
    }
  }
}
