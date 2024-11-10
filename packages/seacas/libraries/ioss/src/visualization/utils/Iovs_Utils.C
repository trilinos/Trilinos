// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_CodeTypes.h"
#include "Ioss_Utils.h"
#include "visualization/utils/Iovs_CatalystLogging.h"
#include "visualization/utils/Iovs_CatalystVersion.h"
#include "visualization/utils/Iovs_Utils.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#if defined(__IOSS_WINDOWS__)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <libgen.h>
#endif

#ifdef IOSS_DLOPEN_ENABLED
#include <dlfcn.h>
#endif

#include <sys/stat.h>

namespace Iovs {
  std::string persistentLdLibraryPathEnvForCatalyst = "";

  Utils::Utils()
  {
    this->dlHandle        = nullptr;
    this->catalystManager = nullptr;
  }

  Utils::~Utils()
  {
    if (this->catalystManager) {
      delete this->catalystManager;
    }
#ifdef IOSS_DLOPEN_ENABLED
    if (this->dlHandle != nullptr) {
      // when loading both libOSMesa.so and the plugin using dlopen,
      // dlclose seg faults. Only call dlclose() if we only made one
      // dlopen call
      if (this->dlHandleLibOSMesa == nullptr) {
        dlclose(this->dlHandle);
      }
    }
#endif
  }

  CatalystManagerBase &Utils::getCatalystManager()
  {
    if (this->catalystManager == nullptr) {
      this->catalystManager = this->createCatalystManagerInstance();
      this->checkCatalystInterfaceAndPluginVersions();
    }
    return *this->catalystManager;
  }

  CatalystManagerBase *Utils::createCatalystManagerInstance()
  {
#ifdef IOSS_DLOPEN_ENABLED
    void *dlh = this->getDlHandle();

    if (!dlh) {
      return nullptr;
    }

    typedef CatalystManagerBase *(*CatalystManagerInstanceFuncType)();

#ifdef __GNUC__
    __extension__
#endif
        CatalystManagerInstanceFuncType mkr = reinterpret_cast<CatalystManagerInstanceFuncType>(
            dlsym(dlh, "CreateCatalystManagerInstance"));
    if (mkr == nullptr) {
      throw std::runtime_error("dlsym call failed to load function "
                               "'CreateCatalystManagerInstance'");
    }
    return (*mkr)();
#else
    return nullptr;
#endif
  }

  void Utils::checkCatalystInterfaceAndPluginVersions()
  {
    CatalystVersion cv;
    std::string     iVer = cv.getIOSSCatalystInterfaceVersion();
    std::string     pVer = this->catalystManager->getCatalystPluginVersion();
    if (!cv.isInterfaceCompatibleWithPlugin(iVer, pVer)) {
      throw std::runtime_error("IOSS Catalyst interface version: " + iVer +
                               ", is not compatible with IOSS Catalyst plugin version: " + pVer);
    }
  }

  std::unique_ptr<Iovs_exodus::CatalystExodusMeshBase>
  Utils::createCatalystExodusMesh(const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props)
  {

    CatalystManagerBase::CatalystExodusMeshInit cmInit;
    cmInit.resultsOutputFilename = dbinfo.databaseFilename;
    this->initMeshFromIOSSProps(cmInit, dbinfo, props);

    cmInit.underScoreVectors = true;
    if (props.exists("CATALYST_UNDERSCORE_VECTORS")) {
      cmInit.underScoreVectors = props.get("CATALYST_UNDERSCORE_VECTORS").get_int();
    }

    cmInit.applyDisplacements = true;
    if (props.exists("CATALYST_APPLY_DISPLACEMENTS")) {
      cmInit.applyDisplacements = props.get("CATALYST_APPLY_DISPLACEMENTS").get_int();
    }

    cmInit.catalystSeparatorCharacter = dbinfo.separatorCharacter;
    cmInit.restartTag                 = this->getRestartTag(dbinfo.databaseFilename);

    return this->getCatalystManager().createCatalystExodusMesh(cmInit);
  }

  std::unique_ptr<Iovs_cgns::CatalystCGNSMeshBase>
  Utils::createCatalystCGNSMesh(const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props)
  {

    CatalystManagerBase::CatalystMeshInit cmInit;
    cmInit.resultsOutputFilename = dbinfo.databaseFilename;
    this->initMeshFromIOSSProps(cmInit, dbinfo, props);
    cmInit.catalystSeparatorCharacter = dbinfo.separatorCharacter;
    cmInit.restartTag                 = this->getRestartTag(dbinfo.databaseFilename);

    return this->getCatalystManager().createCatalystCGNSMesh(cmInit);
  }

  void Utils::initMeshFromIOSSProps(CatalystManagerBase::CatalystMeshInit &cmInit,
                                    const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props)
  {

    if (props.exists("CATALYST_BLOCK_PARSE_JSON_STRING")) {
      cmInit.catalystBlockJSON = props.get("CATALYST_BLOCK_PARSE_JSON_STRING").get_string();
    }
    else if (props.exists("PHACTORI_JSON_SCRIPT")) {
      bool        readOkay             = false;
      std::string phactoriJSONFilePath = props.get("PHACTORI_JSON_SCRIPT").get_string();
      if (dbinfo.parallelUtils->parallel_rank() == 0) {
        std::ifstream f(phactoriJSONFilePath);
        if (f) {
          std::ostringstream ss;
          ss << f.rdbuf();
          cmInit.catalystBlockJSON = ss.str();
          readOkay                 = true;
        }
      }
      this->broadCastStatusCode(readOkay, dbinfo);
      if (!readOkay) {
        std::ostringstream errmsg;
        errmsg << "Unable to read input file: " << phactoriJSONFilePath << "\n";
        IOSS_ERROR(errmsg);
      }
      else {
        this->broadCastString(cmInit.catalystBlockJSON, dbinfo);
      }
    }
    else if (props.exists("PHACTORI_INPUT_SYNTAX_SCRIPT")) {
      std::string phactoriFilePath = props.get("PHACTORI_INPUT_SYNTAX_SCRIPT").get_string();
      CatalystManagerBase::ParseResult pres;
      if (dbinfo.parallelUtils->parallel_rank() == 0) {
        this->getCatalystManager().parsePhactoriFile(phactoriFilePath, pres);
      }
      this->broadCastStatusCode(pres.parseFailed, dbinfo);
      if (pres.parseFailed) {
        std::ostringstream errmsg;
        errmsg << "Unable to parse input file: " << phactoriFilePath << "\n";
        IOSS_ERROR(errmsg);
      }
      else {
        this->broadCastString(pres.jsonParseResult, dbinfo);
        cmInit.catalystBlockJSON = pres.jsonParseResult;
      }
    }

    cmInit.writeCatalystMeshOneFile = false;
    if (props.exists("WRITE_CATALYST_MESH_ONE_FILE_WITH_PREFIX")) {
      cmInit.writeCatalystMeshOneFile = true;
      cmInit.catalystMeshOneFilePrefix =
          props.get("WRITE_CATALYST_MESH_ONE_FILE_WITH_PREFIX").get_string();
    }

    cmInit.writeCatalystMeshFilePerProc = false;
    if (props.exists("WRITE_CATALYST_MESH_FILE_PER_PROC_WITH_PREFIX")) {
      cmInit.writeCatalystMeshFilePerProc = true;
      cmInit.catalystMeshFilePerProcPrefix =
          props.get("WRITE_CATALYST_MESH_FILE_PER_PROC_WITH_PREFIX").get_string();
    }

    if (props.exists("CATALYST_SCRIPT")) {
      cmInit.catalystPythonFilename = props.get("CATALYST_SCRIPT").get_string();
    }
    else {
      cmInit.catalystPythonFilename = this->getCatalystPythonDriverPath();
    }

    if (props.exists("CATALYST_SCRIPT_EXTRA_FILE")) {
      cmInit.catalystData.push_back(props.get("CATALYST_SCRIPT_EXTRA_FILE").get_string());
    }

    if (props.exists("CATALYST_BLOCK_PARSE_INPUT_DECK_NAME")) {
      cmInit.catalystInputDeckName = props.get("CATALYST_BLOCK_PARSE_INPUT_DECK_NAME").get_string();
    }

    cmInit.enableLogging = false;
    if (props.exists("CATALYST_ENABLE_LOGGING")) {
      cmInit.enableLogging = props.get("CATALYST_ENABLE_LOGGING").get_int();
    }

    cmInit.debugLevel = 0;
    if (props.exists("CATALYST_DEBUG_LEVEL")) {
      cmInit.enableLogging = props.get("CATALYST_DEBUG_LEVEL").get_int();
    }

    cmInit.catalystOutputDirectory = CATALYST_OUTPUT_DIRECTORY;
    if (props.exists("CATALYST_OUTPUT_DIRECTORY")) {
      cmInit.catalystOutputDirectory = props.get("CATALYST_OUTPUT_DIRECTORY").get_string();
    }

    cmInit.catalystInputName = "input";
    if (props.exists("CATALYST_INPUT_NAME")) {
      cmInit.catalystInputName = props.get("CATALYST_INPUT_NAME").get_string();
    }

    cmInit.enableCatalystMultiInputPipeline = false;
    if (props.exists("CATALYST_MULTI_INPUT_PIPELINE_NAME")) {
      cmInit.enableCatalystMultiInputPipeline = true;
      cmInit.catalystMultiInputPipelineName =
          props.get("CATALYST_MULTI_INPUT_PIPELINE_NAME").get_string();
    }
  }

  void Utils::writeToCatalystLogFile(const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props)
  {
    if (dbinfo.parallelUtils->parallel_rank() == 0) {
      CatalystLogging catLog = CatalystLogging();
      catLog.setProperties(&props);
      catLog.writeToLogFile();
    }
    dbinfo.parallelUtils->barrier();
  }

  std::string Utils::getRestartTag(const std::string &databaseFilename)
  {
    std::string            restartTag;
    std::string::size_type pos = databaseFilename.rfind(".e-s");
    if (pos != std::string::npos) {
      if (pos + 3 <= databaseFilename.length()) {
        restartTag = databaseFilename.substr(pos + 3, 5);
      }
    }
    return restartTag;
  }

  bool Utils::fileExists(const std::string &filepath)
  {
    struct stat buffer
    {
    };
    return (stat(filepath.c_str(), &buffer) == 0);
  }

  std::string Utils::getDatabaseOutputFilePath(const std::string           &databaseFilename,
                                               const Ioss::PropertyManager &properties)
  {

    int numberOfCatalystBlocks = this->getCatalystManager().getCatalystOutputIDNumber();
    if (!properties.exists("CATALYST_OUTPUT_DIRECTORY")) {
      std::ostringstream s;
      s << databaseFilename << "." << numberOfCatalystBlocks << CATALYST_FILE_SUFFIX;
      return std::string(CATALYST_OUTPUT_DIRECTORY) + "/" + s.str();
    }
    else {
      return databaseFilename;
    }
  }

  void *Utils::getDlHandle()
  {
    if (this->dlHandle == nullptr) {
      loadPluginLibrary();
    }
    return this->dlHandle;
  }

  void Utils::loadPluginLibrary()
  {

    std::string pluginLibraryPath;
    bool        callDlopenLibOSMesa{};
    std::string libOSMesaPath;

    callDlopenLibOSMesa = this->getCatalystPluginPath(pluginLibraryPath, libOSMesaPath);

#ifdef IOSS_DLOPEN_ENABLED
    if (callDlopenLibOSMesa) {
      this->dlHandleLibOSMesa = dlopen(libOSMesaPath.c_str(), RTLD_NOW | RTLD_GLOBAL);
      if (this->dlHandleLibOSMesa == nullptr) {
        throw std::runtime_error(dlerror());
      }
    }
    else {
      this->dlHandleLibOSMesa = nullptr;
    }
    this->dlHandle = dlopen(pluginLibraryPath.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (this->dlHandle == nullptr) {
      throw std::runtime_error(dlerror());
    }
#else
    this->dlHandle          = nullptr;
    this->dlHandleLibOSMesa = nullptr;
#endif
  }

  bool Utils::getCatalystPluginPath(std::string &catalystPluginPath, std::string &libOSMesaPath)
  {
    bool callDlopenLibOSMesa = false;

    if (getenv("CATALYST_PLUGIN") != nullptr) {
      catalystPluginPath  = getenv("CATALYST_PLUGIN");
      callDlopenLibOSMesa = false;
      libOSMesaPath       = CATALYST_LIB_OSMESA;
      return callDlopenLibOSMesa;
    }

    std::string catalystInsDir = this->getCatalystAdapterInstallDirectory();

    if (!catalystInsDir.empty()) {
      catalystPluginPath =
          catalystInsDir + std::string(CATALYST_INSTALL_LIB_DIR) + CATALYST_PLUGIN_DYNAMIC_LIBRARY;
      callDlopenLibOSMesa = false;
      libOSMesaPath       = CATALYST_LIB_OSMESA;
      return callDlopenLibOSMesa;
    }

    catalystPluginPath =
        this->getSierraInstallDirectory() + "/" + std::string(CATALYST_PLUGIN_DYNAMIC_LIBRARY);
    callDlopenLibOSMesa = true;
    libOSMesaPath =
        this->getSierraInstallDirectory() + CATALYST_LIB_OSMESA_DIR + CATALYST_LIB_OSMESA;

    // to avoid depending on name python-3.7, python-3.8, etc., we
    // expect a link to the python directory inside the sierra catalyst
    // install along CATALYST_PARAVIEW_PYTHON_ZIP_FILE
    std::string paraviewPythonZipFile =
        this->getSierraInstallDirectory() + CATALYST_PARAVIEW_PYTHON_ZIP_FILE;
    setPythonPathForParaViewPythonZipFile(paraviewPythonZipFile);
    return callDlopenLibOSMesa;
  }

  void Utils::setPythonPathForParaViewPythonZipFile(std::string &paraviewPythonZipFilePath)
  {
    const char *existingPythonpath = getenv("PYTHONPATH");
    if (existingPythonpath == nullptr) {
      persistentLdLibraryPathEnvForCatalyst = paraviewPythonZipFilePath;
    }
    else {
      persistentLdLibraryPathEnvForCatalyst = paraviewPythonZipFilePath + ":" + existingPythonpath;
    }
#if defined(__IOSS_WINDOWS__)
    SetEnvironmentVariableA("PYTHONPATH", persistentLdLibraryPathEnvForCatalyst.c_str());
#else
    setenv("PYTHONPATH", persistentLdLibraryPathEnvForCatalyst.c_str(), 1);
#endif
  }

  std::string Utils::getCatalystPythonDriverPath()
  {

    std::string catalystInsDir = this->getCatalystAdapterInstallDirectory();

    if (!catalystInsDir.empty()) {
      return catalystInsDir + std::string(CATALYST_INSTALL_PHACTORI_DIR) +
             CATALYST_PLUGIN_PYTHON_MODULE;
    }

    return this->getSierraInstallDirectory() + "/" + std::string(CATALYST_PLUGIN_PYTHON_MODULE);
  }

  std::string Utils::getSierraInstallDirectory()
  {

    std::string sierraInsDir;
    if (getenv("SIERRA_INSTALL_DIR") != nullptr) {
      sierraInsDir = getenv("SIERRA_INSTALL_DIR");
    }
    else {
      std::ostringstream errmsg;
      errmsg << "Environment variable SIERRA_INSTALL_DIR not set.\n"
             << " Unable to find ParaView catalyst dynamic library.\n";
      IOSS_ERROR(errmsg);
    }

    std::string sierraSystem;
    if (getenv("SIERRA_SYSTEM") != nullptr) {
      sierraSystem = getenv("SIERRA_SYSTEM");
    }
    else {
      std::ostringstream errmsg;
      errmsg << "Environment variable SIERRA_SYSTEM not set.\n"
             << " Unable to find ParaView catalyst dynamic library.\n";
      IOSS_ERROR(errmsg);
    }

#if defined(__IOSS_WINDOWS__)
    char *cbuf = _fullpath(nullptr, sierraInsDir.c_str(), _MAX_PATH);
#else
    char *cbuf = realpath(sierraInsDir.c_str(), nullptr);
#endif
    std::string sierraInsPath = cbuf;
    free(cbuf);

    if (!fileExists(sierraInsPath)) {
      std::ostringstream errmsg;
      errmsg << "SIERRA_INSTALL_DIR directory does not exist.\n"
             << "Directory path: " << sierraInsPath << "\n"
             << " Unable to find ParaView catalyst dynamic library.\n";
      IOSS_ERROR(errmsg);
    }

#if defined(__IOSS_WINDOWS__)
    {
      std::ostringstream errmsg;
      errmsg << "This code is not yet supported on windows...\n";
      IOSS_ERROR(errmsg);
    }
#else
    char *cbase = strdup(sierraInsPath.c_str());
    char *cdir  = strdup(sierraInsPath.c_str());
    char *bname = basename(cbase);
    char *dname = dirname(cdir);

    while (strcmp(dname, "/") != 0 && strcmp(dname, ".") != 0 && strcmp(bname, "sierra") != 0) {
      bname = basename(dname);
      dname = dirname(dname);
    }

    if (strcmp(bname, "sierra") == 0) {
      sierraInsPath = dname;
    }

    free(cbase);
    free(cdir);
#endif

    return sierraInsPath + "/" + CATALYST_PLUGIN_PATH + "/" + sierraSystem +
           CATALYST_IOSS_CATALYST_PLUGIN_DIR;
  }

  std::string Utils::getCatalystAdapterInstallDirectory()
  {
    std::string catalystInsDir = "";
    if (getenv("CATALYST_ADAPTER_INSTALL_DIR") != nullptr) {
      catalystInsDir = getenv("CATALYST_ADAPTER_INSTALL_DIR");

      if (!fileExists(catalystInsDir)) {
        std::ostringstream errmsg;
        errmsg << "CATALYST_ADAPTER_INSTALL_DIR directory does\n"
               << "not exist. Directory path: " << catalystInsDir << "\n"
               << "Unable to find ParaView catalyst dynamic library.\n";
        IOSS_ERROR(errmsg);
      }
    }
    return catalystInsDir;
  }

  void Utils::checkDbUsage(Ioss::DatabaseUsage db_usage)
  {
    std::ostringstream errmsg;
    if (db_usage == Ioss::WRITE_HEARTBEAT) {
      errmsg << "ParaView catalyst database type cannot be"
             << " used in a HEARTBEAT block.\n";
      IOSS_ERROR(errmsg);
    }
    else if (db_usage == Ioss::WRITE_HISTORY) {
      errmsg << "ParaView catalyst database type cannot be"
             << " used in a HISTORY block.\n";
      IOSS_ERROR(errmsg);
    }
    else if (db_usage == Ioss::READ_MODEL || db_usage == Ioss::READ_RESTART) {

      errmsg << "ParaView catalyst database type cannot be"
             << " used to read a model.\n";
      IOSS_ERROR(errmsg);
    }
  }

  void Utils::createDatabaseOutputFile(const DatabaseInfo &dbinfo)
  {
    std::ostringstream errmsg;
    if (dbinfo.parallelUtils->parallel_rank() == 0) {
      if (!Utils::fileExists(dbinfo.databaseFilename)) {
        std::ofstream output_file;
        output_file.open(dbinfo.databaseFilename.c_str(), std::ios::out | std::ios::trunc);

        if (!output_file) {
          errmsg << "Unable to create output file: " << dbinfo.databaseFilename << ".\n";
          IOSS_ERROR(errmsg);
        }
        output_file.close();
      }
    }
    dbinfo.parallelUtils->barrier();
  }

  void Utils::reportCatalystErrorMessages(const std::vector<int>         &error_codes,
                                          const std::vector<std::string> &error_messages,
                                          int                             myRank)
  {

    if (!error_codes.empty() && !error_messages.empty() &&
        error_codes.size() == error_messages.size()) {
      for (unsigned int i = 0; i < error_codes.size(); i++) {
        if (error_codes[i] > 0) {
          Ioss::WarnOut() << "\n\n** ParaView Catalyst Plugin Warning Message Severity Level "
                          << error_codes[i] << ", On Processor " << myRank << " **\n\n";
          Ioss::WarnOut() << error_messages[i];
        }
        else {
          std::ostringstream errmsg;
          errmsg << "\n\n** ParaView Catalyst Plugin Error Message Severity Level "
                 << error_codes[i] << ", On Processor " << myRank << " **\n\n"
                 << error_messages[i];
          IOSS_ERROR(errmsg);
        }
      }
    }
  }

  void Utils::broadCastString(IOSS_MAYBE_UNUSED std::string        &s,
                              IOSS_MAYBE_UNUSED const DatabaseInfo &dbinfo)
  {
    IOSS_PAR_UNUSED(s);
    IOSS_PAR_UNUSED(dbinfo);
#ifdef SEACAS_HAVE_MPI
    int size = s.size();
    dbinfo.parallelUtils->broadcast(size);
    if (dbinfo.parallelUtils->parallel_rank() != 0) {
      s.resize(size);
    }
    dbinfo.parallelUtils->broadcast(s);
#endif
  }

  void Utils::broadCastStatusCode(IOSS_MAYBE_UNUSED bool               &statusCode,
                                  IOSS_MAYBE_UNUSED const DatabaseInfo &dbinfo)
  {
    IOSS_PAR_UNUSED(statusCode);
    IOSS_PAR_UNUSED(dbinfo);
#ifdef SEACAS_HAVE_MPI

    int code = statusCode;
    dbinfo.parallelUtils->broadcast(code);
    statusCode = code;
#endif
  }
} // namespace Iovs
