// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_IOVS_UTILS_H
#define IOSS_IOVS_UTILS_H

#include "iovs_export.h"

#include "CatalystManagerBase.h"
#include "Ioss_DBUsage.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_PropertyManager.h"
#include <string>

namespace Iovs {

  class IOVS_EXPORT Utils
  {

  public:
    static Utils &getInstance()
    {
      static Utils instance;
      return instance;
    }

    static bool fileExists(const std::string &filepath);

    std::string getCatalystPythonDriverPath();

    void checkDbUsage(Ioss::DatabaseUsage db_usage);

    struct DatabaseInfo
    {
      std::string                databaseFilename;
      std::string                separatorCharacter;
      const Ioss::ParallelUtils *parallelUtils;
    };

    void createDatabaseOutputFile(const DatabaseInfo &dbinfo);

    std::unique_ptr<Iovs_exodus::CatalystExodusMeshBase>
    createCatalystExodusMesh(const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props);

    std::unique_ptr<Iovs_cgns::CatalystCGNSMeshBase>
    createCatalystCGNSMesh(const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props);

    std::string getDatabaseOutputFilePath(const std::string           &databaseFilename,
                                          const Ioss::PropertyManager &properties);

    void reportCatalystErrorMessages(const std::vector<int>         &error_codes,
                                     const std::vector<std::string> &error_messages, int myRank);

    void writeToCatalystLogFile(const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props);

    CatalystManagerBase &getCatalystManager();

  private:
    Utils();
    ~Utils();
    Utils(const Utils &)            = delete;
    Utils &operator=(const Utils &) = delete;

    CatalystManagerBase *catalystManager = nullptr;

    CatalystManagerBase *createCatalystManagerInstance();
    void                 checkCatalystInterfaceAndPluginVersions();

    void initMeshFromIOSSProps(CatalystManagerBase::CatalystMeshInit &cmInit,
                               const DatabaseInfo &dbinfo, const Ioss::PropertyManager &props);

    std::string getRestartTag(const std::string &databaseFilename);

    void broadCastString(std::string &s, const DatabaseInfo &dbinfo);

    void broadCastStatusCode(bool &statusCode, const DatabaseInfo &dbinfo);

    void loadPluginLibrary();
    void setPythonPathForParaViewPythonZipFile(std::string &paraviewPythonZipFilePath);

    bool getCatalystPluginPath(std::string &catalystPluginPath, std::string &libOSMesaPath);

    std::string getSierraInstallDirectory();

    std::string getCatalystAdapterInstallDirectory();

    void *getDlHandle();

    void *dlHandle          = nullptr;
    void *dlHandleLibOSMesa = nullptr;

#if defined(__APPLE__)
    const char *CATALYST_PLUGIN_DYNAMIC_LIBRARY = "libcatalystioss.dylib";
#else
    const char *CATALYST_PLUGIN_DYNAMIC_LIBRARY = "libcatalystioss.so";
#endif
    const char *CATALYST_PLUGIN_PYTHON_MODULE     = "PhactoriDriver.py";
    const char *CATALYST_PLUGIN_PATH              = "viz/catalyst/install";
    const char *CATALYST_FILE_SUFFIX              = ".dummy.pv.catalyst.e";
    const char *CATALYST_OUTPUT_DIRECTORY         = "CatalystOutput";
    const char *CATALYST_INSTALL_LIB_DIR          = "/lib/";
    const char *CATALYST_INSTALL_PHACTORI_DIR     = "/phactori/";
    const char *CATALYST_IOSS_CATALYST_PLUGIN_DIR = "/current_ioss_catalyst_plugin_version";
    const char *CATALYST_LIB_OSMESA               = "libOSMesa.so";
    const char *CATALYST_LIB_OSMESA_DIR           = "/current_paraview_install/lib/";
    const char *CATALYST_PARAVIEW_PYTHON_ZIP_FILE =
        "/current_paraview_lib_python/site-packages/_paraview.zip";
  };

} // namespace Iovs

#endif
