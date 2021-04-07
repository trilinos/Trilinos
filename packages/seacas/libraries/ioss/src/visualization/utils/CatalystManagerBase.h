// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_MANAGER_BASE_H
#define __CATALYST_MANAGER_BASE_H

#include "CatalystCGNSMeshBase.h"
#include "CatalystExodusMeshBase.h"
#include <memory>
#include <string>
#include <vector>

namespace Iovs {

  class CatalystManagerBase
  {

  public:
    CatalystManagerBase(){};
    virtual ~CatalystManagerBase(){};

    struct ParseResult
    {
      std::string jsonParseResult = "";
      bool        parseFailed     = true;
    };

    virtual void parsePhactoriFile(const std::string &filepath, ParseResult &pres) = 0;

    virtual void parsePhactoriString(const std::string &phactori, ParseResult &pres) = 0;

    virtual int getCatalystOutputIDNumber() = 0;

    // Parameters:
    //   cataystPythonFilename - Python file with instructions for Catalyst.
    //   restartTag - if not empty, contains the current restart iteration string, ie s0001
    //   enableLogging - turn on logging in the adapter. Default is off.
    //   debugLevel - enable catalyst debug output 0, 1, 2. Default is 0.
    //   resultsOutputFilename - filename associated with the Ioss results output block.
    //   catalystOutputDirectory - name of the output directory for storing Catalyst output.
    //                              Default is CatalystOutput.
    //   catalystData - string data vector for development and debugging.
    struct CatalystMeshInit
    {
      std::string              catalystPythonFilename;
      std::string              catalystBlockJSON;
      std::string              catalystSeparatorCharacter;
      std::string              catalystInputDeckName;
      std::string              restartTag;
      bool                     enableLogging;
      int                      debugLevel;
      bool                     writeCatalystMeshOneFile;
      std::string              catalystMeshOneFilePrefix;
      bool                     writeCatalystMeshFilePerProc;
      std::string              catalystMeshFilePerProcPrefix;
      std::string              resultsOutputFilename;
      std::string              catalystOutputDirectory;
      std::vector<std::string> catalystData;
      std::string              catalystInputName;
      std::string              catalystMultiInputPipelineName;
      bool                     enableCatalystMultiInputPipeline;
    };

    virtual std::unique_ptr<Iovs_cgns::CatalystCGNSMeshBase>
    createCatalystCGNSMesh(CatalystMeshInit &cmInit) = 0;

    // Parameters:
    //   underscoreVectors - joined vector variable names end in an underscore.
    //   applyDisplacements - a nodal variable named DISPL or displ is applied to
    //                        the mesh node coordinates each time-step.
    struct CatalystExodusMeshInit : CatalystMeshInit
    {
      bool underScoreVectors;
      bool applyDisplacements;
    };

    virtual std::unique_ptr<Iovs_exodus::CatalystExodusMeshBase>
    createCatalystExodusMesh(CatalystExodusMeshInit &cmInit) = 0;
  };

} // namespace Iovs

#endif // __CATALYST_MANAGER_BASE_H
