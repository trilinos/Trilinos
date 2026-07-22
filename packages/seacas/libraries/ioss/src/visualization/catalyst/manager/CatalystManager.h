// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_MANAGER_H
#define __CATALYST_MANAGER_H

#include "vtkCPDataDescription.h"
#include <map>
#include <time.h>
#include <visualization/utils/CatalystManagerBase.h>
#include <vtkSmartPointer.h>

class coProcessor;
class vtkDoubleArray;
class vtkCPPythonPipeline;
class vtkCPProcessor;
class vtkDataObject;

namespace Iovs {

  class CatalystMeshWriter;

  class CatalystManager : public CatalystManagerBase
  {

  public:
    const std::string catalystPluginVersion = "3.0.0";
    using CatalystPipelineID                = unsigned int;
    using CatalystInputName                 = std::string;
    using CatalystMultiInputPipelineName    = std::string;

    CatalystManager();
    ~CatalystManager();

    std::string getCatalystPluginVersion();

    std::unique_ptr<Iovs_exodus::CatalystExodusMeshBase>
    createCatalystExodusMesh(CatalystExodusMeshInit &cmInit);

    std::unique_ptr<Iovs_cgns::CatalystCGNSMeshBase>
    createCatalystCGNSMesh(CatalystMeshInit &cmInit);

    int getCatalystOutputIDNumber();

    struct CatalystPipelineInfo
    {
      CatalystPipelineID catalystPipelineID;
      CatalystInputName  catalystInputName;
      std::string        getLogFileName() const
      {
        return catalystInputName + "_" + std::to_string(catalystPipelineID) + "_catalyst.log";
      }
    };

    // Description:
    // Deletes pipeline with name results_output_filename and any associated
    // logging data.
    void DeletePipeline(const CatalystPipelineInfo &cpi);

    // Description:
    // Calls the ParaView Catalyst pipeline to run co-processing for this time iteration.
    void PerformCoProcessing(std::vector<int>           &error_and_warning_codes,
                             std::vector<std::string>   &error_and_warning_messages,
                             const CatalystPipelineInfo &cpi);

    // Description:
    // Sets time data for this ParaView Catalyst co-processing iteration.
    // currentTime is the current Ioss simulation time and timeStep is
    // the current time iteration count.
    void SetTimeData(double currentTime, int timeStep, const CatalystPipelineInfo &cpi);

    // Description:
    // Collects memory usage information from all processors and
    // writes the min, max, and mean to the log file.  Also writes the
    // min, max, and mean of the elapsed time since this method was
    // last called.
    void logMemoryUsageAndTakeTimerReading(const CatalystPipelineInfo &cpi);

    void WriteToLogFile(const CatalystPipelineInfo &cpi);

  private:
    void initCatalystPythonSystemPaths();

    class CatalystPipelineState
    {
    public:
      CatalystPipelineState()
      {
        performCoProcessingCount = 0;
        deletePipelineCount      = 0;
        setTimeDataCount         = 0;
      }

      vtkSmartPointer<vtkCPPythonPipeline> &getPipeline() { return pipeline; }

      vtkSmartPointer<vtkCPDataDescription> &getDataDescription() { return dataDescription; }

      std::shared_ptr<CatalystMeshWriter> &getMeshWriter() { return meshWriter; }

      bool canPerformCoProcessing() { return canDoOperation(performCoProcessingCount); }

      bool canDeletePipeline() { return canDoOperation(deletePipelineCount); }

      bool canSetTimeData() { return canDoOperation(setTimeDataCount); }

    private:
      bool canDoOperation(unsigned int &operationCount)
      {
        bool retVal = true;
        operationCount++;
        if (operationCount == getNumberOfInputs()) {
          operationCount = 0;
        }
        else {
          retVal = false;
        }
        return retVal;
      }

      unsigned int getNumberOfInputs()
      {
        unsigned int retVal = 0;
        if (getDataDescription() != nullptr) {
          retVal = getDataDescription()->GetNumberOfInputDescriptions();
        }
        return retVal;
      }

      unsigned int                          performCoProcessingCount;
      unsigned int                          deletePipelineCount;
      unsigned int                          setTimeDataCount;
      vtkSmartPointer<vtkCPPythonPipeline>  pipeline;
      vtkSmartPointer<vtkCPDataDescription> dataDescription;
      std::shared_ptr<CatalystMeshWriter>   meshWriter;
    };

    typedef std::pair<clock_t, clock_t>            TimerPair;
    typedef std::pair<TimerPair, vtkDoubleArray *> LoggingPair;

    CatalystManager(const CatalystManager &)            = delete;
    CatalystManager &operator=(const CatalystManager &) = delete;

    void               initializeIfNeeded();
    void               finalizeIfNeeded();
    bool               canCoProcess();
    void               incrementOutputCounts();
    bool               writeMeshON(const CatalystPipelineInfo &cpi);
    void               writeMesh(const CatalystPipelineInfo &cpi);
    CatalystPipelineID getCatalystPipelineID(CatalystMeshInit &cmInit);

    void                 initCatalystLogging(const CatalystPipelineInfo &cpi);
    void                 initCatalystPipeline(CatalystMeshInit &cmInit, vtkDataObject *vobj,
                                              const CatalystPipelineInfo &cpi);
    void                 addInputToPipeline(vtkDataObject *vobj, const CatalystPipelineInfo &cpi);
    CatalystPipelineInfo createCatalystPipelineInfo(CatalystMeshInit &cmInit);
    void                 registerMeshInPipeline(CatalystMeshInit &cmInit, vtkDataObject *vobj,
                                                const CatalystPipelineInfo &cpi);

    CatalystPipelineID                                           catalystOutputIDNumber;
    CatalystPipelineID                                           catalystOutputReferenceCount;
    vtkCPProcessor                                              *coProcessor;
    std::map<CatalystPipelineID, CatalystPipelineState>          pipelines;
    std::map<CatalystPipelineID, LoggingPair>                    logging;
    std::map<CatalystMultiInputPipelineName, CatalystPipelineID> multiInputPipelines;
  };

  extern "C" {
  CatalystManagerBase *CreateCatalystManagerInstance();
  }

} // namespace Iovs

#endif /* __CATALYST_MANAGER_H */
