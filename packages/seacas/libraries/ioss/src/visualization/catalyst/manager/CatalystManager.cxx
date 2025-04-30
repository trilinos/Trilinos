// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystManager.h"
#include "../cgns/CatalystCGNSMesh.h"
#include "../exodus/CatalystExodusMesh.h"
#include "CatalystMeshWriter.h"
#include "CatalystPythonPaths.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonPipeline.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkMPIController.h"
#include "vtkMultiProcessController.h"
#include "vtkProcessModule.h"
#include "vtkPython.h"
#include "vtkStringArray.h"
#include "vtkPartitionedDataSetCollection.h"
#include <fstream>
#include <sstream>
#include <vtksys/SystemInformation.hxx>

namespace Iovs {

  CatalystManager::CatalystManager()
  {
    coProcessor                  = nullptr;
    catalystOutputIDNumber       = 0;
    catalystOutputReferenceCount = 0;
#ifdef USE_PYTHON_C_API
    initCatalystPythonSystemPaths();
#endif
  }

  CatalystManager::~CatalystManager() {}

  void CatalystManager::initCatalystPythonSystemPaths()
  {
    Py_Initialize();

    PyObject *sys  = PyImport_ImportModule("sys");
    PyObject *path = PyObject_GetAttrString(sys, "path");

    PyList_Insert(path, 0, PyUnicode_FromString(CATALYST_PYTHON_PARAVIEW_ZIP));
    PyList_Insert(path, 0, PyUnicode_FromString(CATALYST_PYTHON_VTK_ZIP));

    Py_DECREF(sys);
    Py_DECREF(path);
  }

  std::string CatalystManager::getCatalystPluginVersion() { return catalystPluginVersion; }

  int CatalystManager::getCatalystOutputIDNumber() { return catalystOutputIDNumber; }

  void CatalystManager::initializeIfNeeded()
  {
    if (!canCoProcess()) {
      coProcessor = vtkCPProcessor::New();
      coProcessor->Initialize();
      catalystOutputReferenceCount = 0;
    }
  }

  void CatalystManager::finalizeIfNeeded()
  {
    if (canCoProcess()) {
      coProcessor->Delete();
      coProcessor = nullptr;
    }
  }

  bool CatalystManager::canCoProcess() { return coProcessor != nullptr; }

  void CatalystManager::incrementOutputCounts()
  {
    catalystOutputIDNumber++;
    catalystOutputReferenceCount++;
  }

  std::unique_ptr<Iovs_exodus::CatalystExodusMeshBase>
  CatalystManager::createCatalystExodusMesh(CatalystExodusMeshInit &cmInit)
  {

    initializeIfNeeded();
    CatalystPipelineInfo cpi = createCatalystPipelineInfo(cmInit);

    Iovs_exodus::CatalystExodusMesh *cem = new Iovs_exodus::CatalystExodusMesh(this, cpi);
    cem->SetUnderscoreVectors(cmInit.underScoreVectors);
    cem->SetApplyDisplacements(cmInit.applyDisplacements);

    registerMeshInPipeline(cmInit, cem->getPartitionedDataSetCollection(), cpi);

    return std::unique_ptr<Iovs_exodus::CatalystExodusMeshBase>(
        dynamic_cast<Iovs_exodus::CatalystExodusMeshBase *>(cem));
  }

  std::unique_ptr<Iovs_cgns::CatalystCGNSMeshBase>
  CatalystManager::createCatalystCGNSMesh(CatalystMeshInit &cmInit)
  {

    initializeIfNeeded();
    CatalystPipelineInfo cpi = createCatalystPipelineInfo(cmInit);

    Iovs_cgns::CatalystCGNSMesh *cgm = new Iovs_cgns::CatalystCGNSMesh(this, cpi);

    registerMeshInPipeline(cmInit, cgm->getPartitionedDataSetCollection(), cpi);

    return std::unique_ptr<Iovs_cgns::CatalystCGNSMeshBase>(
        dynamic_cast<Iovs_cgns::CatalystCGNSMeshBase *>(cgm));
  }

  CatalystManager::CatalystPipelineInfo
  CatalystManager::createCatalystPipelineInfo(CatalystMeshInit &cmInit)
  {

    CatalystPipelineInfo cpi;
    cpi.catalystPipelineID = getCatalystPipelineID(cmInit);
    cpi.catalystInputName  = cmInit.catalystInputName;
    return cpi;
  }

  CatalystManager::CatalystPipelineID
  CatalystManager::getCatalystPipelineID(CatalystMeshInit &cmInit)
  {

    CatalystManager::CatalystPipelineID id                      = catalystOutputIDNumber;
    bool                                doIncrementOutputCounts = true;

    if (cmInit.enableCatalystMultiInputPipeline) {
      std::string pn = cmInit.catalystMultiInputPipelineName;
      if (multiInputPipelines.find(pn) == multiInputPipelines.end()) {
        multiInputPipelines[pn] = catalystOutputIDNumber;
      }
      else {
        doIncrementOutputCounts = false;
      }
      id = multiInputPipelines[pn];
    }

    if (doIncrementOutputCounts) {
      incrementOutputCounts();
    }
    return id;
  }

  void CatalystManager::registerMeshInPipeline(CatalystMeshInit &cmInit, vtkDataObject *vobj,
                                               const CatalystPipelineInfo &cpi)
  {

    if (pipelines.find(cpi.catalystPipelineID) == pipelines.end()) {
      initCatalystPipeline(cmInit, vobj, cpi);
      if (cmInit.enableLogging) {
        initCatalystLogging(cpi);
      }
    }
    else {
      addInputToPipeline(vobj, cpi);
    }
  }

  void CatalystManager::initCatalystLogging(const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id          = cpi.catalystPipelineID;
    std::string        logFileName = cpi.getLogFileName();

    TimerPair       tp = std::make_pair(clock(), clock());
    vtkDoubleArray *da = vtkDoubleArray::New();
    da->SetNumberOfComponents(3);
    LoggingPair lp    = std::make_pair(tp, da);
    this->logging[id] = lp;

    vtkProcessModule *pm   = vtkProcessModule::GetProcessModule();
    vtkMPIController *mpic = vtkMPIController::SafeDownCast(pm->GetGlobalController());
    if (mpic && mpic->GetNumberOfProcesses() > 1) {
      if (mpic->GetLocalProcessId() == 0) {
        std::ofstream logfile;
        logfile.open(logFileName, ios::out | ios::trunc);
        logfile << "# ELAPSED TIME (S)"
                << ",PROC MEM USED - MIN (KiB)"
                << ",PROC MEM USED - MAX (KiB)"
                << ",PROC MEM USED - AVG (KiB)"
                << ",HOST MEM USED - MIN (KiB)"
                << ",HOST MEM USED - MAX (KiB)"
                << ",HOST MEM USED - AVG (KiB)"
                << ",TIME SINCE LAST LOG - MIN (S)"
                << ",TIME SINCE LAST LOG - MAX (S)"
                << ",TIME SINCE LAST LOG - AVG (S)"
                << "\n";
        logfile.close();
      }
    }
    else {
      std::ofstream logfile;
      logfile.open(logFileName, ios::out | ios::trunc);
      logfile << "# ELAPSED TIME (S)"
              << ",PROC MEM USED (KiB)"
              << ",HOST MEM USED (KiB)"
              << ",TIME SINCE LAST LOG (S)"
              << "\n";
      logfile.close();
    }
  }

  void CatalystManager::initCatalystPipeline(CatalystMeshInit &cmInit, vtkDataObject *vobj,
                                             const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    pipelines[id].getPipeline() =
        vtkCPPythonPipeline::CreateAndInitializePipeline(cmInit.catalystPythonFilename.c_str());
    if (pipelines[id].getPipeline() == nullptr) {
      std::cerr << "Unable to initialize ParaView Catalyst with python script "
                << cmInit.catalystPythonFilename << std::endl;
      return;
    }

    pipelines[id].getDataDescription() = vtkSmartPointer<vtkCPDataDescription>::New();
    pipelines[id].getDataDescription()->AddInput(cpi.catalystInputName.c_str());
    pipelines[id]
        .getDataDescription()
        ->GetInputDescriptionByName(cpi.catalystInputName.c_str())
        ->SetGrid(vobj);

    pipelines[id].getMeshWriter() = std::make_shared<CatalystMeshWriter>();
    if (cmInit.writeCatalystMeshOneFile) {
      pipelines[id].getMeshWriter()->setOutputCatalystMeshOneFilePrefix(
          cmInit.catalystMeshOneFilePrefix);
    }
    if (cmInit.writeCatalystMeshFilePerProc) {
      pipelines[id].getMeshWriter()->setOutputCatalystMeshFilePerProcPrefix(
          cmInit.catalystMeshFilePerProcPrefix);
    }

    vtkFieldData   *fd = vtkFieldData::New();
    vtkStringArray *sa = vtkStringArray::New();
    sa->SetName("catalyst_sierra_data");
    vtkIntArray *ec = vtkIntArray::New();
    ec->SetName("catalyst_sierra_error_codes");
    vtkStringArray *em = vtkStringArray::New();
    em->SetName("catalyst_sierra_error_messages");
    sa->InsertNextValue(cmInit.catalystBlockJSON);
    sa->InsertNextValue(cmInit.catalystSeparatorCharacter);
    sa->InsertNextValue(cmInit.catalystInputDeckName);
    sa->InsertNextValue(cmInit.restartTag);
    if (cmInit.enableLogging) {
      sa->InsertNextValue("True");
    }
    else {
      sa->InsertNextValue("");
    }
    std::stringstream ss;
    ss << cmInit.debugLevel;
    sa->InsertNextValue(ss.str().c_str());
    ss.clear();
    sa->InsertNextValue(cmInit.resultsOutputFilename);
    sa->InsertNextValue(cmInit.catalystOutputDirectory);

    for (size_t i = 0; i < cmInit.catalystData.size(); i++) {
      sa->InsertNextValue(cmInit.catalystData[i]);
    }

    fd->AddArray(sa);
    fd->AddArray(ec);
    fd->AddArray(em);
    pipelines[id].getDataDescription()->SetUserData(fd);

    fd->Delete();
    sa->Delete();
    ec->Delete();
    em->Delete();
  }

  void CatalystManager::addInputToPipeline(vtkDataObject *vobj, const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    pipelines[id].getDataDescription()->AddInput(cpi.catalystInputName.c_str());
    pipelines[id]
        .getDataDescription()
        ->GetInputDescriptionByName(cpi.catalystInputName.c_str())
        ->SetGrid(vobj);
  }

  void CatalystManager::DeletePipeline(const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    if (pipelines.find(id) != pipelines.end()) {
      if (!pipelines[id].canDeletePipeline()) {
        return;
      }
      pipelines.erase(id);
    }

    if (logging.find(id) != logging.end()) {
      logging[id].second->Delete();
      logging.erase(id);
    }

    catalystOutputReferenceCount--;
    finalizeIfNeeded();
  }

  void CatalystManager::PerformCoProcessing(std::vector<int>           &error_and_warning_codes,
                                            std::vector<std::string>   &error_and_warning_messages,
                                            const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    if (pipelines.find(id) != pipelines.end()) {

      if (writeMeshON(cpi)) {
        writeMesh(cpi);
        return;
      }

      if (!canCoProcess()) {
        return;
      }

      error_and_warning_codes.clear();
      error_and_warning_messages.clear();

      if (!pipelines[id].canPerformCoProcessing()) {
        return;
      }

      vtkCPPythonPipeline  *pl              = pipelines[id].getPipeline();
      vtkCPDataDescription *dataDescription = pipelines[id].getDataDescription();
      coProcessor->AddPipeline(pl);
      coProcessor->CoProcess(dataDescription);

      vtkFieldData *fd = pipelines[id].getDataDescription()->GetUserData();
      vtkIntArray  *ec =
          vtkIntArray::SafeDownCast(fd->GetAbstractArray("catalyst_sierra_error_codes"));
      vtkStringArray *em =
          vtkStringArray::SafeDownCast(fd->GetAbstractArray("catalyst_sierra_error_messages"));

      if (ec && em && ec->GetNumberOfTuples() > 0 && em->GetNumberOfTuples() > 0 &&
          ec->GetNumberOfTuples() == em->GetNumberOfTuples()) {

        for (int i = 0; i < ec->GetNumberOfTuples(); i++) {
          error_and_warning_codes.push_back(ec->GetValue(i));
          error_and_warning_messages.push_back(em->GetValue(i));
        }
        fd->RemoveArray("catalyst_sierra_error_codes");
        fd->RemoveArray("catalyst_sierra_error_messages");
        vtkIntArray *ec = vtkIntArray::New();
        ec->SetName("catalyst_sierra_error_codes");
        vtkStringArray *em = vtkStringArray::New();
        em->SetName("catalyst_sierra_error_messages");
        fd->AddArray(ec);
        fd->AddArray(em);
        ec->Delete();
        em->Delete();
      }
      coProcessor->RemoveAllPipelines();
    }
  }

  void CatalystManager::SetTimeData(double currentTime, int timeStep,
                                    const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    if (pipelines.find(id) != pipelines.end()) {
      if (!pipelines[id].canSetTimeData()) {
        return;
      }
      pipelines[id].getDataDescription()->SetTimeData(currentTime, timeStep);
    }
  }

  void CatalystManager::logMemoryUsageAndTakeTimerReading(const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    if (this->logging.find(id) != this->logging.end()) {
      vtksys::SystemInformation sysInfo;
      double                    measurements[3];
      measurements[0]                = sysInfo.GetProcMemoryUsed() * (1.0 / 1024.0); // Store in MB
      measurements[1]                = sysInfo.GetHostMemoryUsed() * (1.0 / 1024.0);
      clock_t last_time              = this->logging[id].first.second;
      measurements[2]                = double(clock() - last_time) / (double)CLOCKS_PER_SEC;
      this->logging[id].first.second = clock();
      this->logging[id].second->InsertNextTuple(measurements);
    }
  }

  void CatalystManager::WriteToLogFile(const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id          = cpi.catalystPipelineID;
    std::string        logFileName = cpi.getLogFileName();

    if (this->logging.find(id) != this->logging.end()) {
      vtkProcessModule *pm         = vtkProcessModule::GetProcessModule();
      vtkMPIController *mpic       = vtkMPIController::SafeDownCast(pm->GetGlobalController());
      vtkDoubleArray   *logData    = this->logging[id].second;
      clock_t           begin_time = this->logging[id].first.first;
      if (mpic && mpic->GetNumberOfProcesses() > 1) {
        vtkDoubleArray *recvBufferMin = vtkDoubleArray::New();
        vtkDoubleArray *recvBufferMax = vtkDoubleArray::New();
        vtkDoubleArray *recvBufferSum = vtkDoubleArray::New();
        if (mpic->GetLocalProcessId() == 0) {
          recvBufferMin->SetNumberOfComponents(3);
          recvBufferMin->SetNumberOfTuples(logData->GetNumberOfTuples());

          recvBufferMax->SetNumberOfComponents(3);
          recvBufferMax->SetNumberOfTuples(logData->GetNumberOfTuples());

          recvBufferSum->SetNumberOfComponents(3);
          recvBufferSum->SetNumberOfTuples(logData->GetNumberOfTuples());
        }

        mpic->Reduce(logData, recvBufferMin, vtkCommunicator::MIN_OP, 0);
        mpic->Reduce(logData, recvBufferMax, vtkCommunicator::MAX_OP, 0);
        mpic->Reduce(logData, recvBufferSum, vtkCommunicator::SUM_OP, 0);

        if (mpic->GetLocalProcessId() == 0) {
          std::ofstream logfile;
          logfile.open(logFileName, ios::out | ios::app);
          for (int i = 0; i < logData->GetNumberOfTuples(); i++) {
            double min[3];
            double max[3];
            double sum[3];
            recvBufferMin->GetTuple(i, min);
            recvBufferMax->GetTuple(i, max);
            recvBufferSum->GetTuple(i, sum);
            logfile << double(clock() - begin_time) / (double)CLOCKS_PER_SEC << "," << min[0] << ","
                    << max[0] << "," << sum[0] / (double)mpic->GetNumberOfProcesses() << ","
                    << min[1] << "," << max[1] << ","
                    << sum[1] / (double)mpic->GetNumberOfProcesses() << "," << min[2] << ","
                    << max[2] << "," << sum[2] / (double)mpic->GetNumberOfProcesses() << "\n";
          }
          logfile.close();
        }
        recvBufferMin->Delete();
        recvBufferMax->Delete();
        recvBufferSum->Delete();
      }
      else {
        std::ofstream logfile;
        logfile.open(logFileName, ios::out | ios::app);
        for (int i = 0; i < logData->GetNumberOfTuples(); i++) {
          double data[3];
          logData->GetTuple(i, data);
          logfile << double(clock() - begin_time) / CLOCKS_PER_SEC << "," << data[0] << ","
                  << data[1] << "," << data[2] << "\n";
        }
        logfile.close();
      }
      logData->SetNumberOfTuples(0);
    }
  }

  bool CatalystManager::writeMeshON(const CatalystPipelineInfo &cpi)
  {

    bool retVal = false;

    CatalystPipelineID id = cpi.catalystPipelineID;

    if (pipelines.find(id) != pipelines.end()) {
      auto mw = pipelines[id].getMeshWriter();
      retVal  = mw->outputCatalystMeshOneFileON() || mw->outputCatalystMeshFilePerProcON();
    }

    return retVal;
  }

  void CatalystManager::writeMesh(const CatalystPipelineInfo &cpi)
  {

    CatalystPipelineID id = cpi.catalystPipelineID;

    if (pipelines.find(id) != pipelines.end()) {
      auto           mw   = pipelines[id].getMeshWriter();
      vtkDataObject *vobj = pipelines[id]
                                .getDataDescription()
                                ->GetInputDescriptionByName(cpi.catalystInputName.c_str())
                                ->GetGrid();
      int timeStep = pipelines[id].getDataDescription()->GetTimeStep();
      if (mw->outputCatalystMeshOneFileON()) {
        mw->writeCatalystMeshOneFile(vobj, timeStep);
      }
      if (mw->outputCatalystMeshFilePerProcON()) {
        mw->writeCatalystMeshFilePerProc(vobj, timeStep);
      }
    }
  }

  extern "C" {
  CatalystManagerBase *CreateCatalystManagerInstance()
  {
    CatalystManager *p = new CatalystManager();
    return (CatalystManagerBase *)p;
  }
  }

} // namespace Iovs
