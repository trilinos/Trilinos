
#include "include/ParaViewCatalystSierraAdaptor.h"
#include "vtkExodusIIMultiBlockDataSet.h"
#include "vtkCPProcessor.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"

#ifdef USE_CPP_PIPE
#include "vtkCppPipe.h"
#else
#include "vtkCPPythonScriptPipeline.h"
#endif

#include "vtkVariant.h"
#include "vtkStringArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkFieldData.h"
#include "vtkProcessModule.h"
#include "vtkMPIController.h"
#include <sstream>
#include <vtksys/SystemInformation.hxx>
#include <time.h>
#include <fstream>

#ifdef USE_CPP_PIPE
typedef std::pair <vtkCppPipe*, vtkCPDataDescription*> PipelineDataDescPair;
#else
typedef std::pair <vtkCPPythonScriptPipeline*, vtkCPDataDescription*> PipelineDataDescPair;
#endif

typedef std::pair <clock_t, clock_t> TimerPair;
typedef std::pair <TimerPair, vtkDoubleArray*> LoggingPair;

class ParaViewCatalystSierraAdaptorImplementation
{
public:

  static ParaViewCatalystSierraAdaptorImplementation& getInstance()
  {
    static ParaViewCatalystSierraAdaptorImplementation instance;
    return instance;
  }

  void CleanupCatalyst()
  {
    if(this->coProcessor)
      {
      this->coProcessor->Finalize();
      this->coProcessor->Delete();
      this->coProcessor = 0;
      }
  }

  void DeletePipeline(const char* results_output_filename)
  {
    if (this->pipelines.find(results_output_filename) !=
        this->pipelines.end() )
      {
      this->pipelines[results_output_filename].first->Delete();
      this->pipelines[results_output_filename].second->Delete();
      this->pipelines.erase(results_output_filename);
      }

    if(this->logging.find(results_output_filename) !=
       this->logging.end())
      {
      this->logging[results_output_filename].second->Delete();
      this->logging.erase(results_output_filename);
      }
  }

  void CreateNewPipeline(const char* catalyst_python_filename,
                         const char* catalyst_sierra_block_json,
                         const char* catalyst_sierra_separator_character,
                         const char* catalyst_sierra_input_deck_name,
                         int UnderscoreVectors,
                         int ApplyDisplacements,
                         const char* restart_tag,
                         int enable_logging,
                         int debug_level,
                         const char* results_output_filename,
                         const char* catalyst_output_directory,
                         std::vector<std::string>& catalyst_sierra_data)
  {
    //std::cerr << "CreateNewPipeline enetered\n";
    if (enable_logging)
      {
#if 1
      TimerPair tp = std::make_pair(clock(),clock());
      vtkDoubleArray* da = vtkDoubleArray::New();
      da->SetNumberOfComponents(3);
      LoggingPair lp = std::make_pair(tp,da);
      this->logging[results_output_filename] = lp;

      vtkProcessModule* pm = vtkProcessModule::GetProcessModule();
      vtkMPIController* mpic = vtkMPIController::SafeDownCast(pm->GetGlobalController());
      std::string s(results_output_filename);
      if(mpic && mpic->GetNumberOfProcesses() > 1)
        {
        if(mpic->GetLocalProcessId() == 0)
          {
          ofstream logfile;
          logfile.open((s + ".catalyst.log").c_str(), ios::out | ios::trunc );
          logfile << "# ELAPSED TIME (S)"
                  << ",PROC MEM USED - MIN (KiB)"
                  << ",PROC MEM USED - MAX (KiB)"
                  << ",PROC MEM USED - AVG (KiB)"
                  << ",HOST MEM USED - MIN (KiB)"
                  << ",HOST MEM USED - MAX (KiB)"
                  << ",HOST MEM USED - AVG (KiB)"
                  << ",TIME SINCE LAST LOG - MIN (S)"
                  << ",TIME SINCE LAST LOG - MAX (S)"
                  << ",TIME SINCE LAST LOG - AVG (S)" << "\n";
          logfile.close();
          }
        }
      else
        {
        ofstream logfile;
        logfile.open((s + ".catalyst.log").c_str(), ios::out | ios::trunc );
        logfile << "# ELAPSED TIME (S)"
                << ",PROC MEM USED (KiB)"
                << ",HOST MEM USED (KiB)"
                << ",TIME SINCE LAST LOG (S)" << "\n";
        logfile.close();
        }
#endif
      }

    if ( this->pipelines.find(results_output_filename) == this->pipelines.end() )
      {
      vtkCPDataDescription* dd = vtkCPDataDescription::New();

#ifdef USE_CPP_PIPE
      vtkCppPipe *pl = vtkCppPipe::New();
#else
      vtkCPPythonScriptPipeline *pl = vtkCPPythonScriptPipeline::New();
#endif

      vtkExodusIIMultiBlockDataSet* mbds = vtkExodusIIMultiBlockDataSet::New();
      mbds->SetUnderscoreVectors(UnderscoreVectors);
      mbds->SetApplyDisplacements(ApplyDisplacements);
      dd->AddInput("input");
      dd->GetInputDescriptionByName("input")->SetGrid(mbds);
      mbds->Delete();
      PipelineDataDescPair pddp = std::make_pair(pl,dd);
      this->pipelines[results_output_filename] = pddp;
      }

    if(catalyst_sierra_block_json)
      {
      vtkFieldData* fd = vtkFieldData::New();
      vtkStringArray* sa = vtkStringArray::New();
      sa->SetName("catalyst_sierra_data");
      vtkIntArray* ec = vtkIntArray::New();
      ec->SetName("catalyst_sierra_error_codes");
      vtkStringArray* em = vtkStringArray::New();
      em->SetName("catalyst_sierra_error_messages");
      sa->InsertNextValue(catalyst_sierra_block_json);
      sa->InsertNextValue(catalyst_sierra_separator_character);
      sa->InsertNextValue(catalyst_sierra_input_deck_name);
      sa->InsertNextValue(restart_tag);
      if(enable_logging)
        {
        sa->InsertNextValue("True");
        }
      else
        {
        sa->InsertNextValue("");
        }
      std::stringstream ss;
      ss << debug_level;
      sa->InsertNextValue(ss.str().c_str());
      ss.clear();
      sa->InsertNextValue(results_output_filename);
      sa->InsertNextValue(catalyst_output_directory);

      for(int i=0;i<catalyst_sierra_data.size();i++)
        sa->InsertNextValue(catalyst_sierra_data[i]);

      fd->AddArray(sa);
      fd->AddArray(ec);
      fd->AddArray(em);
      this->pipelines[results_output_filename].second->SetUserData(fd);
      fd->Delete();
      sa->Delete();
      ec->Delete();
      em->Delete();
      }

#ifdef USE_CPP_PIPE
    if(this->pipelines[results_output_filename].first->Initialize() == 0)
      {
      std::cerr << "Unable to initialize ParaView Catalyst CPP Pipeline.\n"
                << "ParaView Catalyst CoProcessing will not be available." << std::endl;
      this->coProcessor->Delete();
      this->coProcessor = 0;
      }
#else
    if(this->pipelines[results_output_filename].first->Initialize(catalyst_python_filename) == 0)
      {
      std::cerr << "Unable to initialize ParaView Catalyst with python script "
                << catalyst_python_filename << std::endl;
      std::cerr << "ParaView Catalyst CoProcessing will not be available." << std::endl;
      this->coProcessor->Delete();
      this->coProcessor = 0;
      }
#endif

    //std::cerr << "CreateNewPipeline returning\n";
  }

  void SetTimeData(double currentTime,
                   int timeStep,
                   const char* results_output_filename)
  {
    if(!this->coProcessor)
      return;

    if(this->pipelines.find(results_output_filename) !=
       this->pipelines.end())
      {
      this->pipelines[results_output_filename].second->SetTimeData(currentTime, timeStep);
      }
  }

  void CoProcess(const char* results_output_filename,
                 std::vector<int>& error_and_warning_codes,
                 std::vector<std::string>& error_and_warning_messages)
  {
    //std::cerr << "CoProcess entered\n";
    if(!this->coProcessor)
      return;

    if(this->pipelines.find(results_output_filename) !=
       this->pipelines.end())
      {
      error_and_warning_codes.clear();
      error_and_warning_messages.clear();

#ifdef USE_CPP_PIPE
      vtkCppPipe* pl = this->pipelines[results_output_filename].first;
#else
      vtkCPPythonScriptPipeline* pl = this->pipelines[results_output_filename].first;
#endif

      vtkCPDataDescription* dd = this->pipelines[results_output_filename].second;
      vtkExodusIIMultiBlockDataSet* mbds = vtkExodusIIMultiBlockDataSet::SafeDownCast(
                                           dd->GetInputDescriptionByName("input")->GetGrid());
      mbds->Register(0);
      pl->Register(0);
      vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::New();
      mb->ShallowCopy(mbds);
      dd->GetInputDescriptionByName("input")->SetGrid(mb);
      this->coProcessor->AddPipeline(pl);

      this->coProcessor->CoProcess(dd);

      vtkFieldData* fd = this->pipelines[results_output_filename].second->GetUserData();
      vtkIntArray* ec = vtkIntArray::SafeDownCast(fd->GetAbstractArray("catalyst_sierra_error_codes"));
      vtkStringArray* em = vtkStringArray::SafeDownCast(fd->GetAbstractArray("catalyst_sierra_error_messages"));

      if(ec && em &&
         ec->GetNumberOfTuples() > 0 &&
         em->GetNumberOfTuples() > 0 &&
         ec->GetNumberOfTuples() == em->GetNumberOfTuples())
        {
        for(int i = 0; i < ec->GetNumberOfTuples(); i++)
          {
          error_and_warning_codes.push_back(ec->GetValue(i));
          error_and_warning_messages.push_back(em->GetValue(i));
          }
        fd->RemoveArray("catalyst_sierra_error_codes");
        fd->RemoveArray("catalyst_sierra_error_messages");
        vtkIntArray* ec = vtkIntArray::New();
        ec->SetName("catalyst_sierra_error_codes");
        vtkStringArray* em = vtkStringArray::New();
        em->SetName("catalyst_sierra_error_messages");
        fd->AddArray(ec);
        fd->AddArray(em);
        ec->Delete();
        em->Delete();
        }

      this->coProcessor->RemoveAllPipelines();
      dd->GetInputDescriptionByName("input")->SetGrid(mbds);
      mbds->Delete();
      pl->Delete();
      mb->Delete();
      }
    //std::cerr << "CoProcess returning\n";
  }

  vtkExodusIIMultiBlockDataSet* GetMultiBlockDataSet(const char* results_output_filename)
  {
    if(this->pipelines.find(results_output_filename) !=
       this->pipelines.end())
      {
      vtkCPDataDescription* dd = this->pipelines[results_output_filename].second;
      return( vtkExodusIIMultiBlockDataSet::SafeDownCast(
              dd->GetInputDescriptionByName("input")->GetGrid()) );
      }
    else
      return 0;
  }

  void logMemoryUsageAndTakeTimerReading(const char* results_output_filename)
  {
    if(this->logging.find(results_output_filename) !=
       this->logging.end())
      {
#if 1
      vtksys::SystemInformation sysInfo;
      vtkProcessModule* pm = vtkProcessModule::GetProcessModule();
      vtkMPIController* mpic = vtkMPIController::SafeDownCast(pm->GetGlobalController());
      double measurements[3];
      measurements[0] = sysInfo.GetProcMemoryUsed()*(1.0/1024.0);  // Store in MB
      measurements[1] = sysInfo.GetHostMemoryUsed()*(1.0/1024.0);
      clock_t last_time = this->logging[results_output_filename].first.second;
      measurements[2] = double( clock () - last_time ) /  (double) CLOCKS_PER_SEC;
      this->logging[results_output_filename].first.second= clock();
      this->logging[results_output_filename].second->InsertNextTuple(measurements);
#endif
      }
  }

  void WriteToLogFile(const char* results_output_filename)
  {
    if(this->logging.find(results_output_filename) !=
       this->logging.end())
      {
#if 1
      vtkProcessModule* pm = vtkProcessModule::GetProcessModule();
      vtkMPIController* mpic = vtkMPIController::SafeDownCast(pm->GetGlobalController());
      vtkDoubleArray* logData = this->logging[results_output_filename].second;
      std::string s(results_output_filename);
      clock_t begin_time = this->logging[results_output_filename].first.first;
      if(mpic && mpic->GetNumberOfProcesses() > 1)
        {
        vtkDoubleArray* recvBufferMin = vtkDoubleArray::New();
        vtkDoubleArray* recvBufferMax = vtkDoubleArray::New();
        vtkDoubleArray* recvBufferSum = vtkDoubleArray::New();
        if(mpic->GetLocalProcessId() == 0)
          {
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

        if(mpic->GetLocalProcessId() == 0)
          {
          ofstream logfile;
          logfile.open((s + ".catalyst.log").c_str(), ios::out | ios::app );
          for(int i = 0;i<logData->GetNumberOfTuples();i++)
            {
            double min[3];
            double max[3];
            double sum[3];
            recvBufferMin->GetTuple(i,min);
            recvBufferMax->GetTuple(i,max);
            recvBufferSum->GetTuple(i,sum);
            logfile << double( clock () - begin_time ) / (double) CLOCKS_PER_SEC
                    << "," << min[0] << "," << max[0] << ","
                    << sum[0]/(double) mpic->GetNumberOfProcesses()
                    << "," << min[1] << "," << max[1] << ","
                    << sum[1]/(double) mpic->GetNumberOfProcesses()
                    << "," << min[2] << "," << max[2] << ","
                    << sum[2]/(double) mpic->GetNumberOfProcesses() << "\n";
            }
          logfile.close();
          }
        recvBufferMin->Delete();
        recvBufferMax->Delete();
        recvBufferSum->Delete();
        }
      else
        {
        ofstream logfile;
        logfile.open((s + ".catalyst.log").c_str(), ios::out | ios::app );
        for(int i = 0;i<logData->GetNumberOfTuples();i++)
          {
          double data[3];
          logData->GetTuple(i,data);
          logfile << double( clock () - begin_time ) / CLOCKS_PER_SEC
                  << "," << data[0] << "," << data[1] << "," << data[2] << "\n";
          }
        logfile.close();
        }
      logData->SetNumberOfTuples(0);
#endif
      }
  }

private:

  ParaViewCatalystSierraAdaptorImplementation()
  {
  this->coProcessor = vtkCPProcessor::New();
  this->coProcessor->Initialize();
  }

  ~ParaViewCatalystSierraAdaptorImplementation()
  {
  for(std::map<std::string, PipelineDataDescPair>::iterator it = this->pipelines.begin();
      it != this->pipelines.end(); it++)
    {
    this->DeletePipeline(it->first.c_str());
    }
  }

  ParaViewCatalystSierraAdaptorImplementation(ParaViewCatalystSierraAdaptorImplementation const&);
  void operator=(ParaViewCatalystSierraAdaptorImplementation const&);

  vtkCPProcessor *coProcessor; 
  std::map<std::string, PipelineDataDescPair> pipelines;
  std::map<std::string, LoggingPair> logging;
};

#if 0
ParaViewCatalystSierraAdaptor::ParaViewCatalystSierraAdaptor()
{
}

ParaViewCatalystSierraAdaptor::~ParaViewCatalystSierraAdaptor()
{
}
#endif

void ParaViewCatalystSierraAdaptor::DeletePipeline(const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();

  pcsai.DeletePipeline(results_output_filename);
}

void ParaViewCatalystSierraAdaptor::CleanupCatalyst()
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();

  pcsai.CleanupCatalyst();
}

void ParaViewCatalystSierraAdaptor::CreateNewPipeline(const char* catalyst_python_filename,
                                                      const char* catalyst_sierra_block_json,
                                                      const char* catalyst_sierra_separator_character,
                                                      const char* catalyst_sierra_input_deck_name,
                                                      int UnderscoreVectors,
                                                      int ApplyDisplacements,
                                                      const char* restart_tag,
                                                      int enable_logging,
                                                      int debug_level,
                                                      const char* results_output_filename,
                                                      const char* catalyst_output_directory,
                                                      std::vector<std::string>& catalyst_sierra_data)
{
  //std::cerr << "CreateNewPipeline enetered 2\n";
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();

  pcsai.CreateNewPipeline(catalyst_python_filename,
                          catalyst_sierra_block_json,
                          catalyst_sierra_separator_character,
                          catalyst_sierra_input_deck_name,
                          UnderscoreVectors,
                          ApplyDisplacements,
                          restart_tag,
                          enable_logging,
                          debug_level,
                          results_output_filename,
                          catalyst_output_directory,
                          catalyst_sierra_data);
  //std::cerr << "CreateNewPipeline returning 2\n";
}

void ParaViewCatalystSierraAdaptor::PerformCoProcessing(const char* results_output_filename,
                                                        std::vector<int>& error_and_warning_codes,
                                                        std::vector<std::string>& error_and_warning_messages)
{
  //std::cerr << "ParaViewCatalystSierraAdaptor::PerformCoProcessing entered\n";
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();

  pcsai.CoProcess(results_output_filename,
                  error_and_warning_codes,
                  error_and_warning_messages);
  //std::cerr << "ParaViewCatalystSierraAdaptor::PerformCoProcessing returning\n";
}

void ParaViewCatalystSierraAdaptor::SetTimeData(double currentTime,
                                                int timeStep,
                                                const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();

  pcsai.SetTimeData(currentTime,
                    timeStep,
                    results_output_filename);
}

void ParaViewCatalystSierraAdaptor::CreateGlobalVariable(std::vector<std::string>& component_names,
                                                         const int* data,
                                                         const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int) 0);
    mbds->CreateGlobalVariable(component_names,
                               v,
                               data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateGlobalVariable(std::vector<std::string>& component_names,
                                                         const double* data,
                                                         const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((double) 0.0);
    mbds->CreateGlobalVariable(component_names,
                               v,
                               data);
    }
}

void ParaViewCatalystSierraAdaptor::InitializeGlobalPoints(int num_points,
                                                           int dimension,
                                                           const double* data,
                                                           const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    mbds->InitializeGlobalPoints(num_points,
                                 dimension,
                                 data);
    }
}

void ParaViewCatalystSierraAdaptor::InitializeElementBlocks(const std::vector<int>& element_block_id_list,
                                                            const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    mbds->InitializeElementBlocks(element_block_id_list);
    }
}

void ParaViewCatalystSierraAdaptor::ReleaseMemory(const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    mbds->ReleaseMemory();
    }
  pcsai.WriteToLogFile(results_output_filename);
}

void ParaViewCatalystSierraAdaptor::CreateElementBlock(const char* elem_block_name,
                                                       int elem_block_id,
                                                       const std::string& elem_type,
                                                       int nodes_per_elem,
                                                       int num_elem,
                                                       const int64_t* global_elem_ids,
                                                       int* connectivity,
                                                       const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int) 0);
    mbds->CreateElementBlock(elem_block_name,
                             elem_block_id,
                             elem_type,
                             nodes_per_elem,
                             num_elem,
                             v,
                             global_elem_ids,
                             connectivity);
    }
}

void ParaViewCatalystSierraAdaptor::CreateElementBlock(const char* elem_block_name,
                                                       int elem_block_id,
                                                       const std::string& elem_type,
                                                       int nodes_per_elem,
                                                       int num_elem,
                                                       const int64_t* global_elem_ids,
                                                       int64_t* connectivity,
                                                       const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int64_t) 0);
    mbds->CreateElementBlock(elem_block_name,
                             elem_block_id,
                             elem_type,
                             nodes_per_elem,
                             num_elem,
                             v,
                             global_elem_ids,
                             connectivity);
    }
}


void ParaViewCatalystSierraAdaptor::CreateNodeSet(const char* node_set_name,
                                                  int node_set_id,
                                                  int num_ids,
                                                  const int* data,
                                                  const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int) 0);
    mbds->CreateNodeSet(node_set_name,
                        node_set_id,
                        num_ids,
                        v,
                        data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateNodeSet(const char* node_set_name,
                                                  int node_set_id,
                                                  int num_ids,
                                                  const int64_t* data,
                                                  const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int64_t) 0);
    mbds->CreateNodeSet(node_set_name,
                        node_set_id,
                        num_ids,
                        v,
                        data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateSideSet(/*const char* side_set_name,*/
                                                  const char* ss_owner_name,
                                                  int side_set_id,
                                                  int num_ids,
                                                  const int* element_ids,
                                                  const int* face_ids,
                                                  const char* results_output_filename)
{
    /*NOTE: Jeff Mauldin JAM 2015Oct8
     CreateSideSet is called once for each block which the sideset
     spans, and the side_set_name for the side set is the ss_owner_name
     with additional characters to indicate which block we are doing.
     The current implementation of the sierra sideset construction
     creates a single independent sideset and collects all the
     nodes and elements from the side set from each block spanned by
     the sideset into that single sideset.  It needs to have the
     ss_owner_name, not the side_set_name, because that is the name
     in the input deck for the sideset for reference for things like
     extractblock.  It may become necessary at a later date to
     pass in both, but for now we
     are just passing in ss_owner_name to give us correct
     functionality while not chaning the function interface*/
    /*
    std::cerr << "ParaViewCatalystSierraAdaptor::CreateSideSet entered (1)\n"
      //"side_set_name: " << side_set_name << "\n"
      "ss_owner_name: " << ss_owner_name << "\n"
      "side_set_id: " << side_set_id << "\n"
      "num_ids: " << num_ids << "\n"
      "element_ids: " << element_ids << "\n"
      "face_ids: " << face_ids << "\n"
      "results_output_filename: " << results_output_filename << "\n";
    */

  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int) 0);
    mbds->CreateSideSet(/*side_set_name,*/
                        ss_owner_name,
                        side_set_id,
                        num_ids,
                        v,
                        element_ids,
                        face_ids);
    }
  //std::cerr << "ParaViewCatalystSierraAdaptor::CreateSideSet returning (1)\n";
}

void ParaViewCatalystSierraAdaptor::CreateSideSet(/*const char* side_set_name,*/
                                                  const char* ss_owner_name,
                                                  int side_set_id,
                                                  int num_ids,
                                                  const int64_t* element_ids,
                                                  const int64_t* face_ids,
                                                  const char* results_output_filename)
{
    /*NOTE: Jeff Mauldin JAM 2015Oct8
     CreateSideSet is called once for each block which the sideset
     spans, and the side_set_name for the side set is the ss_owner_name
     with additional characters to indicate which block we are doing.
     The current implementation of the sierra sideset construction
     creates a single independent sideset and collects all the
     nodes and elements from the side set from each block spanned by
     the sideset into that single sideset.  It needs to have the
     ss_owner_name, not the side_set_name, because that is the name
     in the input deck for the sideset for reference for things like
     extractblock.  It may become necessary at a later date to
     pass in both, but for now we
     are just passing in ss_owner_name to give us correct
     functionality while not chaning the function interface*/
    /*
    std::cerr << "ParaViewCatalystSierraAdaptor::CreateSideSet entered (2)\n"
      //"side_set_name: " << side_set_name << "\n"
      "ss_owner_name: " << ss_owner_name << "\n"
      "side_set_id: " << side_set_id << "\n"
      "num_ids: " << num_ids << "\n"
      "element_ids: " << element_ids << "\n"
      "face_ids: " << face_ids << "\n"
      "results_output_filename: " << results_output_filename << "\n";
      */

  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int64_t) 0);
    mbds->CreateSideSet(/*side_set_name,*/
                        ss_owner_name,
                        side_set_id,
                        num_ids,
                        v,
                        element_ids,
                        face_ids);
    }
  //std::cerr << "ParaViewCatalystSierraAdaptor::CreateSideSet returning (2)\n";
}

void ParaViewCatalystSierraAdaptor::CreateElementVariable(std::vector<std::string>& component_names,
                                                          int elem_block_id,
                                                          const double* data,
                                                          const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((double) 0.0);
    mbds->CreateElementVariable(component_names,
                                elem_block_id,
                                v,
                                data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateElementVariable(std::vector<std::string>& component_names,
                                                          int elem_block_id,
                                                          const int* data,
                                                          const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int) 0);
    mbds->CreateElementVariable(component_names,
                                elem_block_id,
                                v,
                                data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateElementVariable(std::vector<std::string>& component_names,
                                                          int elem_block_id,
                                                          const int64_t* data,
                                                          const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int64_t) 0);
    mbds->CreateElementVariable(component_names,
                                elem_block_id,
                                v,
                                data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateNodalVariable(std::vector<std::string>& component_names,
                                                        const double* data,
                                                        const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((double) 0.0);
    mbds->CreateNodalVariable(component_names,
                              v,
                              data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateNodalVariable(std::vector<std::string>& component_names,
                                                        const int* data,
                                                        const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int) 0);
    mbds->CreateNodalVariable(component_names,
                              v,
                              data);
    }
}

void ParaViewCatalystSierraAdaptor::CreateNodalVariable(std::vector<std::string>& component_names,
                                                        const int64_t* data,
                                                        const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();
  vtkExodusIIMultiBlockDataSet* mbds = pcsai.GetMultiBlockDataSet(results_output_filename);
  if(mbds)
    {
    vtkVariant v((int64_t) 0);
    mbds->CreateNodalVariable(component_names,
                              v,
                              data);
    }
}

void ParaViewCatalystSierraAdaptor::logMemoryUsageAndTakeTimerReading(const char* results_output_filename)
{
  ParaViewCatalystSierraAdaptorImplementation& pcsai =
  ParaViewCatalystSierraAdaptorImplementation::getInstance();

  pcsai.logMemoryUsageAndTakeTimerReading(results_output_filename);
}

extern "C" {
ParaViewCatalystSierraAdaptorBase *ParaViewCatalystSierraAdaptorCreateInstance() {
  //std::cerr << "ParaViewCatalystSierraAdaptorCreateInstance entered\n";
  ParaViewCatalystSierraAdaptorBase * t(new ParaViewCatalystSierraAdaptor());
  //std::cerr << "ParaViewCatalystSierraAdaptorCreateInstance returning\n";
  return t;
}
}

#define USE_STK_DIAG_USER_PLUGIN 0
#if USE_STK_DIAG_USER_PLUGIN
extern "C"
void
dl_register()
{
  ParaViewCatalystSierraAdaptorBaseFactory::Register<ParaViewCatalystSierraAdaptor>("ParaViewCatalystSierraAdaptor",
                                                                                     ParaViewCatalystSierraAdaptorCreateInstance);
}
#endif


