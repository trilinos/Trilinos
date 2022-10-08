// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "include/ParaViewCatalystCGNSAdapter.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "include/vtkCGNSMultiBlockDataSet.h"

#include <vtkXMLPMultiBlockDataWriter.h>
#include <sstream>
#include <vtkStructuredGrid.h>
#include <vtkMultiProcessController.h>
#include <vtkXMLStructuredGridWriter.h>

class ParaViewCatalystCGNSAdapterImplementation
{
public:
  static ParaViewCatalystCGNSAdapterImplementation &getInstance()
  {
    static ParaViewCatalystCGNSAdapterImplementation instance;
    return instance;
  }

  void CleanupCatalyst()
  {
    if (this->dd) {
      this->dd->Delete();
      this->dd = 0;
    }

    if (this->pl) {
      this->pl->Delete();
      this->pl = 0;
    }

    if (this->coProcessor) {
      this->coProcessor->Finalize();
      this->coProcessor->Delete();
      this->coProcessor = 0;
    }
  }

  void CreateNewPipeline(const char *catalyst_python_filename,
                         const char *catalyst_sierra_block_json)
  {
    this->dd = vtkCPDataDescription::New();
    this->pl = vtkCPPythonScriptPipeline::New();

    vtkCGNSMultiBlockDataSet *cgnsmbds = vtkCGNSMultiBlockDataSet::New();
    dd->AddInput("input");
    dd->GetInputDescriptionByName("input")->SetGrid(cgnsmbds);
    cgnsmbds->Delete();

    if (this->pl->Initialize(catalyst_python_filename) == 0) {
      std::cerr << "Unable to initialize ParaView Catalyst with python script "
                << catalyst_python_filename << std::endl;
      std::cerr << "ParaView Catalyst CoProcessing will not be available." << std::endl;
      this->pl->Delete();
      this->pl = 0;
      this->dd->Delete();
      this->dd = 0;
      this->coProcessor->Delete();
      this->coProcessor = 0;
    }
  }

  void SetTimeData(double currentTime,
                   int timeStep)
  {
    if(this->dd) {
      this->dd->SetTimeData(currentTime, timeStep);

    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
    int myrank = controller->GetLocalProcessId();

      vtkCGNSMultiBlockDataSet *cgnsmbds =
       vtkCGNSMultiBlockDataSet::SafeDownCast(dd->GetInputDescriptionByName("input")->GetGrid());

      vtkMultiBlockDataSet *bases = vtkMultiBlockDataSet::SafeDownCast(cgnsmbds->GetBlock(0));

      if(!bases)
        std::cout << "bases is NULL\n";

      std::cout << "number of bases = " << bases->GetNumberOfBlocks() << "\n";

      vtkMultiBlockDataSet *base = vtkMultiBlockDataSet::SafeDownCast(bases->GetBlock(0));

      if(!base)
        std::cout << "base is NULL\n";

      std::cout << "number of base = " << base->GetNumberOfBlocks() << "\n";

      vtkMultiBlockDataSet *zones = vtkMultiBlockDataSet::SafeDownCast(base->GetBlock(0));

      if(!zones)
        std::cout << "zones is NULL\n";

      std::cout << "number of zones = " << zones->GetNumberOfBlocks() << "\n";

      vtkStructuredGrid *sg = vtkStructuredGrid::SafeDownCast(zones->GetBlock(0));

      if(!sg)
        std::cout << "sg is NULL\n";
        //it's okay for sg to be NULL. It means there was no data on
        //this process on this block

/*
      //pts is not used anywhere below, not sure why there was code getting it
      //(which caused core dump when sg was NULL)
      //vtkPoints* pts = sg->GetPoints();

      vtkXMLStructuredGridWriter* writer = vtkXMLStructuredGridWriter::New();
      writer->SetInputData(sg);
      std::ostringstream convert;
      convert << timeStep;
      std::ostringstream convert_rank;
      convert_rank << myrank;
      writer->SetFileName(std::string("test_" + convert.str()
                                        + "_" + convert_rank.str() +".vts").c_str());
      writer->Write();
      writer->Delete();
*/

     // vtkCGNSMultiBlockDataSet *cgnsmbds =
     //  vtkCGNSMultiBlockDataSet::SafeDownCast(dd->GetInputDescriptionByName("input")->GetGrid());
/*      vtkXMLPMultiBlockDataWriter* writer = vtkXMLPMultiBlockDataWriter::New();
      writer->SetInputData(cgnsmbds->GetBlock(0));
      std::ostringstream oss;
      oss << "." << writer->GetDefaultFileExtension();
      std::ostringstream convert;
      convert << timeStep;
      writer->SetFileName(std::string("test_" + convert.str() + oss.str()).c_str());
      if(myrank == 0) {
        writer->SetWriteMetaFile(1);
      }
      writer->Update();
      writer->Delete();
*/
    }
  }

  void PerformCoProcessing()
  {
    vtkCGNSMultiBlockDataSet *cgnsmbds =
    vtkCGNSMultiBlockDataSet::SafeDownCast(dd->GetInputDescriptionByName("input")->GetGrid());

    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
    int myrank = controller->GetLocalProcessId();
      vtkXMLPMultiBlockDataWriter* writer = vtkXMLPMultiBlockDataWriter::New();
      writer->SetInputData(cgnsmbds->GetBlock(0));
      std::ostringstream oss;
      oss << "." << writer->GetDefaultFileExtension();
      std::ostringstream convert;
      convert << 0;
      writer->SetFileName(std::string("test_" + convert.str() + oss.str()).c_str());
      if(myrank == 0) {
        writer->SetWriteMetaFile(1);
      }
      writer->Update();
      writer->Delete();

  }

  void CreateBase(int base_id,
                  const std::string& base_name)
  {
    if(this->dd) {
      vtkCGNSMultiBlockDataSet *cgnsmbds =
      vtkCGNSMultiBlockDataSet::SafeDownCast(dd->GetInputDescriptionByName("input")->GetGrid());
      cgnsmbds->CreateBase(base_id, base_name);
    }
  }

  void AddStructuredZoneData(int base_id,
                             int zone_id,
                             const std::string& zone_name,
                             const std::string& data_name,
                             int ni,
                             int nj,
                             int nk,
                             int comp_count,
                             bool is_cell_field,
                             char field_suffix_separator,
                             double* data,
                             int size)
  {
    if(this->dd) {
      vtkCGNSMultiBlockDataSet *cgnsmbds =
      vtkCGNSMultiBlockDataSet::SafeDownCast(dd->GetInputDescriptionByName("input")->GetGrid());
      cgnsmbds->AddStructuredZoneData(base_id,
                                      zone_id,
                                      zone_name,
                                      data_name,
                                      ni,
                                      nj,
                                      nk,
                                      comp_count,
                                      is_cell_field,
                                      field_suffix_separator,
                                      data,
                                      size);
    }
  }

  int parseFile(const std::string &filepath,
                CatalystParserInterface::parse_info &pinfo)
  {
  }

  int parseString(const std::string &s, CatalystParserInterface::parse_info &pinfo)
  {
  }

private:
  ParaViewCatalystCGNSAdapterImplementation()
  {
    this->coProcessor = vtkCPProcessor::New();
    this->coProcessor->Initialize();
  }

  ~ParaViewCatalystCGNSAdapterImplementation()
  {
  }

  ParaViewCatalystCGNSAdapterImplementation(ParaViewCatalystCGNSAdapterImplementation const &);
  void operator=(ParaViewCatalystCGNSAdapterImplementation const &);

  vtkCPProcessor *coProcessor;
  vtkCPPythonScriptPipeline *pl;
  vtkCPDataDescription *dd;
};

void ParaViewCatalystCGNSAdapter::CleanupCatalyst()
{
  ParaViewCatalystCGNSAdapterImplementation &pcsai =
      ParaViewCatalystCGNSAdapterImplementation::getInstance();

  pcsai.CleanupCatalyst();
}

void ParaViewCatalystCGNSAdapter::CreateNewPipeline(const char *catalyst_python_filename,
                                                    const char *catalyst_sierra_block_json)
{
  ParaViewCatalystCGNSAdapterImplementation &pcsai =
      ParaViewCatalystCGNSAdapterImplementation::getInstance();

  pcsai.CreateNewPipeline(catalyst_python_filename, catalyst_sierra_block_json);
}

void ParaViewCatalystCGNSAdapter::PerformCoProcessing()
{
  ParaViewCatalystCGNSAdapterImplementation &pcsai =
      ParaViewCatalystCGNSAdapterImplementation::getInstance();

  pcsai.PerformCoProcessing();
}

void ParaViewCatalystCGNSAdapter::CreateBase(int base_id,
                                             const std::string& base_name)
{
  ParaViewCatalystCGNSAdapterImplementation &pcsai =
      ParaViewCatalystCGNSAdapterImplementation::getInstance();

  pcsai.CreateBase(base_id,
                   base_name);
}

void ParaViewCatalystCGNSAdapter::AddStructuredZoneData(int base_id,
                                                        int zone_id,
                                                        const std::string& zone_name,
                                                        const std::string& data_name,
                                                        int ni,
                                                        int nj,
                                                        int nk,
                                                        int comp_count,
                                                        bool is_cell_field,
                                                        char field_suffix_separator,
                                                        double* data,
                                                        int size)
{
  ParaViewCatalystCGNSAdapterImplementation &pcsai =
      ParaViewCatalystCGNSAdapterImplementation::getInstance();

  pcsai.AddStructuredZoneData(base_id,
                              zone_id,
                              zone_name,
                              data_name,
                              ni,
                              nj,
                              nk,
                              comp_count,
                              is_cell_field,
                              field_suffix_separator,
                              data,
                              size);
}

void ParaViewCatalystCGNSAdapter::SetTimeData(double currentTime, int timeStep)
{
  ParaViewCatalystCGNSAdapterImplementation &pcsai =
      ParaViewCatalystCGNSAdapterImplementation::getInstance();

  pcsai.SetTimeData(currentTime, timeStep);
}

int ParaViewCatalystCGNSAdapter::parseFile(const std::string &filepath,
                                           CatalystParserInterface::parse_info &pinfo)
{
  return CatalystParserInterface::parseFile(filepath, pinfo);
}

int ParaViewCatalystCGNSAdapter::parseString(const std::string &s,
                                             CatalystParserInterface::parse_info &pinfo)
{
  return CatalystParserInterface::parseString(s, pinfo);
}

extern "C" {
ParaViewCatalystCGNSAdapterBase *ParaViewCatalystCGNSAdapterCreateInstance()
{
  ParaViewCatalystCGNSAdapterBase *t(new ParaViewCatalystCGNSAdapter());
  return t;
}
}
