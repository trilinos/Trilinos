#include "vtkCppPipe.h"
#include "CppPipeProxies.h"
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkObjectFactory.h>
#include <vtkSMSessionProxyManager.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>
#include <vtksys/SystemInformation.hxx>

#include <stdio.h>
#include <stdlib.h>

vtkStandardNewMacro(vtkCppPipe);

//----------------------------------------------------------------------------
vtkCppPipe::vtkCppPipe()
{
}

//----------------------------------------------------------------------------
vtkCppPipe::~vtkCppPipe()
{
}

//----------------------------------------------------------------------------
int vtkCppPipe::Initialize()
{
  return 1;
}

//----------------------------------------------------------------------------
int vtkCppPipe::RequestDataDescription(
  vtkCPDataDescription* dataDescription)
{
  if(!dataDescription)
    {
    vtkWarningMacro("dataDescription is NULL.");
    return 0;
    }

  if(dataDescription->GetForceOutput() == true)
    {
    dataDescription->GetInputDescriptionByName("input")->AllFieldsOn();
    dataDescription->GetInputDescriptionByName("input")->GenerateMeshOn();
    }
  return 1;
}

//----------------------------------------------------------------------------
int vtkCppPipe::CoProcess(
  vtkCPDataDescription* dataDescription)
{
  if(!dataDescription)
    {
    vtkWarningMacro("DataDescription is NULL");
    return 0;
    }

  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
    dataDescription->GetInputDescriptionByName("input")->GetGrid());

  if(!grid)
    {
    vtkWarningMacro("DataDescription is missing input unstructured grid.");
    return 0;
    }
  if(this->RequestDataDescription(dataDescription) == 0)
    {
    return 1;
    }

  vtkFieldData* fd = dataDescription->GetUserData();
  std::string image_basename = "";

  Json::Value json;
  if(fd)
    {
    vtkStringArray* sa = vtkStringArray::SafeDownCast(
        fd->GetAbstractArray("catalyst_sierra_data"));

    if(sa)
      {
      Json::Reader reader;
      bool rc = reader.parse(sa->GetValue(0), json, false);
      if(!rc)
        {
        vtkWarningMacro("Parse of catalyst JSON failed.");
        return 0;
        }
      image_basename = sa->GetValue(2);
      }
    else 
      {
      vtkWarningMacro("DataDescription is missing catalyst_sierra_data string array.");
      return 0;
      }
    }
  else
    {
    vtkWarningMacro("DataDescription is missing user data.");
    return 0;
    }
  
  CppPipeProxies proxies(dataDescription);

  for(Json::ValueIterator itr = json["imageset blocks"].begin() ;
                          itr != json["imageset blocks"].end() ; itr++)
    {
    
    proxies.RenderImageSet(itr.key().asString(),
                           image_basename,
                           json);
    }
  
  return 1;
}

//----------------------------------------------------------------------------
void vtkCppPipe::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
