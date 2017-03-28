
#include "CppPipeProxies.h"
#include "CppPipeSlice.h"
#include "CppPipeClip.h"
#include "CppPipeBoxClip.h"
#include "CppPipeThreshold.h"
#include "CppPipeContour.h"

#include <vtkDoubleArray.h>
#include <vtkPVDataInformation.h>
#include <vtkCPDataDescription.h>
#include <vtkSMSourceProxy.h>
#include <vtkSMSessionProxyManager.h>
#include <vtkSMProxyManager.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPVTrivialProducer.h>
#include <vtkCPInputDataDescription.h>
#include <vtkSMRenderViewProxy.h>
#include <vtkMultiProcessController.h>
#include <vtkSMPropertyHelper.h>
#include <vtkSMRepresentationProxy.h>
#include <vtksys/SystemTools.hxx>
#include <sstream>

CppPipeProxies::CppPipeProxies(vtkCPDataDescription* dataDescription)
{
  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
    dataDescription->GetInputDescriptionByName("input")->GetGrid());

  vtkSMProxyManager* proxyManager = vtkSMProxyManager::GetProxyManager();
  vtkSMSessionProxyManager* sessionProxyManager =
    proxyManager->GetActiveSessionProxyManager();

  this->producer = vtkSMSourceProxy::SafeDownCast(
    sessionProxyManager->NewProxy("sources", "PVTrivialProducer"));
  this->producer->UpdateVTKObjects();
  vtkObjectBase* clientSideObject = this->producer->GetClientSideObject();
  vtkPVTrivialProducer* realProducer =
    vtkPVTrivialProducer::SafeDownCast(clientSideObject);
  realProducer->SetOutput(grid);

  this->render = vtkSMRenderViewProxy::SafeDownCast(
    sessionProxyManager->NewProxy("views", "RenderView"));

  this->dataDescription = dataDescription;
  this->globalBounds = vtkDoubleArray::New();
  this->GetGlobalDataBoundsParallel(this->globalBounds);

  this->slice = 0;
  this->clip_plane = 0;
  this->slice_plane = 0;
  this->clip = 0;
  this->threshold = 0;
  this->contour = 0;
  this->box = 0;
}

CppPipeProxies::~CppPipeProxies()
{
  if(this->producer)
    this->producer->Delete();

  if(this->render)
    this->render->Delete();
 
  if(this->globalBounds)
    this->globalBounds->Delete();

  if(this->slice)
    this->slice->Delete();

  if(this->clip)
    this->clip->Delete();

  if(this->threshold)
    this->threshold->Delete();

  if(this->contour)
    this->contour->Delete();

  if(this->clip_plane)
    this->clip_plane->Delete();

  if(this->slice_plane)
    this->slice_plane->Delete();

  if(this->box)
    this->box->Delete();
}

void CppPipeProxies::RenderImageSet(const std::string& ims_name,
                                     const std::string& image_basename,
                                     const Json::Value& json)
{
  CppPipeImageSet ims(json["imageset blocks"][ims_name]);
  std::string bd = ims.GetBaseDirectoryName();
  if(vtkMultiProcessController::GetGlobalController()->GetLocalProcessId() == 0)
    vtksys::SystemTools::MakeDirectory(bd);

  vtkSMRepresentationProxy* representation = this->SetupOperations(ims_name,
                                                                   json);

  if(ims_name == "\\#ImageSet\\$ Shortcut")
    {
    CppPipeMultiCamera8 mc(json["imageset blocks"][ims_name],
                            this->globalBounds);
    CppPipeRepresentation phrep(json["imageset blocks"][ims_name],
                                 this->dataDescription,
                                 this->render);
    this->RenderMultiCamera8(ims, mc, phrep, image_basename, representation);
    }
  else
    {
    CppPipeRepresentation phrep(this->dataDescription,
                                 this->render);
    if(json["imageset blocks"][ims_name].isMember("representation"))
      {
      std::string representation_name = json["imageset blocks"][ims_name]["representation"].asString();
      phrep.SetJSON(json["representation blocks"][representation_name]);
      }

    if(json["imageset blocks"][ims_name].isMember("camera"))
      {
      std::string camera_name = json["imageset blocks"][ims_name]["camera"].asString();
      if(json["camera blocks"][camera_name]["camera type"].asString() == "camera")
        {
        CppPipeCamera c(json["camera blocks"][camera_name],
                         this->globalBounds);
        this->RenderCamera(ims, c, phrep, image_basename, representation);
        }
      else
        {
        CppPipeMultiCamera8 mc(json["camera blocks"][camera_name],
                                this->globalBounds);
        this->RenderMultiCamera8(ims, mc, phrep, image_basename, representation);
        }
      }
    else
      {
      CppPipeMultiCamera8 mc(this->globalBounds);
      this->RenderMultiCamera8(ims, mc, phrep, image_basename, representation);
      }
    }

  representation->Delete();
}

void CppPipeProxies::RenderMultiCamera8(CppPipeImageSet& ims,
                                         CppPipeMultiCamera8& mc,
                                         CppPipeRepresentation& phrep, 
                                         const std::string& image_basename,
                                         vtkSMRepresentationProxy* representation)
{
  ims.InitializeProxy(this->render);
  ims.InitializeProxy(representation);
  phrep.InitializeProxy(representation);
 
  std::string base = ims.GetImageBaseName();
  if(base.empty())
    base = image_basename;

  std::ostringstream o;
  o << dataDescription->GetTimeStep();

  for (int i=0; i<8; i++)
    {
    CppPipeMultiCamera8::Direction d = static_cast<CppPipeMultiCamera8::Direction>(i);
    mc.SetCameraDirection(d);
    mc.InitializeProxy(this->render);
    if(ims.GetImageFormat() == "png")
      {
      std::string filepath = ims.GetBaseDirectoryName() + "/" + base + "."
                             + std::string(mc.DirectionStrings[i]) + "." + o.str() + ".png"; 
      this->render->WriteImage(filepath.c_str(), "vtkPNGWriter");
      }
    else
      {
      std::string filepath = ims.GetBaseDirectoryName() + "/"  + base + "."
                             + std::string(mc.DirectionStrings[i]) + "." + o.str() + ".jpg"; 
      this->render->WriteImage(filepath.c_str(), "vtkJPEGWriter");
      }
    }
}

void CppPipeProxies::RenderCamera(CppPipeImageSet& ims,
                                   CppPipeCamera& c,
                                   CppPipeRepresentation& phrep,
                                   const std::string& image_basename,
                                   vtkSMRepresentationProxy* representation)
{
  ims.InitializeProxy(this->render);
  ims.InitializeProxy(representation);
  phrep.InitializeProxy(representation);

  std::string base = ims.GetImageBaseName();
  if(base.empty())
    base = image_basename;

  std::ostringstream o;
  o << dataDescription->GetTimeStep();

  c.InitializeProxy(this->render);

  if(ims.GetImageFormat() == "png")
    {
    std::string filepath = ims.GetBaseDirectoryName() + "/" + base + "."
                           + o.str() + ".png";
    this->render->WriteImage(filepath.c_str(), "vtkPNGWriter");
    }
  else
    {
    std::string filepath = ims.GetBaseDirectoryName() + "/"  + base + "."
                           + o.str() + ".jpg";
    this->render->WriteImage(filepath.c_str(), "vtkJPEGWriter");
    }
}

void CppPipeProxies::GetGlobalDataBoundsParallel(vtkDoubleArray* globalBounds)
{
  double* localDataBounds = this->producer->GetDataInformation()->GetBounds();
  vtkDoubleArray* localarray = vtkDoubleArray::New();
  localarray->SetNumberOfTuples(6);
  localarray->SetValue(0, -localDataBounds[0]);
  localarray->SetValue(1,  localDataBounds[1]);
  localarray->SetValue(2, -localDataBounds[2]);
  localarray->SetValue(3,  localDataBounds[3]);
  localarray->SetValue(4, -localDataBounds[4]);
  localarray->SetValue(5,  localDataBounds[5]);

  globalBounds->SetNumberOfTuples(6);

  vtkMultiProcessController::GetGlobalController()->AllReduce(localarray, globalBounds, 0);

  globalBounds->SetTuple1(0, -globalBounds->GetTuple1(0));
  globalBounds->SetTuple1(2, -globalBounds->GetTuple1(2));
  globalBounds->SetTuple1(4, -globalBounds->GetTuple1(4));
 
  localarray->Delete();
}

vtkSMSourceProxy* CppPipeProxies::SetupContour(const Json::Value& json,
                                                vtkSMSourceProxy* p)
{
  vtkSMProxyManager* proxyManager = vtkSMProxyManager::GetProxyManager();
  vtkSMSessionProxyManager* sessionProxyManager =
    proxyManager->GetActiveSessionProxyManager();

  CppPipeContour pc(json);

  if(!this->contour)
    this->contour = vtkSMSourceProxy::SafeDownCast(
                       sessionProxyManager->NewProxy("filters", "Contour"));

  vtkSMPropertyHelper(this->contour, "Input").Add(p);
  pc.InitializeProxy(this->contour);
  return this->contour;
}

vtkSMSourceProxy* CppPipeProxies::SetupThreshold(const Json::Value& json,
                                                  vtkSMSourceProxy* p)
{
  vtkSMProxyManager* proxyManager = vtkSMProxyManager::GetProxyManager();
  vtkSMSessionProxyManager* sessionProxyManager =
    proxyManager->GetActiveSessionProxyManager();

  CppPipeThreshold pt(json);

  if(!this->threshold)
    this->threshold = vtkSMSourceProxy::SafeDownCast(
                         sessionProxyManager->NewProxy("filters", "Threshold"));

  vtkSMPropertyHelper(this->threshold, "Input").Add(p);
  pt.InitializeProxy(this->threshold);
  return this->threshold;
}

vtkSMSourceProxy* CppPipeProxies::SetupBoxClip(const Json::Value& json,
                                                vtkSMSourceProxy* p)
{
  vtkSMProxyManager* proxyManager = vtkSMProxyManager::GetProxyManager();
  vtkSMSessionProxyManager* sessionProxyManager =
    proxyManager->GetActiveSessionProxyManager();

  CppPipeBoxClip bc(json,
                     this->globalBounds);
  if(!this->clip)
    this->clip = vtkSMSourceProxy::SafeDownCast(
                    sessionProxyManager->NewProxy("filters", "Clip"));
  if(!this->box)
    this->box = sessionProxyManager->NewProxy("implicit_functions", "Box");

  vtkSMPropertyHelper(this->clip, "Input").Add(p);
  vtkSMPropertyHelper(this->clip, "ClipFunction").Add(this->box);
  this->clip->UpdatePropertyInformation();
  bc.InitializeProxy(this->clip);
  return this->clip;
}

vtkSMSourceProxy* CppPipeProxies::SetupClip(const Json::Value& json,
                                             vtkSMSourceProxy* p)
{
  vtkSMProxyManager* proxyManager = vtkSMProxyManager::GetProxyManager();
  vtkSMSessionProxyManager* sessionProxyManager =
    proxyManager->GetActiveSessionProxyManager();

  CppPipeClip pc(json,
                  this->globalBounds);
  if(!this->clip)
    this->clip = vtkSMSourceProxy::SafeDownCast(
                    sessionProxyManager->NewProxy("filters", "Clip"));
  if(!this->clip_plane)
    this->clip_plane = sessionProxyManager->NewProxy("implicit_functions", "Plane");

  vtkSMPropertyHelper(this->clip, "Input").Add(p);
  vtkSMPropertyHelper(this->clip, "ClipFunction").Add(this->clip_plane);
  this->clip->UpdatePropertyInformation();
  pc.InitializeProxy(this->clip);
  return this->clip;
}

vtkSMSourceProxy* CppPipeProxies::SetupSlice(const Json::Value& json,
                                              vtkSMSourceProxy* p)
{
  vtkSMProxyManager* proxyManager = vtkSMProxyManager::GetProxyManager();
  vtkSMSessionProxyManager* sessionProxyManager =
    proxyManager->GetActiveSessionProxyManager();

  CppPipeSlice ps(json,
                   this->globalBounds);
  if(!this->slice)
    this->slice = vtkSMSourceProxy::SafeDownCast(
                    sessionProxyManager->NewProxy("filters", "Cut"));
  if(!this->slice_plane)
    this->slice_plane = sessionProxyManager->NewProxy("implicit_functions", "Plane");

  vtkSMPropertyHelper(this->slice, "Input").Add(p);
  vtkSMPropertyHelper(this->slice, "CutFunction").Add(this->slice_plane);
  this->slice->UpdatePropertyInformation();
  ps.InitializeProxy(this->slice);
  return this->slice;
}

vtkSMRepresentationProxy* CppPipeProxies::SetupOperations(const std::string& ims_name,
                                                           const Json::Value& json)
{
  vtkSMSourceProxy* p = this->producer;

  Json::Value op = json["imageset blocks"][ims_name];
  while(op.isMember("operation") ||
        op.isMember("input") )
    {
    if(op.isMember("operation"))
      op = json["operation blocks"][op["operation"].asString()];  
    else if(op.isMember("input"))
      op = json["operation blocks"][op["input"].asString()];  

    if(op["type"].asString() == "clip")
      p = this->SetupClip(op,
                          p);
    else if(op["type"].asString() == "slice")
      p = this->SetupSlice(op,
                           p);
    else if(op["type"].asString() == "threshold")
      p = this->SetupThreshold(op,
                               p);
    else if(op["type"].asString() == "boxclip")
      p = this->SetupBoxClip(op,
                             p);
    else if(op["type"].asString() == "contour")
      p = this->SetupContour(op,
                             p);
    }

  if(json["imageset blocks"][ims_name].isMember("threshold"))
    {
    Json::Value t;
    if(json["imageset blocks"][ims_name]["threshold"][0].asString() == "scalar") 
      {
      t["variable scalar"] = json["imageset blocks"][ims_name]["threshold"][1].asString();
      if(json["imageset blocks"][ims_name]["threshold"][2].asString() == "keep above")
        t["keep above"] = json["imageset blocks"][ims_name]["threshold"][3][0];
      else if(json["imageset blocks"][ims_name]["threshold"][2].asString() == "keep below")
        t["keep below"] = json["imageset blocks"][ims_name]["threshold"][3][0];
      else if(json["imageset blocks"][ims_name]["threshold"][2].asString() == "keep between")
        {
        t["keep between"] = Json::Value(Json::arrayValue);
        t["keep between"].append(json["imageset blocks"][ims_name]["threshold"][3][0]);
        t["keep between"].append(json["imageset blocks"][ims_name]["threshold"][3][1]);
        }
      
      p = this->SetupThreshold(t,
                               p);
      }
    }

  if(json["imageset blocks"][ims_name].isMember("box clip"))
    {
    Json::Value c;
    c["center at absolute point"] = json["imageset blocks"][ims_name]["box clip"][0];
    c["absolute extents"] = json["imageset blocks"][ims_name]["box clip"][1];
    if(json["imageset blocks"][ims_name]["box clip"][2].asString() == "keep inside")
      c["keep inside box"] = true;
    else
      c["keep inside box"] = false;

    p = this->SetupBoxClip(c,
                           p);
    }

  if(json["imageset blocks"][ims_name].isMember("clip"))
    {
    Json::Value c;
    c["absolute point on plane"] = json["imageset blocks"][ims_name]["clip"][0];
    c["plane normal"] = json["imageset blocks"][ims_name]["clip"][1];

    p = this->SetupClip(c,
                        p);
    }

  if(json["imageset blocks"][ims_name].isMember("contour"))
    {
    Json::Value c;
    if(json["imageset blocks"][ims_name]["contour"][0].asString() == "scalar") 
      {
      c["variable scalar"] = json["imageset blocks"][ims_name]["contour"][1].asString();
      if(json["imageset blocks"][ims_name]["contour"][2].asString() == "value list")
        {
        c["contour value"] = Json::Value(Json::arrayValue);
        for(int i=0;i<json["imageset blocks"][ims_name]["contour"][3].size();i++)
          c["contour value"].append(json["imageset blocks"][ims_name]["contour"][3][i].asDouble());
        }
      else if(json["imageset blocks"][ims_name]["contour"][2].asString() == "value sequence")
        {
        c["contour value sequence"] = Json::Value(Json::arrayValue);
        for(int i=0;i<json["imageset blocks"][ims_name]["contour"][3].size();i++)
          c["contour value sequence"].append(json["imageset blocks"][ims_name]["contour"][3][i].asDouble());
        }

      p = this->SetupContour(c,
                             p);
      }
    }

  if(json["imageset blocks"][ims_name].isMember("slice"))
    {
    Json::Value s;
    s["absolute point on plane"] = json["imageset blocks"][ims_name]["slice"][0];
    s["plane normal"] = json["imageset blocks"][ims_name]["slice"][1];

    p = this->SetupSlice(s,
                         p);
    }
  
  vtkSMRepresentationProxy* representation = render->CreateDefaultRepresentation(p, 0);
  vtkSMPropertyHelper(representation, "Input").Add(p);
  vtkSMPropertyHelper(this->render, "Representations").Add(representation);
  return representation;
}
