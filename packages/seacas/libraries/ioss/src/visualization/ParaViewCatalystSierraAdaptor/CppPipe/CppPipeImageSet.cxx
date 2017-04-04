
#include <vtkSMRenderViewProxy.h>
#include <vtkSMRepresentationProxy.h>
#include <vtkSMPropertyHelper.h>
#include "CppPipeImageSet.h"

CppPipeImageSet::CppPipeImageSet() : CppPipeObject()
{
  this->SetDefaults();
}

CppPipeImageSet::CppPipeImageSet(const Json::Value& jv) : CppPipeObject(jv)
{
  this->SetDefaults();
}

std::string CppPipeImageSet::GetImageFormat() const
{
  if(this->settings.isMember("image format"))
    return(this->settings["image format"].asString());
  else
    return(this->defaults["image format"].asString());
}

std::string CppPipeImageSet::GetBaseDirectoryName() const
{
  if(this->settings.isMember("image basedirectory"))
    return(this->settings["image basedirectory"].asString());
  else
    return(this->defaults["image basedirectory"].asString());
}

std::string CppPipeImageSet::GetImageBaseName() const
{
  if(this->settings.isMember("image basename"))
    return(this->settings["image basename"].asString());
  else
    return(this->defaults["image basename"].asString());
}

void CppPipeImageSet::SetDefaults()
{
  this->defaults["image basedirectory"] = "CatalystOutput";
  this->defaults["image basename"] = "";
  this->defaults["image format"] = "png";
  this->defaults["image digit count"] = 0;
  this->defaults["image size"] = Json::Value(Json::arrayValue);
  this->defaults["image size"].append(1920);
  this->defaults["image size"].append(1080);
  this->defaults["background color"] = Json::Value(Json::arrayValue);
  this->defaults["background color"].append(0.32);
  this->defaults["background color"].append(0.34);
  this->defaults["background color"].append(0.43);
  this->defaults["surface color"] = Json::Value(Json::arrayValue);
  this->defaults["surface color"].append(1.0);
  this->defaults["surface color"].append(1.0);
  this->defaults["surface color"].append(1.0);
  this->defaults["edge color"] = Json::Value(Json::arrayValue);
  this->defaults["edge color"].append(0.0);
  this->defaults["edge color"].append(0.0);
  this->defaults["edge color"].append(0.5);
  this->defaults["text color"] = Json::Value(Json::arrayValue);
  this->defaults["text color"].append(1.0);
  this->defaults["text color"].append(1.0);
  this->defaults["text color"].append(1.0);
  this->defaults["axes color"] = Json::Value(Json::arrayValue);
  this->defaults["axes color"].append(1.0);
  this->defaults["axes color"].append(1.0);
  this->defaults["axes color"].append(1.0);
}

void CppPipeImageSet::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMRenderViewProxy::SafeDownCast(proxy) ) 
    {
    vtkSMRenderViewProxy* rp = vtkSMRenderViewProxy::SafeDownCast(proxy);

    if(this->settings.isMember("image size"))
      {
      vtkSMPropertyHelper(rp, "ViewSize").Set(0, this->settings["image size"][0].asInt());
      vtkSMPropertyHelper(rp, "ViewSize").Set(1, this->settings["image size"][1].asInt());
      }
    else
      {
      vtkSMPropertyHelper(rp, "ViewSize").Set(0, this->defaults["image size"][0].asInt());
      vtkSMPropertyHelper(rp, "ViewSize").Set(1, this->defaults["image size"][1].asInt());
      }

    if(this->settings.isMember("background color"))
      {
      vtkSMPropertyHelper(rp, "Background").Set(0, this->settings["background color"][0].asDouble());
      vtkSMPropertyHelper(rp, "Background").Set(1, this->settings["background color"][1].asDouble());
      vtkSMPropertyHelper(rp, "Background").Set(2, this->settings["background color"][2].asDouble());
      }
    else
      {
      vtkSMPropertyHelper(rp, "Background").Set(0, this->defaults["background color"][0].asDouble());
      vtkSMPropertyHelper(rp, "Background").Set(1, this->defaults["background color"][1].asDouble());
      vtkSMPropertyHelper(rp, "Background").Set(2, this->defaults["background color"][2].asDouble());
      }

    if(this->settings.isMember("text color"))
      {
      vtkSMPropertyHelper(rp, "OrientationAxesLabelColor").Set(0, this->settings["text color"][0].asDouble());
      vtkSMPropertyHelper(rp, "OrientationAxesLabelColor").Set(1, this->settings["text color"][1].asDouble());
      vtkSMPropertyHelper(rp, "OrientationAxesLabelColor").Set(2, this->settings["text color"][2].asDouble());
      }
    else
      {
      vtkSMPropertyHelper(rp, "OrientationAxesLabelColor").Set(0, this->defaults["text color"][0].asDouble());
      vtkSMPropertyHelper(rp, "OrientationAxesLabelColor").Set(1, this->defaults["text color"][1].asDouble());
      vtkSMPropertyHelper(rp, "OrientationAxesLabelColor").Set(2, this->defaults["text color"][2].asDouble());
      }

    rp->UpdatePropertyInformation();
    rp->UpdateVTKObjects();
    }

  if( proxy &&
      vtkSMRepresentationProxy::SafeDownCast(proxy) ) 
    {
    vtkSMRepresentationProxy* rep = vtkSMRepresentationProxy::SafeDownCast(proxy);

    if(this->settings.isMember("surface color"))
      {
      vtkSMPropertyHelper(rep, "DiffuseColor").Set(0, this->settings["surface color"][0].asDouble());
      vtkSMPropertyHelper(rep, "DiffuseColor").Set(1, this->settings["surface color"][1].asDouble());
      vtkSMPropertyHelper(rep, "DiffuseColor").Set(2, this->settings["surface color"][2].asDouble());
      }
    else
      {
      vtkSMPropertyHelper(rep, "DiffuseColor").Set(0, this->defaults["surface color"][0].asDouble());
      vtkSMPropertyHelper(rep, "DiffuseColor").Set(1, this->defaults["surface color"][1].asDouble());
      vtkSMPropertyHelper(rep, "DiffuseColor").Set(2, this->defaults["surface color"][2].asDouble());
      }

    if(this->settings.isMember("edge color"))
      {
      vtkSMPropertyHelper(rep, "EdgeColor").Set(0, this->settings["edge color"][0].asDouble());
      vtkSMPropertyHelper(rep, "EdgeColor").Set(1, this->settings["edge color"][1].asDouble());
      vtkSMPropertyHelper(rep, "EdgeColor").Set(2, this->settings["edge color"][2].asDouble());
      }
    else
      {
      vtkSMPropertyHelper(rep, "EdgeColor").Set(0, this->defaults["edge color"][0].asDouble());
      vtkSMPropertyHelper(rep, "EdgeColor").Set(1, this->defaults["edge color"][1].asDouble());
      vtkSMPropertyHelper(rep, "EdgeColor").Set(2, this->defaults["edge color"][2].asDouble());
      }

    if(this->settings.isMember("axes color"))
      {
      vtkSMPropertyHelper(rep, "CubeAxesColor").Set(0, this->settings["axes color"][0].asDouble());
      vtkSMPropertyHelper(rep, "CubeAxesColor").Set(1, this->settings["axes color"][1].asDouble());
      vtkSMPropertyHelper(rep, "CubeAxesColor").Set(2, this->settings["axes color"][2].asDouble());
      }
    else
      {
      vtkSMPropertyHelper(rep, "CubeAxesColor").Set(0, this->defaults["axes color"][0].asDouble());
      vtkSMPropertyHelper(rep, "CubeAxesColor").Set(1, this->defaults["axes color"][1].asDouble());
      vtkSMPropertyHelper(rep, "CubeAxesColor").Set(2, this->defaults["axes color"][2].asDouble());
      }

    rep->UpdatePropertyInformation();
    rep->UpdateVTKObjects();
    }
}
