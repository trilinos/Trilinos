
#include <vtkSMRenderViewProxy.h>
#include <vtkSMPropertyHelper.h>
#include <vtkDoubleArray.h>
#include <vtkMultiProcessController.h>
#include <vtkMath.h>
#include "CppPipeCamera.h"

CppPipeCamera::CppPipeCamera(vtkDoubleArray* globalBounds) : CppPipeObject()
{
  this->SetGlobalBounds(globalBounds);
  this->SetDefaults();
}

CppPipeCamera::CppPipeCamera(const Json::Value& jv,
                               vtkDoubleArray* globalBounds) : CppPipeObject(jv)
{
  this->SetGlobalBounds(globalBounds);
  this->SetDefaults();
}

CppPipeCamera::~CppPipeCamera()
{
  if(this->gb)
    this->gb->Delete();
}

void CppPipeCamera::SetGlobalBounds(vtkDoubleArray* globalBounds)
{
  this->gb = vtkDoubleArray::New();
  this->gb->SetNumberOfTuples(6);
  if(globalBounds->GetNumberOfTuples() == 6);
    for(int i=0; i<6 ; i++)
      this->gb->SetTuple1(i, globalBounds->GetTuple1(i));

  this->inFocalPoint[0] = (this->gb->GetTuple1(1) + this->gb->GetTuple1(0))/0.5;
  this->inFocalPoint[1] = (this->gb->GetTuple1(3) + this->gb->GetTuple1(2))/0.5;
  this->inFocalPoint[2] = (this->gb->GetTuple1(5) + this->gb->GetTuple1(4))/0.5;

  this->intX = this->gb->GetTuple1(1) - this->gb->GetTuple1(0);
  this->intY = this->gb->GetTuple1(3) - this->gb->GetTuple1(2);
  this->intZ = this->gb->GetTuple1(5) - this->gb->GetTuple1(4);
}

void CppPipeCamera::SetDefaults()
{
  this->defaults["look at relative point"] = Json::Value(Json::arrayValue);
  this->defaults["look at relative point"].append(0.0);
  this->defaults["look at relative point"].append(0.0);
  this->defaults["look at relative point"].append(0.0);

  this->defaults["look direction"] = Json::Value(Json::arrayValue);
  this->defaults["look direction"].append(-1.0);
  this->defaults["look direction"].append(-1.0);
  this->defaults["look direction"].append(-1.0);

  this->defaults["look at relative distance"] = 1.0;

  this->defaults["up vector"] = Json::Value(Json::arrayValue);
  this->defaults["up vector"].append(0.0);
  this->defaults["up vector"].append(1.0);
  this->defaults["up vector"].append(0.0);
  this->defaults["camera fov"] = 30.0;
  this->defaults["image name addon"] = "";
}

void CppPipeCamera::SetLookPoint(double* lookPoint)
{
  if(this->settings.isMember("look at absolute point"))
    {
    lookPoint[0] = this->settings["look at absolute point"][0].asDouble();
    lookPoint[1] = this->settings["look at absolute point"][1].asDouble();
    lookPoint[2] = this->settings["look at absolute point"][2].asDouble();
    }
  else
    {
    if(this->settings.isMember("look at relative point"))
      {
      lookPoint[0] = this->inFocalPoint[0] + this->settings["look at relative point"][0].asDouble()*this->intX;
      lookPoint[1] = this->inFocalPoint[1] + this->settings["look at relative point"][1].asDouble()*this->intY;
      lookPoint[2] = this->inFocalPoint[2] + this->settings["look at relative point"][2].asDouble()*this->intZ;
      }
    else
      {
      lookPoint[0] = this->inFocalPoint[0] + this->defaults["look at relative point"][0].asDouble()*this->intX;
      lookPoint[1] = this->inFocalPoint[1] + this->defaults["look at relative point"][1].asDouble()*this->intY;
      lookPoint[2] = this->inFocalPoint[2] + this->defaults["look at relative point"][2].asDouble()*this->intZ;
      }
    }
}

void CppPipeCamera::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMRenderViewProxy::SafeDownCast(proxy) ) 
    {
    vtkSMRenderViewProxy* rv = vtkSMRenderViewProxy::SafeDownCast(proxy);

    bool hasCameraAt = false;
    double EyePosition[3];
    double lookPoint[3];
    if(this->settings.isMember("camera at absolute point"))
      {
      EyePosition[0] = this->settings["camera at absolute point"][0].asDouble();
      EyePosition[1] = this->settings["camera at absolute point"][1].asDouble();
      EyePosition[2] = this->settings["camera at absolute point"][2].asDouble();
      hasCameraAt = true;
      }

    if(this->settings.isMember("camera at relative point"))
      {
      EyePosition[0] = this->inFocalPoint[0] + this->settings["camera at relative point"][0].asDouble()*this->intX;
      EyePosition[1] = this->inFocalPoint[1] + this->settings["camera at relative point"][1].asDouble()*this->intY;
      EyePosition[2] = this->inFocalPoint[2] + this->settings["camera at relative point"][2].asDouble()*this->intZ;
      hasCameraAt = true;
      }

    if(hasCameraAt)
      {
      if(this->settings.isMember("look direction"))
        {
        lookPoint[0] = EyePosition[0] + this->settings["look direction"][0].asDouble();
        lookPoint[1] = EyePosition[1] + this->settings["look direction"][1].asDouble();
        lookPoint[2] = EyePosition[2] + this->settings["look direction"][2].asDouble();
        }
      else
        this->SetLookPoint(lookPoint);
      }
    else
      {
      this->SetLookPoint(lookPoint);
      double mag = sqrt(this->intX*this->intX + 
                        this->intY*this->intY + 
                        this->intZ*this->intZ);
      double distance = 0.0;
      if(this->settings.isMember("look at absolute distance"))
        distance = this->settings["look at absolute distance"].asDouble();
      else if(this->settings.isMember("look at relative distance"))
        distance = mag*this->settings["look at relative distance"].asDouble();
      else
        distance = mag*this->defaults["look at relative distance"].asDouble();
 
      if(this->settings.isMember("look direction"))
        {
        EyePosition[0] = lookPoint[0] + distance*(-this->settings["look direction"][0].asDouble());
        EyePosition[1] = lookPoint[1] + distance*(-this->settings["look direction"][1].asDouble());
        EyePosition[2] = lookPoint[2] + distance*(-this->settings["look direction"][2].asDouble());
        }
      else
        {
        EyePosition[0] = lookPoint[0] + distance*(-this->defaults["look direction"][0].asDouble());
        EyePosition[1] = lookPoint[1] + distance*(-this->defaults["look direction"][1].asDouble());
        EyePosition[2] = lookPoint[2] + distance*(-this->defaults["look direction"][2].asDouble());
        }
      }

    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(0, lookPoint[0]);
    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(1, lookPoint[1]);
    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(2, lookPoint[2]);

    vtkSMPropertyHelper(rv, "CameraPosition").Set(0, EyePosition[0]);
    vtkSMPropertyHelper(rv, "CameraPosition").Set(1, EyePosition[1]);
    vtkSMPropertyHelper(rv, "CameraPosition").Set(2, EyePosition[2]);

    double ViewAngle;
    if(this->settings.isMember("camera fov"))
      {
      ViewAngle = this->settings["camera fov"].asDouble();
      }
    else
      {
      ViewAngle = this->defaults["camera fov"].asDouble();
      }

   vtkSMPropertyHelper(rv, "CameraViewAngle").Set(ViewAngle);

    double ViewUp[3];
    if(this->settings.isMember("up vector"))
      {
      ViewUp[0] = this->settings["up vector"][0].asDouble();
      ViewUp[1] = this->settings["up vector"][1].asDouble();
      ViewUp[2] = this->settings["up vector"][2].asDouble();
      }
    else
      {
      ViewUp[0] = this->defaults["up vector"][0].asDouble();
      ViewUp[1] = this->defaults["up vector"][1].asDouble();
      ViewUp[2] = this->defaults["up vector"][2].asDouble();
      }

    double inFocalPoint[3];
    inFocalPoint[0] = (this->gb->GetTuple1(1) + this->gb->GetTuple1(0))/0.5;
    inFocalPoint[1] = (this->gb->GetTuple1(3) + this->gb->GetTuple1(2))/0.5;
    inFocalPoint[2] = (this->gb->GetTuple1(5) + this->gb->GetTuple1(4))/0.5;

    double Xdim = this->gb->GetTuple1(1) - this->gb->GetTuple1(0);
    double Ydim = this->gb->GetTuple1(3) - this->gb->GetTuple1(2);
    double Zdim = this->gb->GetTuple1(5) - this->gb->GetTuple1(4);

    double Mdim = 0.0;
    if (Xdim >= Zdim && Xdim >= Ydim)
      Mdim = Xdim;
    else if (Ydim >= Xdim && Ydim >= Zdim)
      Mdim = Ydim;
    else if (Zdim >= Xdim && Zdim >= Ydim)
      Mdim = Zdim;
 
    vtkSMPropertyHelper(rv, "CameraViewUp").Set(0, ViewUp[0]);
    vtkSMPropertyHelper(rv, "CameraViewUp").Set(1, ViewUp[1]);
    vtkSMPropertyHelper(rv, "CameraViewUp").Set(2, ViewUp[2]);

    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(0, inFocalPoint[0]);
    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(1, inFocalPoint[1]);
    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(2, inFocalPoint[2]);

    rv->UpdatePropertyInformation();
    rv->UpdateVTKObjects();
    }
}
