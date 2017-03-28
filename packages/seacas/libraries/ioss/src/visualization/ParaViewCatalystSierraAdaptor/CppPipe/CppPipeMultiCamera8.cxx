
#include <vtkSMRenderViewProxy.h>
#include <vtkSMPropertyHelper.h>
#include <vtkDoubleArray.h>
#include <vtkMultiProcessController.h>
#include <vtkMath.h>
#include "CppPipeMultiCamera8.h"

const char* CppPipeMultiCamera8::DirectionStrings[8] = {"x1", "x2", "y1", "y2",
                                                         "z1", "z2", "xyz1", "xyz2"};

CppPipeMultiCamera8::CppPipeMultiCamera8(vtkDoubleArray* globalBounds) : CppPipeObject()
{
  this->SetGlobalBounds(globalBounds);
  this->SetDefaults();
}

CppPipeMultiCamera8::CppPipeMultiCamera8(const Json::Value& jv,
                                           vtkDoubleArray* globalBounds) : CppPipeObject(jv)
{
  this->SetGlobalBounds(globalBounds);
  this->SetDefaults();
}

CppPipeMultiCamera8::~CppPipeMultiCamera8()
{
  if(this->gb)
    this->gb->Delete();
}

void CppPipeMultiCamera8::SetGlobalBounds(vtkDoubleArray* globalBounds)
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

  if (this->intX >= this->intZ && this->intX >= this->intY)
    this->maxDim = this->intX;
  else if (this->intY >= this->intX && this->intY >= this->intZ)
    this->maxDim = this->intY;
  else if (this->intZ >= this->intX && this->intZ >= this->intY)
    this->maxDim = this->intZ;
}

void CppPipeMultiCamera8::SetDefaults()
{
  this->defaults["up vector"] = Json::Value(Json::arrayValue);
  this->defaults["up vector"].append(0.0);
  this->defaults["up vector"].append(1.0);
  this->defaults["up vector"].append(0.0);

  this->defaults["look at relative point"] = Json::Value(Json::arrayValue);
  this->defaults["look at relative point"].append(0.0);
  this->defaults["look at relative point"].append(0.0);
  this->defaults["look at relative point"].append(0.0);

  this->defaults["look at relative distance"] = 1.0;

  this->defaults["camera fov"] = 30.0;
  this->defaults["image name addon"] = "";
  this->d = X1;
}

void CppPipeMultiCamera8::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMRenderViewProxy::SafeDownCast(proxy) ) 
    {
    vtkSMRenderViewProxy* rv = vtkSMRenderViewProxy::SafeDownCast(proxy);

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

    double lookPoint[3];
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

    double EyePosition[3];
    EyePosition[0] = lookPoint[0];
    EyePosition[1] = lookPoint[1];
    EyePosition[2] = lookPoint[2];

    double mag = sqrt(this->intX*this->intX +
                      this->intY*this->intY +
                      this->intZ*this->intZ);
    double distance = 0.0;
    if(this->settings.isMember("look at absolute distance"))
      distance = this->settings["look at absolute distance"].asDouble();
    else if(this->settings.isMember("look at relative distance"))
      distance = mag*this->settings["look at relative distance"].asDouble();
    else
      distance = 2.0*mag*this->defaults["look at relative distance"].asDouble();

    if(this->d == X1)
      {
      EyePosition[0] += distance;
      double dir[3] = {-1.0, 0.0, 0.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == X2)
      {
      EyePosition[0] -= distance;
      double dir[3] = {1.0, 0.0, 0.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == Y1)
      {
      EyePosition[1] += distance;
      double dir[3] = {0.0, -1.0, 0.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == Y2)
      {
      EyePosition[1] -= distance;
      double dir[3] = {0.0, 1.0, 0.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == Z1)
      {
      EyePosition[2] += distance;
      double dir[3] = {0.0, 0.0, -1.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == Z2)
      {
      EyePosition[2] -= distance;
      double dir[3] = {0.0, 0.0, 1.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == XYZ1)
      {
      EyePosition[0] += distance;
      EyePosition[1] += distance;
      EyePosition[2] += distance;
      double dir[3] = {-1.0, -1.0, -1.0};
      this->CheckParallelVectors(dir, ViewUp);
      }
    else if(this->d == XYZ2)
      {
      EyePosition[0] -= distance;
      EyePosition[1] -= distance;
      EyePosition[2] -= distance;
      double dir[3] = {1.0, 1.0, 1.0};
      this->CheckParallelVectors(dir, ViewUp);
      }

    vtkSMPropertyHelper(rv, "CameraViewUp").Set(0, ViewUp[0]);
    vtkSMPropertyHelper(rv, "CameraViewUp").Set(1, ViewUp[1]);
    vtkSMPropertyHelper(rv, "CameraViewUp").Set(2, ViewUp[2]);

    vtkSMPropertyHelper(rv, "CameraPosition").Set(0, EyePosition[0]);
    vtkSMPropertyHelper(rv, "CameraPosition").Set(1, EyePosition[1]);
    vtkSMPropertyHelper(rv, "CameraPosition").Set(2, EyePosition[2]);

    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(0, lookPoint[0]);
    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(1, lookPoint[1]);
    vtkSMPropertyHelper(rv, "CameraFocalPoint").Set(2, lookPoint[2]);

    rv->UpdatePropertyInformation();
    rv->UpdateVTKObjects();
    }
}

void CppPipeMultiCamera8::SetCameraDirection(Direction d)
{
  this->d = d;
}

void CppPipeMultiCamera8::CheckParallelVectors(double* look_direction, 
                                                double* view_up)
{
  if(fabs(vtkMath::Dot(look_direction, view_up)) > 0.99)
    {
    view_up[0] = 0.0;
    view_up[1] = 0.0;
    view_up[2] = 1.0;
    if(fabs(vtkMath::Dot(look_direction, view_up)) > 0.99)
      {
      view_up[0] = 0.0;
      view_up[1] = 1.0;
      view_up[2] = 0.0;
      }
    }
}
