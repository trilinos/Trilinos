
#include "CppPipeSlice.h"
#include <vtkSMSourceProxy.h>
#include <vtkDoubleArray.h>
#include <vtkSMPropertyHelper.h>

CppPipeSlice::CppPipeSlice(vtkDoubleArray* globalBounds) : CppPipeObject()
{
  this->SetDefaults();
  this->SetGlobalBounds(globalBounds);
}

CppPipeSlice::CppPipeSlice(const Json::Value& jv,
                             vtkDoubleArray* globalBounds) : CppPipeObject(jv)
{
  this->SetDefaults();
  this->SetGlobalBounds(globalBounds);
}

CppPipeSlice::~CppPipeSlice()
{
  if(this->gb)
    this->gb->Delete();
}

void CppPipeSlice::SetDefaults()
{
  this->defaults["relative point on plane"] = Json::Value(Json::arrayValue);
  this->defaults["relative point on plane"].append(0.0);
  this->defaults["relative point on plane"].append(0.0);
  this->defaults["relative point on plane"].append(0.0);

  this->defaults["plane normal"] = Json::Value(Json::arrayValue);
  this->defaults["plane normal"].append(0.0);
  this->defaults["plane normal"].append(1.0);
  this->defaults["plane normal"].append(0.0);

  this->defaults["cut type"] = "smooth";
}

void CppPipeSlice::SetGlobalBounds(vtkDoubleArray* globalBounds)
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

void CppPipeSlice::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMSourceProxy::SafeDownCast(proxy) ) 
    {
    vtkSMSourceProxy* sp = vtkSMSourceProxy::SafeDownCast(proxy);

    vtkSMProperty* planeProperty = sp->GetProperty("CutFunction");
    vtkSMPropertyHelper planePropertyHelper(planeProperty);
    vtkSMProxy* planeProxy = planePropertyHelper.GetAsProxy();

    double planePoint[3];
    bool absolute_point = false;
    if(this->settings.isMember("absolute point on plane"))
      {
      planePoint[0] = this->settings["absolute point on plane"][0].asDouble();
      planePoint[1] = this->settings["absolute point on plane"][1].asDouble();
      planePoint[2] = this->settings["absolute point on plane"][2].asDouble();
      absolute_point = true;
      }

    if(this->settings.isMember("relative point on plane"))
      {
      planePoint[0] = this->inFocalPoint[0] + this->settings["relative point on plane"][0].asDouble()*this->intX;
      planePoint[1] = this->inFocalPoint[1] + this->settings["relative point on plane"][1].asDouble()*this->intY;
      planePoint[2] = this->inFocalPoint[2] + this->settings["relative point on plane"][2].asDouble()*this->intZ;
      }
    else if(!absolute_point)
      {
      planePoint[0] = this->inFocalPoint[0] + this->defaults["relative point on plane"][0].asDouble()*this->intX;
      planePoint[1] = this->inFocalPoint[1] + this->defaults["relative point on plane"][1].asDouble()*this->intY;
      planePoint[2] = this->inFocalPoint[2] + this->defaults["relative point on plane"][2].asDouble()*this->intZ;
      }

    vtkSMPropertyHelper(planeProxy, "Origin").Set(0, planePoint[0]);
    vtkSMPropertyHelper(planeProxy, "Origin").Set(1, planePoint[1]);
    vtkSMPropertyHelper(planeProxy, "Origin").Set(2, planePoint[2]);

    double planeNormal[3];
    if(this->settings.isMember("plane normal"))
      {
      planeNormal[0] = this->settings["plane normal"][0].asDouble();
      planeNormal[1] = this->settings["plane normal"][1].asDouble();
      planeNormal[2] = this->settings["plane normal"][2].asDouble();
      }
    else
      {
      planeNormal[0] = this->defaults["plane normal"][0].asDouble();
      planeNormal[1] = this->defaults["plane normal"][1].asDouble();
      planeNormal[2] = this->defaults["plane normal"][2].asDouble();
      }

    vtkSMPropertyHelper(planeProxy, "Normal").Set(0, planeNormal[0]);
    vtkSMPropertyHelper(planeProxy, "Normal").Set(1, planeNormal[1]);
    vtkSMPropertyHelper(planeProxy, "Normal").Set(2, planeNormal[2]);

    if(this->settings.isMember("cut type"))
      {
      if(this->settings["cut type"].asString() == "crinkle")
        vtkSMPropertyHelper(sp, "PreserveInputCells").Set(1);
      else
        vtkSMPropertyHelper(sp, "PreserveInputCells").Set(0);
      }
    else
      {
      if(this->defaults["cut type"].asString() == "crinkle")
        vtkSMPropertyHelper(sp, "PreserveInputCells").Set(1);
      else
        vtkSMPropertyHelper(sp, "PreserveInputCells").Set(0);
      }

    vtkSMPropertyHelper(sp, "Triangulate the slice").Set(0);

    planeProxy->UpdatePropertyInformation();
    planeProxy->UpdateVTKObjects();

    sp->UpdatePropertyInformation();
    sp->UpdateVTKObjects();
    }
}
