
#include "CppPipeBoxClip.h"
#include <vtkSMSourceProxy.h>
#include <vtkDoubleArray.h>
#include <vtkSMPropertyHelper.h>

CppPipeBoxClip::CppPipeBoxClip(vtkDoubleArray* globalBounds) : CppPipeObject()
{
  this->SetDefaults();
  this->SetGlobalBounds(globalBounds);
}

CppPipeBoxClip::CppPipeBoxClip(const Json::Value& jv,
                           vtkDoubleArray* globalBounds) : CppPipeObject(jv)
{
  this->SetDefaults();
  this->SetGlobalBounds(globalBounds);
}

CppPipeBoxClip::~CppPipeBoxClip()
{
  if(this->gb)
    this->gb->Delete();
}

void CppPipeBoxClip::SetDefaults()
{
  this->defaults["center at relative point"] = Json::Value(Json::arrayValue);
  this->defaults["center at relative point"].append(0.0);
  this->defaults["center at relative point"].append(0.0);
  this->defaults["center at relative point"].append(0.0);

  this->defaults["relative extents"] = Json::Value(Json::arrayValue);
  this->defaults["relative extents"].append(0.5);
  this->defaults["relative extents"].append(0.5);
  this->defaults["relative extents"].append(0.5);

  this->defaults["rotations"] = Json::Value(Json::arrayValue);
  this->defaults["rotations"].append(0.0);
  this->defaults["rotations"].append(0.0);
  this->defaults["rotations"].append(0.0);

  this->defaults["keep inside box"] = true;
  this->defaults["cut type"] = "smooth";
}

void CppPipeBoxClip::SetGlobalBounds(vtkDoubleArray* globalBounds)
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

void CppPipeBoxClip::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMSourceProxy::SafeDownCast(proxy) ) 
    {
    vtkSMSourceProxy* sp = vtkSMSourceProxy::SafeDownCast(proxy);

    vtkSMProperty* boxProperty = sp->GetProperty("ClipFunction");
    vtkSMPropertyHelper boxPropertyHelper(boxProperty);
    vtkSMProxy* boxProxy = boxPropertyHelper.GetAsProxy();

    double boxCenterPoint[3];
    bool absolute_center_point = false;
    if(this->settings.isMember("center at absolute point"))
      {
      boxCenterPoint[0] = this->settings["center at absolute point"][0].asDouble();
      boxCenterPoint[1] = this->settings["center at absolute point"][1].asDouble();
      boxCenterPoint[2] = this->settings["center at absolute point"][2].asDouble();
      absolute_center_point = true;
      }

    if(this->settings.isMember("center at relative point"))
      {
      boxCenterPoint[0] = this->inFocalPoint[0] + this->settings["center at relative point"][0].asDouble()*this->intX;
      boxCenterPoint[1] = this->inFocalPoint[1] + this->settings["center at relative point"][1].asDouble()*this->intY;
      boxCenterPoint[2] = this->inFocalPoint[2] + this->settings["center at relative point"][2].asDouble()*this->intZ;
      }
    else if(!absolute_center_point)
      {
      boxCenterPoint[0] = this->inFocalPoint[0] + this->defaults["center at relative point"][0].asDouble()*this->intX;
      boxCenterPoint[1] = this->inFocalPoint[1] + this->defaults["center at relative point"][1].asDouble()*this->intY;
      boxCenterPoint[2] = this->inFocalPoint[2] + this->defaults["center at relative point"][2].asDouble()*this->intZ;
      }

    vtkSMPropertyHelper(boxProxy, "Position").Set(0, boxCenterPoint[0]);
    vtkSMPropertyHelper(boxProxy, "Position").Set(1, boxCenterPoint[1]);
    vtkSMPropertyHelper(boxProxy, "Position").Set(2, boxCenterPoint[2]);

    double boxExtents[6];
    bool absolute_box_extents = false;
    if(this->settings.isMember("absolute extents"))
      {
      boxExtents[0] = boxCenterPoint[0] - this->settings["absolute extents"][0].asDouble();
      boxExtents[1] = boxCenterPoint[0] + this->settings["absolute extents"][0].asDouble();

      boxExtents[2] = boxCenterPoint[1] - this->settings["absolute extents"][1].asDouble();
      boxExtents[3] = boxCenterPoint[1] + this->settings["absolute extents"][1].asDouble();

      boxExtents[4] = boxCenterPoint[2] - this->settings["absolute extents"][2].asDouble();
      boxExtents[5] = boxCenterPoint[2] + this->settings["absolute extents"][2].asDouble();
      absolute_box_extents = true;
      }

    if(this->settings.isMember("relative extents"))
      {
      boxExtents[0] = this->inFocalPoint[0] - this->settings["relative extents"][0].asDouble()*this->intX;
      boxExtents[1] = this->inFocalPoint[0] + this->settings["relative extents"][0].asDouble()*this->intX;

      boxExtents[2] = this->inFocalPoint[1] - this->settings["relative extents"][1].asDouble()*this->intY;
      boxExtents[3] = this->inFocalPoint[1] + this->settings["relative extents"][1].asDouble()*this->intY;

      boxExtents[4] = this->inFocalPoint[2] - this->settings["relative extents"][2].asDouble()*this->intZ;
      boxExtents[5] = this->inFocalPoint[2] + this->settings["relative extents"][2].asDouble()*this->intZ;
      }
    else if(!absolute_box_extents)
      {
      boxExtents[0] = this->inFocalPoint[0] - this->defaults["relative extents"][0].asDouble()*this->intX;
      boxExtents[1] = this->inFocalPoint[0] + this->defaults["relative extents"][0].asDouble()*this->intX;

      boxExtents[2] = this->inFocalPoint[1] - this->defaults["relative extents"][1].asDouble()*this->intY;
      boxExtents[3] = this->inFocalPoint[1] + this->defaults["relative extents"][1].asDouble()*this->intY;

      boxExtents[4] = this->inFocalPoint[2] - this->defaults["relative extents"][2].asDouble()*this->intZ;
      boxExtents[5] = this->inFocalPoint[2] + this->defaults["relative extents"][2].asDouble()*this->intZ;
      }

    vtkSMPropertyHelper(boxProxy, "Bounds").Set(0, boxExtents[0]);
    vtkSMPropertyHelper(boxProxy, "Bounds").Set(1, boxExtents[1]);
    vtkSMPropertyHelper(boxProxy, "Bounds").Set(2, boxExtents[2]);
    vtkSMPropertyHelper(boxProxy, "Bounds").Set(3, boxExtents[3]);
    vtkSMPropertyHelper(boxProxy, "Bounds").Set(4, boxExtents[4]);
    vtkSMPropertyHelper(boxProxy, "Bounds").Set(5, boxExtents[5]);

    double Rotations[3];
    if(this->settings.isMember("rotations"))
      {
      Rotations[0] = this->settings["rotations"][0].asDouble();
      Rotations[1] = this->settings["rotations"][1].asDouble();
      Rotations[2] = this->settings["rotations"][2].asDouble();
      }
    else
      {
      Rotations[0] = this->defaults["rotations"][0].asDouble();
      Rotations[1] = this->defaults["rotations"][1].asDouble();
      Rotations[2] = this->defaults["rotations"][2].asDouble();
      }

    vtkSMPropertyHelper(boxProxy, "Rotation").Set(0, Rotations[0]);
    vtkSMPropertyHelper(boxProxy, "Rotation").Set(1, Rotations[1]);
    vtkSMPropertyHelper(boxProxy, "Rotation").Set(2, Rotations[2]);

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

    if(this->settings.isMember("keep inside box"))
      {
      if(this->settings["keep inside box"].asBool())
        vtkSMPropertyHelper(sp, "InsideOut").Set(1);
      else
        vtkSMPropertyHelper(sp, "InsideOut").Set(0);
      }
    else
      {
      if(this->defaults["keep inside box"].asBool())
        vtkSMPropertyHelper(sp, "InsideOut").Set(1);
      else
        vtkSMPropertyHelper(sp, "InsideOut").Set(0);
      }

    boxProxy->UpdatePropertyInformation();
    boxProxy->UpdateVTKObjects();

    sp->UpdatePropertyInformation();
    sp->UpdateVTKObjects();
    }
}
