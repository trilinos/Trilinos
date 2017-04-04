
#include "CppPipeThreshold.h"
#include <vtkSMSourceProxy.h>
#include <vtkSMPropertyHelper.h>
#include <vtkMath.h>

CppPipeThreshold::CppPipeThreshold() : CppPipeObject()
{
  this->SetDefaults();
}

CppPipeThreshold::CppPipeThreshold(const Json::Value& jv) : CppPipeObject(jv)
{
  this->SetDefaults();
}

CppPipeThreshold::~CppPipeThreshold()
{

}

void CppPipeThreshold::SetDefaults()
{

}

void CppPipeThreshold::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMSourceProxy::SafeDownCast(proxy) ) 
    {
    vtkSMSourceProxy* sp = vtkSMSourceProxy::SafeDownCast(proxy);

    if(this->settings.isMember("variable scalar"))
      {
      vtkSMPropertyHelper(sp, "SelectInputScalars").Set(4, this->settings["variable scalar"].asString().c_str());
      }

    if(this->settings.isMember("keep between"))
      {
      vtkSMPropertyHelper(sp, "ThresholdBetween").Set(0, this->settings["keep between"][0].asDouble());
      vtkSMPropertyHelper(sp, "ThresholdBetween").Set(1, this->settings["keep between"][1].asDouble());
      }
    else if(this->settings.isMember("keep above"))
      {
      vtkSMPropertyHelper(sp, "ThresholdBetween").Set(0, this->settings["keep above"].asDouble());
      vtkSMPropertyHelper(sp, "ThresholdBetween").Set(1, vtkMath::Inf());
      }
    else if(this->settings.isMember("keep below"))
      {
      vtkSMPropertyHelper(sp, "ThresholdBetween").Set(0, vtkMath::NegInf());
      vtkSMPropertyHelper(sp, "ThresholdBetween").Set(1, this->settings["keep below"].asDouble());
      }

    sp->UpdatePropertyInformation();
    sp->UpdateVTKObjects();
    }
}
