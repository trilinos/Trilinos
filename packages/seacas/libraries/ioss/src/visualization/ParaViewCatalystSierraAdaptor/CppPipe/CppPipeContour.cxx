
#include "CppPipeContour.h"
#include <vtkSMSourceProxy.h>
#include <vtkSMPropertyHelper.h>

CppPipeContour::CppPipeContour() : CppPipeObject()
{
  this->SetDefaults();
}

CppPipeContour::CppPipeContour(const Json::Value& jv) : CppPipeObject(jv)
{
  this->SetDefaults();
}

CppPipeContour::~CppPipeContour()
{

}

void CppPipeContour::SetDefaults()
{

}

void CppPipeContour::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMSourceProxy::SafeDownCast(proxy) ) 
    {
    vtkSMSourceProxy* sp = vtkSMSourceProxy::SafeDownCast(proxy);

    if(this->settings.isMember("variable scalar"))
      {
      vtkSMPropertyHelper(sp, "SelectInputScalars").Set(4, this->settings["variable scalar"].asString().c_str());
      }

    if(this->settings.isMember("contour value"))
      {
      double* d = new double[this->settings["contour value"].size()];
      for(int i=0;i<this->settings["contour value"].size();i++)
        d[i] = this->settings["contour value"][i].asDouble();
      
      vtkSMPropertyHelper(sp, "ContourValues").Set(d, this->settings["contour value"].size());
      delete [] d;
      }

    if(this->settings.isMember("contour value sequence"))
      {
      double start = this->settings["contour value sequence"][0].asDouble();
      double step = this->settings["contour value sequence"][1].asDouble();
      double stop = this->settings["contour value sequence"][2].asDouble();
      int n = int((stop-start)/step) + 1;

      double* d = new double[n];
      for(int i=0;i<n;i++)
        d[i] = start + i*step;

      vtkSMPropertyHelper(sp, "ContourValues").Set(d, n);
      delete [] d;
      }

    sp->UpdatePropertyInformation();
    sp->UpdateVTKObjects();
    }
}
