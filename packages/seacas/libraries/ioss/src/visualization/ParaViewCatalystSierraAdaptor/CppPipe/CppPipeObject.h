#ifndef CPPPIPEOBJECT_H
#define CPPPIPEOBJECT_H

#ifdef USE_VTK_JSON
#include "vtk_jsoncpp.h"
#else
#include <json/json.h>
#endif

#include <vtkSMProxy.h>

class CppPipeObject
{
public:

  CppPipeObject() {}

  CppPipeObject(const Json::Value& jv)
    {
    this->settings = jv;
    }

  virtual ~CppPipeObject() {}

  virtual void InitializeProxy(vtkSMProxy* p) = 0;

  void SetJSON(const Json::Value& jv)
    {
    this->settings = jv;
    }

  const Json::Value& GetDefaults() { return this->defaults; }
  const Json::Value& GetSettings() { return this->settings; }

protected:
  virtual void SetDefaults() = 0;
  Json::Value settings;
  Json::Value defaults;

private:
};
#endif
