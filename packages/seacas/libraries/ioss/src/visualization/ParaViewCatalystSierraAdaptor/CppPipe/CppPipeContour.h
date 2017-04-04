#ifndef CPPPIPECONTOUR_H
#define CPPPIPECONTOUR_H

#include "CppPipeObject.h"

class CppPipeContour : public CppPipeObject
{
public:
  CppPipeContour();

  CppPipeContour(const Json::Value& jv);

  virtual ~CppPipeContour();

  void InitializeProxy(vtkSMProxy* p);

private:
  void SetDefaults();
};
#endif
