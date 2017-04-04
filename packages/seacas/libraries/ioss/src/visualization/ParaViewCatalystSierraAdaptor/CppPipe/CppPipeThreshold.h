#ifndef CPPPIPETHRESHOLD_H
#define CPPPIPETHRESHOLD_H

#include "CppPipeObject.h"

class CppPipeThreshold : public CppPipeObject
{
public:
  CppPipeThreshold();

  CppPipeThreshold(const Json::Value& jv);

  virtual ~CppPipeThreshold();

  void InitializeProxy(vtkSMProxy* p);

private:
  void SetDefaults();
};
#endif
