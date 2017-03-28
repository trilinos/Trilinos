#ifndef CPPPIPECLIP_H
#define CPPPIPECLIP_H

#include "CppPipeObject.h"

class vtkDoubleArray;

class CppPipeClip : public CppPipeObject
{
public:
  CppPipeClip(vtkDoubleArray* globalBounds);

  CppPipeClip(const Json::Value& jv,
                vtkDoubleArray* globalBounds);

  virtual ~CppPipeClip();

  void InitializeProxy(vtkSMProxy* p);

private:
  void SetDefaults();
  void SetGlobalBounds(vtkDoubleArray* globalBounds);
  vtkDoubleArray* gb;
  double inFocalPoint[3];
  double intX;
  double intY;
  double intZ;
};
#endif
