#ifndef CPPPIPEBOXCLIP_H
#define CPPPIPEBOXCLIP_H

#include "CppPipeObject.h"

class vtkDoubleArray;

class CppPipeBoxClip : public CppPipeObject
{
public:
  CppPipeBoxClip(vtkDoubleArray* globalBounds);

  CppPipeBoxClip(const Json::Value& jv,
                vtkDoubleArray* globalBounds);

  virtual ~CppPipeBoxClip();

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
