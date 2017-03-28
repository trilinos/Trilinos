#ifndef CPPPIPESLICE_H
#define CPPPIPESLICE_H

#include "CppPipeObject.h"

class vtkDoubleArray;

class CppPipeSlice : public CppPipeObject
{
public:
  CppPipeSlice(vtkDoubleArray* globalBounds);

  CppPipeSlice(const Json::Value& jv,
                vtkDoubleArray* globalBounds);

  virtual ~CppPipeSlice();

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
