#ifndef CPPPIPECAMERA_H
#define CPPPIPECAMERA_H

#include "CppPipeObject.h"

class vtkDoubleArray;

class CppPipeCamera : public CppPipeObject
{
public:
  CppPipeCamera(vtkDoubleArray* globalBounds);

  CppPipeCamera(const Json::Value& jv,
                 vtkDoubleArray* globalBounds);

  virtual ~CppPipeCamera();

  void InitializeProxy(vtkSMProxy* p);

private:
  void SetGlobalBounds(vtkDoubleArray* globalBounds);
  void SetLookPoint(double* lookPoint);
  void SetDefaults();
  vtkDoubleArray* gb;
  double inFocalPoint[3];
  double intX;
  double intY;
  double intZ;
};
#endif
