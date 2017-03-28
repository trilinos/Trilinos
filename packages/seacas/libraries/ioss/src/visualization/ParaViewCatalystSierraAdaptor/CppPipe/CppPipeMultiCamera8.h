#ifndef CPPPIPEMULTICAMERA8_H
#define CPPPIPEMULTICAMERA8_H

#include "CppPipeObject.h"

class vtkDoubleArray;

class CppPipeMultiCamera8 : public CppPipeObject
{
public:
  CppPipeMultiCamera8(vtkDoubleArray* globalBounds);

  CppPipeMultiCamera8(const Json::Value& jv, 
                       vtkDoubleArray* globalBounds);

  virtual ~CppPipeMultiCamera8();

  void InitializeProxy(vtkSMProxy* p);

  enum Direction {X1, X2, Y1, Y2, 
                  Z1, Z2, XYZ1, XYZ2};

  static const char* DirectionStrings[];

  void SetCameraDirection(Direction d);

private:
  void SetGlobalBounds(vtkDoubleArray* globalBounds);
  void SetDefaults();
  void CheckParallelVectors(double* look_direction, 
                            double* view_up);
  Direction d;
  vtkDoubleArray* gb;
  double maxDim;
  double intX;
  double intY;
  double intZ;
  double inFocalPoint[3];
};
#endif
