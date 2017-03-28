#ifndef CPPPIPEREPRESENTATION_H
#define CPPPIPEREPRESENTATION_H

#include <vtkDataObject.h>
#include "CppPipeObject.h"

class vtkCPDataDescription;
class vtkCompositeDataSet;
class vtkDoubleArray;

class CppPipeRepresentation : public CppPipeObject
{
public:
  CppPipeRepresentation(vtkCPDataDescription* dataDescription,
                         vtkSMProxy* view);

  CppPipeRepresentation(const Json::Value& jv,
                         vtkCPDataDescription* dataDescription,
                         vtkSMProxy* view);

  virtual ~CppPipeRepresentation() {};

  void InitializeProxy(vtkSMProxy* p);

private:
  vtkDataObject::AttributeTypes GetAttributeType(vtkCompositeDataSet* cds,
                                                 const char* array_name,
                                                 vtkDoubleArray* localArrayBounds);

  void GetGlobalArrayBoundsParallel(vtkDoubleArray* globalBounds,
                                    vtkDoubleArray* localBounds);

  void InitializeLookupTable(vtkSMProxy* lut,
                             double min,
                             double max);

  void SetDefaults();
  Json::Value colorLegendCollection;
  vtkCPDataDescription* dataDescription;
  vtkSMProxy* view;
};
#endif
