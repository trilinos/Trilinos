#ifndef CPPPIPEIMAGESET_H
#define CPPPIPEIMAGESET_H

#include "CppPipeObject.h"

class CppPipeImageSet : public CppPipeObject
{
public:
  CppPipeImageSet();

  CppPipeImageSet(const Json::Value& jv);

  virtual ~CppPipeImageSet() {};

  void InitializeProxy(vtkSMProxy* proxy);

  std::string GetBaseDirectoryName() const; 

  std::string GetImageFormat() const;

  std::string GetImageBaseName() const;

private:
  void SetDefaults();
};
#endif
