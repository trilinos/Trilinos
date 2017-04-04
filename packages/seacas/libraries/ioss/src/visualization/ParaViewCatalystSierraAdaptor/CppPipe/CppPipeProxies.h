#ifndef CPPPIPEPROXIES_H
#define CPPPIPEPROXIES_H

#include "CppPipeObject.h"
#include "CppPipeImageSet.h"
#include "CppPipeMultiCamera8.h"
#include "CppPipeCamera.h"
#include "CppPipeRepresentation.h"

class vtkCPDataDescription;
class vtkSMSourceProxy;
class vtkSMRenderViewProxy;
class vtkCPDataDescription;
class vtkDoubleArray;
class vtkSMRepresentationProxy;

class CppPipeProxies
{
public:

  CppPipeProxies(vtkCPDataDescription* dataDescription);

  virtual ~CppPipeProxies();

  void RenderImageSet(const std::string& ims_name,
                      const std::string& image_basename,
                      const Json::Value& json);

private:
  void GetGlobalDataBoundsParallel(vtkDoubleArray* globalBounds);
  void RenderMultiCamera8(CppPipeImageSet& ims,
                          CppPipeMultiCamera8& mc,
                          CppPipeRepresentation& phrep,
                          const std::string& image_basename,
                          vtkSMRepresentationProxy* representation);
  void RenderCamera(CppPipeImageSet& ims,
                    CppPipeCamera& c,
                    CppPipeRepresentation& phrep,
                    const std::string& image_basename,
                    vtkSMRepresentationProxy* representation);
  CppPipeProxies();

  vtkSMRepresentationProxy* SetupOperations(const std::string& ims_name, 
                                            const Json::Value& json);

  vtkSMSourceProxy* SetupSlice(const Json::Value& json,
                               vtkSMSourceProxy* p);

  vtkSMSourceProxy* SetupClip(const Json::Value& json,
                              vtkSMSourceProxy* p);

  vtkSMSourceProxy* SetupBoxClip(const Json::Value& json,
                                 vtkSMSourceProxy* p);

  vtkSMSourceProxy* SetupThreshold(const Json::Value& json,
                                   vtkSMSourceProxy* p);

  vtkSMSourceProxy* SetupContour(const Json::Value& json,
                                 vtkSMSourceProxy* p);

  vtkSMSourceProxy* producer;
  vtkSMSourceProxy* slice;
  vtkSMSourceProxy* clip;
  vtkSMSourceProxy* threshold;
  vtkSMSourceProxy* contour;
  vtkSMProxy* clip_plane;
  vtkSMProxy* slice_plane;
  vtkSMProxy* box;
  vtkSMRenderViewProxy* render;
  vtkCPDataDescription* dataDescription;
  vtkDoubleArray* globalBounds;
};
#endif
