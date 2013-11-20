#ifndef UNITTESTUTILS_HPP
#define UNITTESTUTILS_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

const int spatialDim = 3;
typedef stk::search::ident::IdentProc<unsigned, unsigned> MyBoxId;
typedef stk::search::box::SphereBoundingBox<MyBoxId, double, spatialDim> My3dSphereBoundingBox;
typedef stk::search::box::PointBoundingBox<MyBoxId, double, spatialDim> My3dPointBoundingBox;
typedef stk::search::box::AxisAlignedBoundingBox<MyBoxId, double, spatialDim> My3dAxisAlignedBoundingBox;

template<class BoundingBoxType>
BoundingBoxType generate3dBoundingBox(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId);

template<>
inline My3dPointBoundingBox generate3dBoundingBox<My3dPointBoundingBox>(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId)
{
    MyBoxId id(entityId, procId);
    double center[spatialDim];
    center[0] = centerX;
    center[1] = centerY;
    center[2] = centerZ;
    return My3dPointBoundingBox(center, id);
}

template<>
inline My3dSphereBoundingBox generate3dBoundingBox<My3dSphereBoundingBox>(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId)
{
    MyBoxId id(entityId, procId);
    double center[spatialDim];
    center[0] = centerX;
    center[1] = centerY;
    center[2] = centerZ;
    return My3dSphereBoundingBox(center, radius, id);
}

//       ------------
//      |            |
//      |      radius|
//      |      ------|
//      |            |
//      |            |
//       ------------
// width = 2*radius
template<>
inline My3dAxisAlignedBoundingBox generate3dBoundingBox<My3dAxisAlignedBoundingBox>(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId)
{
    MyBoxId id(entityId, procId);
    double corners[2*spatialDim];
    corners[0] = centerX-radius;
    corners[1] = centerY-radius;
    corners[2] = centerZ-radius;
    corners[3] = centerX+radius;
    corners[4] = centerY+radius;
    corners[5] = centerZ+radius;
    return My3dAxisAlignedBoundingBox(corners, id);
}

#endif
