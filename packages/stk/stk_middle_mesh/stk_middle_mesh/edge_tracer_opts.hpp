#ifndef EDGE_TRACER_OPTS_H
#define EDGE_TRACER_OPTS_H

#include <algorithm>

namespace stk {
namespace middle_mesh {
namespace impl {

struct EdgeTracerTolerances
{
    explicit EdgeTracerTolerances(double eps = 1e-12, double edgeErrorTol = 2)
      : inRangeTol(std::max(eps, 1e-9))
      , endpointRedistributionDistance(eps)
      , vertexClassificationTol(eps)
      , edgeIntersectionPrimitiveZeroTol(eps)
      , edgeIntersectionPrimitiveErrorTol(edgeErrorTol)
      , rangeNarrowingTol(0.1)
    {}

    double inRangeTol; // tolerance for determining if alpha and beta are in the range [0, 1]

    double endpointRedistributionDistance; // if an alpha values are within this distance of 0 or 1,
                                           // redistribute them so they are evenly spaced within this distance.
                                           // This avoids problems when an alpha that is slightly greater than 1 or
                                           // less than zero is chosen as the "best" intersection

    double vertexClassificationTol; // If beta is within this distance of 0 or 1, classify it is a vertex
                                    // intersection

    double edgeIntersectionPrimitiveZeroTol; // tolerance for edge intersection primitive to consider values to be zero

    double edgeIntersectionPrimitiveErrorTol; // tolerance for edge intersection primitive to reject solutions if the
                                              // angle between the direction vector and the displacement vector is
                                              // too large.  Value is in radians

    double rangeNarrowingTol; // If no intersections are in the range [0, 1] (to within m_inRangeTol), use this larger
                              // tolerance to try to find an intersection
};

} // namespace impl
} // namespace middle_mesh
} // namespace stk

#endif