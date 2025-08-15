#ifndef KRINO_KRINO_GEOMETRY_AKRI_SEGMENTWITHSENSITIVITIES_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_SEGMENTWITHSENSITIVITIES_HPP_

#include <stk_math/StkVector.hpp>

namespace krino {
namespace SegmentWithSens {

double length_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dLength);

stk::math::Vector3d normal2d_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dNormal);

} // namespace SegmentWithSens
} // namespace krino



#endif /* KRINO_KRINO_GEOMETRY_AKRI_SEGMENTWITHSENSITIVITIES_HPP_ */
