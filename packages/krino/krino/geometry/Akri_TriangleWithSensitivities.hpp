#ifndef KRINO_KRINO_GEOMETRY_AKRI_TRIANGLEWITHSENSITIVITIES_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_TRIANGLEWITHSENSITIVITIES_HPP_
#include <stk_math/StkVector.hpp>

namespace krino {
namespace TriangleWithSens {

double area_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dArea);

stk::math::Vector3d normal_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dNormal);

stk::math::Vector3d area_vector_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dAreaVector);

} // namespace TriangleWithSens
} // namespace krino


#endif /* KRINO_KRINO_GEOMETRY_AKRI_TRIANGLEWITHSENSITIVITIES_HPP_ */
