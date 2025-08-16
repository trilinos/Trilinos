#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNITMESHUTILS_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNITMESHUTILS_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>
#include <string>

namespace krino {

class FieldRef;
class BoundingBoxMesh;

stk::mesh::Entity find_local_node_closest_to_location(const stk::mesh::BulkData& mesh, const stk::mesh::Selector& nodeSelector, const FieldRef coordsField, const stk::math::Vector3d& location);

void populate_bounding_box_mesh_and_activate(BoundingBoxMesh & bboxMesh, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize);
void generate_and_write_bounding_box_mesh(const stk::topology elemTopology, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize, const std::string & filename);

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNITMESHUTILS_HPP_ */
