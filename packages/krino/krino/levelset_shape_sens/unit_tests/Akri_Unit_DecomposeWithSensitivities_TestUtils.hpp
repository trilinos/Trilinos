#ifndef KRINO_KRINO_LEVELSET_SHAPE_SENS_UNIT_TESTS_AKRI_UNIT_DECOMPOSEWITHSENSITIVITIES_TESTUTILS_HPP_
#define KRINO_KRINO_LEVELSET_SHAPE_SENS_UNIT_TESTS_AKRI_UNIT_DECOMPOSEWITHSENSITIVITIES_TESTUTILS_HPP_

#include <memory>
#include <string>
#include <vector>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_LS_Field.hpp> // LS_Field
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; class MetaData; class Part; } }

// Shared helpers for the DecomposeWithSensitivities test suite. The suite is
// split by geometry family across several translation units
// (_LineOrPlane, _CircleOrSphere, _Phases); the helpers declared here are the
// ones used by more than one of those files and therefore live in a single
// shared translation unit to avoid duplicate definitions.

stk::mesh::Selector build_output_selector(const stk::mesh::MetaData & meta, const stk::mesh::Part & activePart);

void output_nonvoid_mesh(stk::mesh::BulkData & mesh, const std::string & fileName, int step, double time);

void output_full_mesh(stk::mesh::BulkData & mesh, const std::string & fileName, int step, double time);

std::unique_ptr<krino::BoundingBoxMesh> build_circle_or_sphere_conforming_bounding_box_mesh(
    const int dim,
    const stk::math::Vector3d & minCorner,
    const stk::math::Vector3d & maxCorner,
    const double meshSize,
    const stk::math::Vector3d & sphereCenter,
    const double sphereRadius,
    const bool doComputeClosestPointSensitivities,
    const bool doSnapping,
    const bool doWriteMesh,
    std::vector<krino::LS_Field> & lsFields);

#endif /* KRINO_KRINO_LEVELSET_SHAPE_SENS_UNIT_TESTS_AKRI_UNIT_DECOMPOSEWITHSENSITIVITIES_TESTUTILS_HPP_ */
