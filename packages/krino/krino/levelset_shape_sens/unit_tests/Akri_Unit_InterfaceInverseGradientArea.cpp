#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceInverseGradientArea.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SurfaceNodalNormal.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <cmath>
#include <memory>
#include <vector>

// Tests of the interface inverse gradient area on planar interfaces.  Every interface side lies in
// the plane, so its area vector is |A_s| times the plane's normal, and if the level set gradient of
// every side's negative-phase parent element has magnitude s, the value at each interface node is
// |N|/s, where N is the node's surface nodal normal (the sum of its sides' area vectors).

namespace
{

struct PlaneCase
{
  stk::math::Vector3d normal;
  double offset{0.};
  double slopeOnPositiveSide{1.};
  double slopeOnNegativeSide{1.};
  double meshSize{0.25};
  bool doSnapping{false};
  bool isOnMeshFaces{false}; // non-flattened BCC in 3D, so the plane x+y=0 lies on tet faces
};

double levelset_for_plane(const PlaneCase & plane, const stk::math::Vector3d & x)
{
  const double dist = Dot(plane.normal, x) - plane.offset;
  if (std::abs(dist) < 1.e-12) return 0.;
  return dist * ((dist > 0.) ? plane.slopeOnPositiveSide : plane.slopeOnNegativeSide);
}

// Returns the number of interface nodes checked on all processors.
size_t check_interface_inverse_gradient_area_for_plane(const unsigned dim, const PlaneCase & plane)
{
  const stk::topology elemTopo =
      (dim == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4;
  auto bboxMesh = std::make_unique<krino::BoundingBoxMesh>(elemTopo, MPI_COMM_WORLD);
  if (plane.isOnMeshFaces && dim == 3)
    bboxMesh->set_mesh_structure_type(krino::BCC_BOUNDING_BOX_MESH);
  stk::mesh::MetaData & meta = bboxMesh->meta_data();
  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::
      setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(meta, 1);
  krino::ConformingDecompositionWithSensitivities decomposition(meta, lsFields, plane.doSnapping);

  const krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);
  const stk::mesh::Selector interfaceSelector =
      phaseSupport.get_negative_levelset_interface_selector(lsFields[0].identifier);
  const krino::FieldRef normalField =
      krino::register_surface_nodal_normal_field(meta, "TEST_NORMAL", interfaceSelector);

  krino::populate_bounding_box_mesh_and_activate(
      *bboxMesh, {-1., -1., -1.}, {1., 1., 1.}, plane.meshSize);
  stk::mesh::BulkData & mesh = bboxMesh->bulk_data();

  const krino::FieldRef coordsField(meta.coordinate_field());
  auto lsData = lsFields[0].isovar.data<double, stk::mesh::ReadWrite>();
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(lsFields[0].isovar)))
    for (auto node : *bucketPtr)
      lsData.entity_values(node)() =
          levelset_for_plane(plane, get_vector_field(mesh, coordsField, node, dim));

  decomposition.decompose();

  const krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meta);
  krino::update_surface_nodal_normal(mesh,
      cdfemSupport.get_coords_field(),
      normalField,
      interfaceSelector,
      phaseSupport.get_negative_levelset_block_selector(lsFields[0].identifier) &
          krino::AuxMetaData::get(meta).active_part());
  const krino::FieldRef invGradAreaField =
      krino::update_interface_inverse_gradient_area(mesh, lsFields[0]);

  size_t numChecked = 0;
  size_t numChildNodesOnInterface = 0;
  auto invGradAreaData = invGradAreaField.data<double>();
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::NODE_RANK, interfaceSelector & meta.locally_owned_part()))
  {
    for (auto node : *bucketPtr)
    {
      if (mesh.bucket(node).member(cdfemSupport.get_child_node_part())) ++numChildNodesOnInterface;
      const double gold =
          get_vector_field(mesh, normalField, node, dim).length() / plane.slopeOnNegativeSide;
      EXPECT_NEAR(gold, invGradAreaData.entity_values(node)(), 1.e-12 * gold);
      ++numChecked;
    }
  }
  krino::all_reduce_sum(mesh.parallel(), numChecked, numChildNodesOnInterface);
  if (plane.isOnMeshFaces)
  {
    EXPECT_EQ(0u, numChildNodesOnInterface) << "The interface is not on mesh faces.";
  }
  return numChecked;
}

PlaneCase off_axis_plane(const unsigned dim)
{
  PlaneCase plane;
  plane.normal =
      (dim == 2) ? stk::math::Vector3d(0.6, 0.8, 0.) : stk::math::Vector3d(0.48, 0.6, 0.64);
  plane.offset = 0.1;
  plane.slopeOnPositiveSide = 2.;
  plane.slopeOnNegativeSide = 2.;
  return plane;
}

// 2D: y = 0 is a row of the triangular lattice for mesh size 0.2 on [-1,1] (12 rows).  3D: x+y = 0
// contains faces of the BCC tets.  The level set's slope is halved on the negative side, so each
// side's negative-phase element (gradient 0.5) gives a different value than its positive-phase one.
PlaneCase plane_on_mesh_faces_with_gradient_jump(const unsigned dim)
{
  PlaneCase plane;
  plane.normal = (dim == 2) ? stk::math::Vector3d(0., 1., 0.) :
                              stk::math::Vector3d(1., 1., 0.) / std::sqrt(2.);
  plane.offset = 0.;
  plane.slopeOnPositiveSide = 1.;
  plane.slopeOnNegativeSide = 0.5;
  plane.meshSize = (dim == 2) ? 0.2 : 1. / 3.;
  plane.isOnMeshFaces = true;
  return plane;
}

} // namespace

TEST(InterfaceInverseGradientArea, plane2d_offsetFromNodes)
{
  EXPECT_GT(check_interface_inverse_gradient_area_for_plane(2, off_axis_plane(2)), 0u);
}

TEST(InterfaceInverseGradientArea, plane3d_offsetFromNodes)
{
  EXPECT_GT(check_interface_inverse_gradient_area_for_plane(3, off_axis_plane(3)), 0u);
}

// Snapped sides stay in the plane, and the gradients use the unsnapped parent element coordinates.
TEST(InterfaceInverseGradientArea, plane2d_offsetFromNodes_snapping)
{
  PlaneCase plane = off_axis_plane(2);
  plane.doSnapping = true;
  EXPECT_GT(check_interface_inverse_gradient_area_for_plane(2, plane), 0u);
}

TEST(InterfaceInverseGradientArea, plane3d_offsetFromNodes_snapping)
{
  PlaneCase plane = off_axis_plane(3);
  plane.doSnapping = true;
  EXPECT_GT(check_interface_inverse_gradient_area_for_plane(3, plane), 0u);
}

TEST(InterfaceInverseGradientArea, plane2d_onMeshFaces_gradientJump_usesNegativeSideElement)
{
  EXPECT_GT(
      check_interface_inverse_gradient_area_for_plane(2, plane_on_mesh_faces_with_gradient_jump(2)),
      0u);
}

TEST(InterfaceInverseGradientArea, plane3d_onMeshFaces_gradientJump_usesNegativeSideElement)
{
  EXPECT_GT(
      check_interface_inverse_gradient_area_for_plane(3, plane_on_mesh_faces_with_gradient_jump(3)),
      0u);
}

TEST(InterfaceInverseGradientArea, noInterface_noInterfaceNodes)
{
  PlaneCase plane = off_axis_plane(3);
  plane.offset = 5.;
  EXPECT_EQ(0u, check_interface_inverse_gradient_area_for_plane(3, plane));
}
