// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

using LevelSetFunction = std::function<double(const stk::math::Vector3d &)>;

// Set KRINO_UNIT_TEST_OUTPUT (to anything) to write each decomposition to <test suite>_<test>.e in
// the working directory, and the mesh before it to <test suite>_<test>_background.e, to view them.
bool write_test_output() { return nullptr != std::getenv("KRINO_UNIT_TEST_OUTPUT"); }

std::string test_output_file_name(const std::string & suffix = "")
{
  const ::testing::TestInfo * testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
  return std::string(testInfo->test_suite_name()) + "_" + testInfo->name() + suffix + ".e";
}

// A bounding box mesh, a flat-walled BCC mesh in 3D or triangular lattice in 2D, with one level set
// per phase, set at the nodes from linear functions, so that the interfaces are planes and the
// phase volumes of an exact decomposition are known.
class PlanarPhasesDecomposition
{
public:
  PlanarPhasesDecomposition(const BoundingBoxMesh::BoundingBoxType & domain,
      const double meshSize,
      const stk::topology elemTopology,
      const unsigned numPhases,
      const Edge_Degeneracy_Handling edgeDegeneracyHandling,
      const std::optional<bool> cutSharpFeaturesBeforeSurfaces,
      const BoundingBoxMeshStructureType meshStructure,
      const std::optional<double> edgeTol = std::nullopt)
      : myBBoxMesh(elemTopology, MPI_COMM_WORLD)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    myLSFields = LSPerPhasePolicy::setup_levelsets_on_all_blocks(meta, numPhases);

    CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
    cdfemSupport.set_coords_field(meta.coordinate_field());
    cdfemSupport.set_cdfem_edge_degeneracy_handling(edgeDegeneracyHandling);
    if (edgeTol) cdfemSupport.set_cdfem_edge_tol(*edgeTol);
    if (cutSharpFeaturesBeforeSurfaces)
      cdfemSupport.set_cut_sharp_features_before_surfaces(*cutSharpFeaturesBeforeSurfaces);
    cdfemSupport.register_cdfem_snap_displacements_field();
    cdfemSupport.register_parent_node_ids_field();
    for (auto && lsField : myLSFields)
      cdfemSupport.add_interpolation_field(lsField.isovar);
    cdfemSupport.finalize_fields();

    myBBoxMesh.set_mesh_structure_type(meshStructure);
    populate_bounding_box_mesh_and_activate(
        myBBoxMesh, domain.get_min(), domain.get_max(), meshSize);
  }

  // Decompose for the level sets, one per phase, and measure the volume of each phase.
  std::vector<double>
  decompose_and_measure_phase_volumes(const std::vector<LevelSetFunction> & levelSets)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    STK_ThrowRequire(levelSets.size() == myLSFields.size());

    const unsigned dim = meta.spatial_dimension();
    const FieldRef coordsField(meta.coordinate_field());
    auto coordsData = coordsField.data<double>();
    for (auto && bucket : mesh.buckets(stk::topology::NODE_RANK))
    {
      for (auto node : *bucket)
      {
        const stk::math::Vector3d x(coordsData.entity_values(node), dim);
        for (size_t i = 0; i < levelSets.size(); ++i)
          myLSFields[i].isovar.data<double, stk::mesh::ReadWrite>().entity_values(node)() =
              levelSets[i](x);
      }
    }

    AuxMetaData & auxMeta = AuxMetaData::get(meta);
    const std::unique_ptr<InterfaceGeometry> geometry = create_levelset_geometry(
        dim, auxMeta.active_part(), CDFEM_Support::get(meta), Phase_Support::get(meta), myLSFields);
    if (write_test_output())
      output_composed_mesh_with_fields(
          mesh, auxMeta.active_part(), test_output_file_name("_background"), 1, 0.);
    CDMesh::decompose_mesh(mesh, *geometry);
    myNumFakeCrossingsCut = CDMesh::get_decomposed_mesh(mesh).num_fake_crossings_cut();
    if (write_test_output())
      output_composed_mesh_with_fields(mesh, auxMeta.active_part(), test_output_file_name(), 1, 0.);

    const stk::mesh::Selector ownedActive = meta.locally_owned_part() & auxMeta.active_part();
    std::vector<double> volumes;
    for (size_t phase = 0; phase < levelSets.size(); ++phase)
    {
      const stk::mesh::Selector phaseSelector =
          auxMeta.get_part("block_1_P" + std::to_string(phase)) & ownedActive;
      double volume =
          compute_volume_or_surface_area(mesh, stk::topology::ELEMENT_RANK, phaseSelector);
      all_reduce_sum(mesh.parallel(), volume);
      volumes.push_back(volume);
    }

    myInterfaceAreas.assign(levelSets.size(), std::vector<double>(levelSets.size(), 0.));
    const stk::mesh::Selector ownedSides = meta.locally_owned_part();
    for (size_t phase0 = 0; phase0 < levelSets.size(); ++phase0)
    {
      for (size_t phase1 = phase0 + 1; phase1 < levelSets.size(); ++phase1)
      {
        const stk::mesh::Part * interfacePart = Phase_Support::get(meta).find_interface_part(
            phase_part(meta, phase0), phase_part(meta, phase1));
        double area = (nullptr == interfacePart) ?
            0. :
            compute_volume_or_surface_area(mesh, meta.side_rank(), *interfacePart & ownedSides);
        all_reduce_sum(mesh.parallel(), area);
        myInterfaceAreas[phase0][phase1] = myInterfaceAreas[phase1][phase0] = area;
      }
    }
    return volumes;
  }

  // The leaf edges cut where an interface's plane crosses them but another phase is smallest, in
  // the last decomposition (CDMesh::num_fake_crossings_cut).
  size_t num_fake_crossings_cut() const { return myNumFakeCrossingsCut; }

  // The area (length in 2D) of the interface between two phases, of the last decomposition.
  double interface_area(const size_t phase0, const size_t phase1) const
  {
    return myInterfaceAreas[phase0][phase1];
  }

private:
  static const stk::mesh::Part & phase_part(const stk::mesh::MetaData & meta, const size_t phase)
  {
    return AuxMetaData::get(meta).get_part("block_1_P" + std::to_string(phase));
  }

  BoundingBoxMesh myBBoxMesh;
  std::vector<LS_Field> myLSFields;
  std::vector<std::vector<double>> myInterfaceAreas;
  size_t myNumFakeCrossingsCut = 0;
};

// An octagon of apothem 0.3, centered in the unit square [-0.5,0.5]^2 (or the unit cube, as a prism
// along z), cut from it by 8 planes, each n.x + 0.3 with unit normal n: the phase of each plane is
// beyond it, and the middle phase, with level set 0, is the octagon.  The planes x=+-0.3 and
// y=+-0.3 go through nodes of the meshes, so that faces of leaves lie in them.  Each phase volume
// is the area of its region of the square (times the unit depth in 3D).
std::vector<LevelSetFunction> octagon_level_sets()
{
  const std::array<std::array<double, 2>, 8> normals{{{{0., 1.}},
      {{1., 1.}},
      {{1., 0.}},
      {{1., -1.}},
      {{0., -1.}},
      {{-1., -1.}},
      {{-1., 0.}},
      {{-1., 1.}}}};
  std::vector<LevelSetFunction> levelSets;
  for (auto && normal : normals)
  {
    const double length = std::hypot(normal[0], normal[1]);
    const double nx = normal[0] / length;
    const double ny = normal[1] / length;
    levelSets.push_back(
        [nx, ny](const stk::math::Vector3d & x) { return nx * x[0] + ny * x[1] + 0.3; });
  }
  levelSets.push_back([](const stk::math::Vector3d &) { return 0.; });
  return levelSets;
}

const std::array<double, 9> octagonPhaseAreas{0.06627416997969524,
    0.10916738879314766,
    0.06627416997969524,
    0.10916738879314766,
    0.06627416997969524,
    0.10916738879314766,
    0.06627416997969524,
    0.10916738879314766,
    0.2982337649086284};

// The area (length in 2D) of the interface between two phases of the octagon, per unit depth: a
// side of the octagon, between the middle phase and an outer one; the ray from a corner of the
// octagon to the wall of the square, between neighboring outer phases (all at odd multiples of
// 22.5 degrees from the axes); none between other phases.
double octagon_interface_area(const size_t phase0, const size_t phase1)
{
  constexpr size_t middle = 8;
  const double halfAngle = std::acos(-1.) / 8.;
  if (phase0 == middle || phase1 == middle) return 2. * 0.3 * std::tan(halfAngle);
  const size_t gap = (phase0 > phase1) ? phase0 - phase1 : phase1 - phase0;
  if (gap == 1 || gap == 7) return (0.5 - 0.3) / std::cos(halfAngle);
  return 0.;
}

constexpr Edge_Degeneracy_Handling onlyCutting = SNAP_TO_NODE;
constexpr Edge_Degeneracy_Handling snapping =
    SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE;

void expect_octagon_phase_volumes(const stk::topology elemTopology,
    const Edge_Degeneracy_Handling edgeDegeneracyHandling,
    const std::optional<bool> cutSharpFeaturesBeforeSurfaces,
    const BoundingBoxMeshStructureType meshStructure,
    const double meshSize,
    const double tol,
    const std::optional<double> edgeTol = std::nullopt)
{
  const bool is3D = elemTopology == stk::topology::TETRAHEDRON_4;
  const double halfDepth = is3D ? 0.5 : 0.;
  const BoundingBoxMesh::BoundingBoxType domain(
      stk::math::Vector3d(-0.5, -0.5, -halfDepth), stk::math::Vector3d(0.5, 0.5, halfDepth));
  PlanarPhasesDecomposition decomposition(domain,
      meshSize,
      elemTopology,
      9,
      edgeDegeneracyHandling,
      cutSharpFeaturesBeforeSurfaces,
      meshStructure,
      edgeTol);
  const std::vector<double> volumes =
      decomposition.decompose_and_measure_phase_volumes(octagon_level_sets());
  for (size_t phase = 0; phase < volumes.size(); ++phase)
  {
    EXPECT_NEAR(octagonPhaseAreas[phase], volumes[phase], tol) << "phase " << phase;
  }
  for (size_t phase0 = 0; phase0 < volumes.size(); ++phase0)
  {
    for (size_t phase1 = phase0 + 1; phase1 < volumes.size(); ++phase1)
    {
      EXPECT_NEAR(
          octagon_interface_area(phase0, phase1), decomposition.interface_area(phase0, phase1), tol)
          << "interface between phases " << phase0 << " and " << phase1;
    }
  }
  // With the sharp features cut first, no plane cuts a leaf where it is no interface.
  if (cutSharpFeaturesBeforeSurfaces.value_or(true))
  {
    EXPECT_EQ(0u, decomposition.num_fake_crossings_cut());
  }
}

void expect_octagon_phase_volumes(const stk::topology elemTopology,
    const Edge_Degeneracy_Handling edgeDegeneracyHandling,
    const std::optional<bool> cutSharpFeaturesBeforeSurfaces)
{
  const bool is3D = elemTopology == stk::topology::TETRAHEDRON_4;
  expect_octagon_phase_volumes(elemTopology,
      edgeDegeneracyHandling,
      cutSharpFeaturesBeforeSurfaces,
      is3D ? FLAT_WALLED_BCC_BOUNDING_BOX_MESH : FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH,
      is3D ? 0.2 : 0.1,
      1.e-10);
}


} // namespace

// The phases of the octagon are exact whether or not the sharp features (where three phases meet)
// are cut first, and whether only cutting or snapping.  With the sharp features cut first, the
// leaves are cut only by the interfaces that cross them; without, every plane crossing a leaf
// cuts it, also where it is no interface.

TEST(OctagonTet4Decomposition, onlyCuttingWithoutSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4, onlyCutting, false);
}

// Without the sharp features cut first, every plane crossing a leaf cuts it, fake crossings too:
// the count of them sees them.
TEST(OctagonTet4Decomposition, onlyCuttingWithoutSharpFeaturesFirst_cutsFakeCrossings)
{
  PlanarPhasesDecomposition decomposition(
      BoundingBoxMesh::BoundingBoxType(
          stk::math::Vector3d(-0.5, -0.5, -0.5), stk::math::Vector3d(0.5, 0.5, 0.5)),
      0.2,
      stk::topology::TETRAHEDRON_4,
      9,
      onlyCutting,
      false,
      FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  decomposition.decompose_and_measure_phase_volumes(octagon_level_sets());
  EXPECT_LT(0u, decomposition.num_fake_crossings_cut());
}

TEST(OctagonTet4Decomposition, onlyCuttingWithSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4, onlyCutting, true);
}

TEST(OctagonTet4Decomposition, snapping_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4, snapping, std::nullopt);
}

TEST(OctagonTet4Decomposition, snappingWithoutSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4, snapping, false);
}

TEST(OctagonTri3Decomposition, onlyCuttingWithoutSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TRIANGLE_3_2D, onlyCutting, false);
}

TEST(OctagonTri3Decomposition, onlyCuttingWithSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TRIANGLE_3_2D, onlyCutting, true);
}

TEST(OctagonTri3Decomposition, snapping_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TRIANGLE_3_2D, snapping, std::nullopt);
}

TEST(OctagonTri3Decomposition, snappingWithoutSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TRIANGLE_3_2D, snapping, false);
}

// On cubic meshes of tets, where the planes x=+-0.3 and y=+-0.3 cut through the cells.  With the
// sharp features cut first, a sharp feature is moved within the edge tolerance onto an entity it is
// near, so that the volumes are exact only to a small fraction of that.  Without the sharp features
// cut first, setting the signs of the leaves from their crossings, as with them, would miss where
// the interfaces cross leaves beyond their uncaptured sharp features (by 6e-3 in a phase volume).
constexpr double cubicMeshTol = 1.e-6;

TEST(OctagonTet4CubicMeshDecomposition, snappingWithoutSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4,
      snapping,
      false,
      CUBIC_BOUNDING_BOX_MESH,
      1. / 3.,
      cubicMeshTol);
}

TEST(OctagonTet4CubicMeshDecomposition, snapping_phaseVolumesNearlyExact)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4,
      snapping,
      std::nullopt,
      CUBIC_BOUNDING_BOX_MESH,
      0.25,
      cubicMeshTol);
}

TEST(OctagonTet4CubicMeshDecomposition, onlyCuttingWithoutSharpFeaturesFirst_exactPhaseVolumes)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4,
      onlyCutting,
      false,
      CUBIC_BOUNDING_BOX_MESH,
      0.25,
      cubicMeshTol);
}

// With the edge tolerance that snapping uses in an input file (0.05).  A sharp-edge crossing found
// again beside a node that captures it, as the sharp edge through the node crosses a leaf face
// beside it, is that node: cut, it was moved off the sharp edge by the tolerance (by 3e-5 in a
// phase volume here).
TEST(OctagonTet4CubicMeshDecomposition, snappingWithSnapEdgeTolerance_phaseVolumesNearlyExact)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4,
      snapping,
      std::nullopt,
      CUBIC_BOUNDING_BOX_MESH,
      0.25,
      cubicMeshTol,
      0.05);
}

// On a coarser mesh, a cut of a sharp-edge crossing splits leaves into new faces, one of which
// another sharp edge crosses inside it: that crossing is cut in a later round (left uncut, a leaf
// with three phases took the volume of one, by 2e-3).
TEST(OctagonTet4CubicMeshDecomposition,
    snappingWithSnapEdgeToleranceOnCoarserMesh_phaseVolumesNearlyExact)
{
  expect_octagon_phase_volumes(stk::topology::TETRAHEDRON_4,
      snapping,
      std::nullopt,
      CUBIC_BOUNDING_BOX_MESH,
      1. / 3.,
      cubicMeshTol,
      0.05);
}

TEST(OctagonTet4CubicMeshDecomposition, onlyCuttingWithSharpFeaturesFirst_phaseVolumesNearlyExact)
{
  expect_octagon_phase_volumes(
      stk::topology::TETRAHEDRON_4, onlyCutting, true, CUBIC_BOUNDING_BOX_MESH, 0.25, cubicMeshTol);
}

} // namespace krino
