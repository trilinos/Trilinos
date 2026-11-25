// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <array>
#include <memory>
#include <optional>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshBlockSkinSurface.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_ScopedTestFile.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_Unit_MatingFacetedSurfaces.hpp>
#include <Akri_Unit_MultiBlockMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace krino
{

namespace
{

// The volume of each surface's phase and the area of its interface (its whole boundary, with the other
// phases and the elsewhere phase), the volume of the elsewhere phase, and the area of all the interfaces.
struct MatingDecompositionMeasures
{
  std::vector<double> phaseVolumes;
  std::vector<double> phaseInterfaceAreas;
  double elsewhereVolume;
  double allInterfacesArea;
  double worstQuality; // mean ratio
};

// A flat-walled BCC tet4 mesh of a box (or a flat-walled triangular lattice of tri3 in 2D) with
// faceted bounding surfaces as a region sets them up: each surface's phase on its negative side, and
// the elsewhere phase on the positive side of all of them.
class MatingSurfacesDecomposition
{
public:
  MatingSurfacesDecomposition(const BoundingBoxMesh::BoundingBoxType & domain,
      const double meshSize,
      const unsigned numSurfaces,
      const Edge_Degeneracy_Handling edgeDegeneracyHandling,
      const stk::topology elemTopology = stk::topology::TETRAHEDRON_4)
      : myBBoxMesh(elemTopology, MPI_COMM_WORLD), myNumSurfaces(numSurfaces)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    myLSFields = LSPerBoundingSurfacePolicy::setup_levelsets_on_all_blocks(meta, numSurfaces);

    CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
    cdfemSupport.set_coords_field(meta.coordinate_field());
    cdfemSupport.set_cdfem_edge_degeneracy_handling(edgeDegeneracyHandling);
    cdfemSupport.register_cdfem_snap_displacements_field();
    cdfemSupport.register_parent_node_ids_field();
    cdfemSupport.finalize_fields();

    myBBoxMesh.set_mesh_structure_type(elemTopology == stk::topology::TETRAHEDRON_4 ?
            FLAT_WALLED_BCC_BOUNDING_BOX_MESH :
            FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH);
    populate_bounding_box_mesh_and_activate(myBBoxMesh, domain.get_min(), domain.get_max(), meshSize);
  }

  MatingDecompositionMeasures decompose_and_measure(const std::vector<const Surface *> & surfaces)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    const AuxMetaData & auxMeta = AuxMetaData::get(meta);
    const Phase_Support & phaseSupport = Phase_Support::get(meta);

    const std::vector<Surface_Identifier> surfaceIdentifiers = surface_identifiers();
    const AnalyticSurfaceInterfaceGeometry geometry(
        surfaceIdentifiers, surfaces, auxMeta.active_part(), CDFEM_Support::get(meta), phaseSupport);
    CDMesh::decompose_mesh(mesh, geometry);

    const stk::mesh::Selector ownedActive = meta.locally_owned_part() & auxMeta.active_part();
    MatingDecompositionMeasures measures;
    for (const Surface_Identifier surfaceIdentifier : surfaceIdentifiers)
    {
      measures.phaseVolumes.push_back(compute_volume_or_surface_area(mesh,
          stk::topology::ELEMENT_RANK,
          phaseSupport.get_negative_levelset_block_selector(surfaceIdentifier) & ownedActive));
      measures.phaseInterfaceAreas.push_back(compute_volume_or_surface_area(mesh,
          meta.side_rank(),
          phaseSupport.get_negative_levelset_interface_selector(surfaceIdentifier) & ownedActive));
    }
    for (std::vector<double> * values : {&measures.phaseVolumes, &measures.phaseInterfaceAreas})
    {
      const std::vector<double> local = *values;
      stk::all_reduce_sum(mesh.parallel(), local.data(), values->data(), local.size());
    }
    measures.elsewhereVolume = compute_volume_or_surface_area(mesh,
        stk::topology::ELEMENT_RANK,
        phaseSupport.get_selector_of_blocks_with_phase(LSPerBoundingSurfacePolicy::elsewhere_phase_tag(myNumSurfaces)) &
            ownedActive);
    all_reduce_sum(mesh.parallel(), measures.elsewhereVolume);
    measures.allInterfacesArea = compute_volume_or_surface_area(
        mesh, meta.side_rank(), phaseSupport.get_all_interface_surfaces_selector() & ownedActive);
    all_reduce_sum(mesh.parallel(), measures.allInterfacesArea);
    measures.worstQuality =
        compute_mesh_quality(mesh, auxMeta.active_part(), meta.coordinate_field(), MeanRatioQualityMetric());
    return measures;
  }

  // Whether CDFEM would decompose the mesh again for the surfaces, as an application asks before each
  // step.
  bool does_decomposition_need_update(const std::vector<const Surface *> & surfaces)
  {
    const stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    const AnalyticSurfaceInterfaceGeometry geometry(surface_identifiers(),
        surfaces,
        AuxMetaData::get(meta).active_part(),
        CDFEM_Support::get(meta),
        Phase_Support::get(meta));
    return CDMesh::decomposition_needs_update(myBBoxMesh.bulk_data(), geometry);
  }

private:
  std::vector<Surface_Identifier> surface_identifiers() const
  {
    std::vector<Surface_Identifier> surfaceIdentifiers;
    for (unsigned i = 0; i < myNumSurfaces; ++i)
      surfaceIdentifiers.push_back(myLSFields[i].identifier);
    return surfaceIdentifiers;
  }

  BoundingBoxMesh myBBoxMesh;
  unsigned myNumSurfaces;
  std::vector<LS_Field> myLSFields;
};

const BoundingBoxMesh::BoundingBoxType matingDomain(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
const stk::math::Vector3d matingCuboidDimensions(1.0, 0.7, 0.8);
const Cuboid matingCuboid(stk::math::Vector3d(0.03, -0.02, 0.05), matingCuboidDimensions, stk::math::Vector3d(0.3, 0.5, 0.2));
const BoundingBoxMesh::BoundingBoxType matingSquareDomain(stk::math::Vector3d(-1., -1., 0.), stk::math::Vector3d(1., 1., 0.));
const stk::math::Vector3d matingRectangleCenter(0.03, -0.02, 0.);
const std::array<double, 2> matingRectangleDimensions{{1.0, 0.7}};
constexpr double matingRectangleAngle = 0.3;

// The volume of each half, the area of its interface (its whole boundary), the volume of the elsewhere
// phase, and the area of all the interfaces (with the middle face between the halves once), within relTol
// of the exact ones.
void expect_halves_of_cuboid(const MatingDecompositionMeasures & measures, const double relTol)
{
  const stk::math::Vector3d & d = matingCuboidDimensions;
  const double halfVolume = 0.5 * d[0] * d[1] * d[2];
  const double middleFaceArea = d[0] * d[1];
  const double halfAreaWithElsewhere = d[0] * d[1] + d[0] * d[2] + d[1] * d[2];
  const double halfArea = halfAreaWithElsewhere + middleFaceArea;
  const double domainVolume = 8.;
  for (unsigned i = 0; i < 2; ++i)
  {
    EXPECT_NEAR(halfVolume, measures.phaseVolumes[i], relTol * halfVolume) << "half " << i;
    EXPECT_NEAR(halfArea, measures.phaseInterfaceAreas[i], relTol * halfArea) << "half " << i;
  }
  EXPECT_NEAR(domainVolume - 2. * halfVolume, measures.elsewhereVolume, relTol * domainVolume);
  const double allInterfacesArea = 2. * halfAreaWithElsewhere + middleFaceArea;
  EXPECT_NEAR(allInterfacesArea, measures.allInterfacesArea, relTol * allInterfacesArea);
  EXPECT_LT(0., measures.worstQuality);
}

// As expect_halves_of_cuboid, for the halves of the rectangle in 2D (areas and lengths).
void expect_halves_of_rectangle(const MatingDecompositionMeasures & measures, const double relTol)
{
  const std::array<double, 2> & d = matingRectangleDimensions;
  const double halfArea = 0.5 * d[0] * d[1];
  const double middleSegmentLength = d[0];
  const double halfLengthWithElsewhere = d[0] + d[1];
  const double halfLength = halfLengthWithElsewhere + middleSegmentLength;
  const double domainArea = 4.;
  for (unsigned i = 0; i < 2; ++i)
  {
    EXPECT_NEAR(halfArea, measures.phaseVolumes[i], relTol * halfArea) << "half " << i;
    EXPECT_NEAR(halfLength, measures.phaseInterfaceAreas[i], relTol * halfLength) << "half " << i;
  }
  EXPECT_NEAR(domainArea - 2. * halfArea, measures.elsewhereVolume, relTol * domainArea);
  const double allInterfacesLength = 2. * halfLengthWithElsewhere + middleSegmentLength;
  EXPECT_NEAR(allInterfacesLength, measures.allInterfacesArea, relTol * allInterfacesLength);
  EXPECT_LT(0., measures.worstQuality);
}

// The halves of the cuboid (or of the rectangle in 2D), with sharp features at 135 degrees.
std::array<std::unique_ptr<MeshFacetedSurface>, 2> mating_halves_with_sharp_features(const bool is2D, const bool matchingMiddle)
{
  auto halves = is2D ?
      mating_halves_of_rectangle(matingRectangleCenter, matingRectangleDimensions, matingRectangleAngle, matchingMiddle) :
      mating_halves_of_cuboid(matingCuboid, matchingMiddle);
  for (auto & half : halves)
    half->set_sharp_feature_angle(135.);
  return halves;
}

MatingDecompositionMeasures decompose_halves_of_cuboid(const Edge_Degeneracy_Handling edgeDegeneracyHandling, const bool matchingMiddleFaces)
{
  const auto halves = mating_halves_with_sharp_features(false, matchingMiddleFaces);
  MatingSurfacesDecomposition decomposition(matingDomain, 0.2, 2, edgeDegeneracyHandling);
  return decomposition.decompose_and_measure({halves[0].get(), halves[1].get()});
}

MatingDecompositionMeasures decompose_halves_of_rectangle(const Edge_Degeneracy_Handling edgeDegeneracyHandling, const bool matchingMiddleSegments)
{
  const auto halves = mating_halves_with_sharp_features(true, matchingMiddleSegments);
  MatingSurfacesDecomposition decomposition(matingSquareDomain, 0.2, 2, edgeDegeneracyHandling, stk::topology::TRIANGLE_3_2D);
  return decomposition.decompose_and_measure({halves[0].get(), halves[1].get()});
}

// The halves of the cuboid as two hexes, block_1 the bottom half and block_2 the top, sharing the middle face.
MultiBlockMesh hex_halves_of_cuboid(const Cuboid & cuboid, const stk::ParallelMachine comm)
{
  std::vector<stk::math::Vector3d> nodeLocs;
  for (unsigned i = 0; i < 4; ++i)
    nodeLocs.push_back(cuboid.vertex_location(i));
  for (unsigned i = 0; i < 4; ++i)
    nodeLocs.push_back(0.5 * (cuboid.vertex_location(i) + cuboid.vertex_location(4 + i)));
  for (unsigned i = 0; i < 4; ++i)
    nodeLocs.push_back(cuboid.vertex_location(4 + i));
  return MultiBlockMesh(stk::topology::HEXAHEDRON_8, nodeLocs, {{1, 2, 3, 4, 5, 6, 7, 8}, {5, 6, 7, 8, 9, 10, 11, 12}}, comm);
}

// The halves of the cuboid as the skins of the two hex blocks of a mesh file, with sharp features at 135
// degrees.
MatingDecompositionMeasures decompose_block_skin_halves_of_cuboid(const Edge_Degeneracy_Handling edgeDegeneracyHandling)
{
  const ScopedTestFile file("MatingHalvesOfRotatedCuboidBlocks.e");
  write_mesh_file_on_processor_0(file.name(),
      [](const stk::ParallelMachine comm) { return hex_halves_of_cuboid(matingCuboid, comm); });
  const auto meshFile = get_block_skin_mesh_file(MPI_COMM_WORLD, file.name());
  MeshBlockSkinSurface bottom(3, meshFile, {"block_1"}, 1);
  MeshBlockSkinSurface top(3, meshFile, {"block_2"}, 1);
  for (auto * half : {&bottom, &top})
    half->set_sharp_feature_angle(135.);
  MatingSurfacesDecomposition decomposition(matingDomain, 0.2, 2, edgeDegeneracyHandling);
  return decomposition.decompose_and_measure({&bottom, &top});
}

// After decomposing the mesh for the halves, whether it would be decomposed again for them unchanged.
bool does_decomposition_of_halves_need_update(const bool is2D, const Edge_Degeneracy_Handling edgeDegeneracyHandling)
{
  const auto halves = mating_halves_with_sharp_features(is2D, true);
  const std::vector<const Surface *> surfaces{halves[0].get(), halves[1].get()};
  MatingSurfacesDecomposition decomposition(is2D ? matingSquareDomain : matingDomain,
      0.2,
      2,
      edgeDegeneracyHandling,
      is2D ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4);
  decomposition.decompose_and_measure(surfaces);
  return decomposition.does_decomposition_need_update(surfaces);
}

} // namespace

// Each half of the cuboid with sharp features at 135 degrees, which include the edges of the middle
// face where the halves and the elsewhere phase meet.  Where the halves mate, the crossings of the two
// surfaces are one point on both, so a node snapped there is on both: each half's volume and areas are
// exact, as for the whole cuboid.
TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedCuboid_snapping_exactVolumesAndAreas)
{
  expect_halves_of_cuboid(decompose_halves_of_cuboid(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, true), 1.e-10);
}

// The same, with the halves read as the skins of two hex blocks of a mesh file (mesh_block_skin), which
// share the middle face.
TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedCuboidAsBlockSkinsOfMeshFile_snapping_exactVolumesAndAreas)
{
  expect_halves_of_cuboid(decompose_block_skin_halves_of_cuboid(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE), 1.e-10);
}

TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedCuboidAsBlockSkinsOfMeshFile_cuttingOnly_volumesAndAreasAsForOneSurface)
{
  expect_halves_of_cuboid(decompose_block_skin_halves_of_cuboid(SNAP_TO_NODE), 1.e-5);
}

// The same, with the middle face split along different diagonals by the two halves: the surfaces
// coincide there without matching facets.
TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedCuboidWithoutMatchingFacets_snapping_exactVolumesAndAreas)
{
  expect_halves_of_cuboid(decompose_halves_of_cuboid(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, false), 1.e-10);
}

// Only cutting, with the sharp features cut first: as accurate as cutting only each half alone (or the
// rotated cuboid, RotatedCuboidDecomposition), with no element of a half given the elsewhere phase.
TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedCuboid_cuttingOnly_volumesAndAreasAsForOneSurface)
{
  expect_halves_of_cuboid(decompose_halves_of_cuboid(SNAP_TO_NODE, true), 1.e-5);
}

// In 2D, the halves of a rotated rectangle meeting along the middle segment: exact areas and lengths
// with snapping, with matching middle segments and with the top half splitting it in two.
TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedRectangle_snapping_exactAreasAndLengths)
{
  expect_halves_of_rectangle(decompose_halves_of_rectangle(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, true), 1.e-10);
}

TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedRectangleWithoutMatchingSegments_snapping_exactAreasAndLengths)
{
  expect_halves_of_rectangle(decompose_halves_of_rectangle(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, false), 1.e-10);
}

TEST(MatingFacetedSurfacesDecomposition, halvesOfRotatedRectangle_cuttingOnly_areasAndLengthsAsForOneSurface)
{
  expect_halves_of_rectangle(decompose_halves_of_rectangle(SNAP_TO_NODE, true), 1.e-5);
}

// Mating halves that have not moved do not need the mesh decomposed again: a node where the halves
// mate is on both surfaces, though the phase of each half has a sign only for its own surface.
TEST(MatingFacetedSurfacesDecomposition, unchangedHalvesOfRotatedCuboid_noUpdateNeeded)
{
  EXPECT_FALSE(does_decomposition_of_halves_need_update(false, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE));
  EXPECT_FALSE(does_decomposition_of_halves_need_update(false, SNAP_TO_NODE));
}

TEST(MatingFacetedSurfacesDecomposition, unchangedHalvesOfRotatedRectangle_noUpdateNeeded)
{
  EXPECT_FALSE(does_decomposition_of_halves_need_update(true, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE));
  EXPECT_FALSE(does_decomposition_of_halves_need_update(true, SNAP_TO_NODE));
}

} // namespace krino

