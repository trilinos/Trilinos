// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// A decomposition must not depend on how the mesh is decomposed across processors, and, when its
// global ids are not parallel consistent (CDFEM_Support::set_global_ids_are_NOT_parallel_consistent,
// as when refinement assigns them, which changes them with the decomposition and the build), on its
// ids either.  Each test here decomposes the same mesh six ways, and requires the same decomposition,
// element for element, within each of two groups:
// - with the global ids parallel consistent (the default), where ties may be broken by ids: as
//   generated, and decomposed across the processors differently;
// - with them not parallel consistent: as generated, with its ids reversed, decomposed across the
//   processors differently, and both.
// On one processor only the ids change.

#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_Unit_DecompositionSignature.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace krino
{

namespace
{

struct MeshVariation
{
  bool globalIdsAreParallelConsistent;
  bool reverseEntityIds;
  bool alternateProcessorDecomposition;
};

// Within each group, the first is the one the others are compared with.  Ids are varied only where
// they are not declared parallel consistent.
const std::array<std::vector<MeshVariation>, 2> invarianceMeshVariationGroups{{
    {{true, false, false}, {true, false, true}},
    {{false, false, false}, {false, true, false}, {false, false, true}, {false, true, true}}}};

// Whether each level set (or surface) is an interface between its negative and positive sides, or
// each is a phase, where it is the smallest.
enum class InvarianceLevelSetPolicy { PER_INTERFACE, PER_PHASE };

// A generated mesh of a box (a flat-walled BCC tet mesh, or a flat-walled triangular lattice in 2D),
// with level sets as given, set up for decomposition, and generated with the given variation of its
// ids and its decomposition across processors.
class VariedMeshDecomposition
{
public:
  VariedMeshDecomposition(const stk::topology elemTopology,
      const BoundingBoxMesh::BoundingBoxType & domain,
      const double meshSize,
      const unsigned numSurfaces,
      const InvarianceLevelSetPolicy levelSetPolicy,
      const Edge_Degeneracy_Handling edgeDegeneracyHandling,
      const MeshVariation & variation)
      : myBBoxMesh(elemTopology, MPI_COMM_WORLD)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    myLSFields = (levelSetPolicy == InvarianceLevelSetPolicy::PER_PHASE) ?
        LSPerPhasePolicy::setup_levelsets_on_all_blocks(meta, numSurfaces) :
        LSPerInterfacePolicy::setup_levelsets_on_all_blocks(meta, numSurfaces);

    CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
    cdfemSupport.set_coords_field(meta.coordinate_field());
    cdfemSupport.set_cdfem_edge_degeneracy_handling(edgeDegeneracyHandling);
    cdfemSupport.register_cdfem_snap_displacements_field();
    cdfemSupport.register_parent_node_ids_field();
    cdfemSupport.finalize_fields();
    if (!variation.globalIdsAreParallelConsistent)
      cdfemSupport.set_global_ids_are_NOT_parallel_consistent();

    myBBoxMesh.set_mesh_structure_type(elemTopology.base() == stk::topology::TETRAHEDRON_4 ?
            FLAT_WALLED_BCC_BOUNDING_BOX_MESH :
            FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH);
    myBBoxMesh.set_reverse_entity_ids(variation.reverseEntityIds);
    myBBoxMesh.set_alternate_processor_decomposition(variation.alternateProcessorDecomposition);
    populate_bounding_box_mesh_and_activate(
        myBBoxMesh, domain.get_min(), domain.get_max(), meshSize);
  }

  std::vector<ElementSignature> decompose_and_gather_signature(const std::vector<const Surface *> & surfaces)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    AuxMetaData & auxMeta = AuxMetaData::get(meta);

    std::vector<Surface_Identifier> surfaceIdentifiers;
    for (size_t i = 0; i < surfaces.size(); ++i)
      surfaceIdentifiers.push_back(myLSFields[i].identifier);
    const AnalyticSurfaceInterfaceGeometry geometry(
        surfaceIdentifiers, surfaces, auxMeta.active_part(), CDFEM_Support::get(meta), Phase_Support::get(meta));
    CDMesh::decompose_mesh(mesh, geometry);

    return gather_decomposition_signature(mesh, auxMeta.active_part());
  }

  // The level sets given as functions of position, at every node, decomposed.
  std::vector<ElementSignature> decompose_level_sets_and_gather_signature(
      const std::vector<std::function<double(const stk::math::Vector3d &)>> & levelSets)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    AuxMetaData & auxMeta = AuxMetaData::get(meta);

    const FieldRef coordsField(meta.coordinate_field());
    const unsigned dim = meta.spatial_dimension();
    for (const stk::mesh::Bucket * bucket : mesh.buckets(stk::topology::NODE_RANK))
      for (const stk::mesh::Entity node : *bucket)
        for (size_t i = 0; i < levelSets.size(); ++i)
          get_scalar_field(mesh, myLSFields[i].isovar, node) = levelSets[i](get_vector_field(mesh, coordsField, node, dim));

    const std::unique_ptr<InterfaceGeometry> geometry = create_levelset_geometry(
        dim, auxMeta.active_part(), CDFEM_Support::get(meta), Phase_Support::get(meta), myLSFields);
    CDMesh::decompose_mesh(mesh, *geometry);

    return gather_decomposition_signature(mesh, auxMeta.active_part());
  }

private:
  BoundingBoxMesh myBBoxMesh;
  std::vector<LS_Field> myLSFields;
};

// The number of different sets of element parts in the signature (on every processor): at least
// two when the surfaces cut the mesh into more than one phase.
size_t global_num_element_part_sets(stk::ParallelMachine comm, const std::vector<ElementSignature> & signature)
{
  std::set<std::vector<std::string>> partSets;
  for (const ElementSignature & element : signature)
    partSets.insert(element.partNames);
  const size_t localNum = partSets.size();
  size_t globalNum = 0;
  stk::all_reduce_max(comm, &localNum, &globalNum, 1);
  return globalNum;
}

void expect_decomposition_of_each_variation_the_same(const stk::topology elemTopology,
    const BoundingBoxMesh::BoundingBoxType & domain,
    const double meshSize,
    const unsigned numLevelSets,
    const InvarianceLevelSetPolicy levelSetPolicy,
    const Edge_Degeneracy_Handling edgeDegeneracyHandling,
    const std::function<std::vector<ElementSignature>(VariedMeshDecomposition &)> & decompose_and_gather_signature)
{
  for (const std::vector<MeshVariation> & variations : invarianceMeshVariationGroups)
  {
    std::vector<ElementSignature> firstSignature;
    for (const MeshVariation & variation : variations)
    {
      VariedMeshDecomposition decomposition(
          elemTopology, domain, meshSize, numLevelSets, levelSetPolicy, edgeDegeneracyHandling, variation);
      const std::vector<ElementSignature> signature = decompose_and_gather_signature(decomposition);
      if (&variation == &variations.front())
      {
        firstSignature = signature;
        ASSERT_LE(2u, global_num_element_part_sets(MPI_COMM_WORLD, signature));
        continue;
      }
      EXPECT_EQ(0u, global_num_elements_differing_between_signatures(MPI_COMM_WORLD, firstSignature, signature))
          << "global ids parallel consistent " << variation.globalIdsAreParallelConsistent << ", reversed ids "
          << variation.reverseEntityIds << ", alternate processor decomposition "
          << variation.alternateProcessorDecomposition << "\n"
          << describe_elements_differing_between_signatures(firstSignature, signature, 20);
    }
  }
}

void expect_decomposition_independent_of_ids_and_processor_decomposition(const stk::topology elemTopology,
    const BoundingBoxMesh::BoundingBoxType & domain,
    const double meshSize,
    const Edge_Degeneracy_Handling edgeDegeneracyHandling,
    const std::vector<const Surface *> & surfaces)
{
  expect_decomposition_of_each_variation_the_same(elemTopology, domain, meshSize, surfaces.size(),
      InvarianceLevelSetPolicy::PER_INTERFACE, edgeDegeneracyHandling,
      [&surfaces](VariedMeshDecomposition & decomposition)
      { return decomposition.decompose_and_gather_signature(surfaces); });
}

const BoundingBoxMesh::BoundingBoxType invarianceCubeDomain(stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.));
const BoundingBoxMesh::BoundingBoxType invarianceSquareDomain(stk::math::Vector3d(-1., -1., 0.), stk::math::Vector3d(1., 1., 0.));
constexpr double invarianceMeshSize = 0.3;

Cuboid rotated_cuboid_for_invariance_test()
{
  return Cuboid(stk::math::Vector3d(0.03, -0.02, 0.05), stk::math::Vector3d(1.0, 0.7, 0.8), stk::math::Vector3d(0.3, 0.5, 0.2));
}

Cuboid rotated_rectangle_for_invariance_test()
{
  return Cuboid(stk::math::Vector3d(0.03, -0.02, 0.), stk::math::Vector3d(1.0, 0.7, 10.), stk::math::Vector3d(0., 0., 0.3));
}

// A slab thin against the mesh (its two faces cross many mesh edges close together, and some
// within the edge tolerance of each other), tilted to the mesh, as a connected faceted surface
// without sharp features.  The facets are all owned by processor 0.
std::unique_ptr<MeshFacetedSurface> tilted_thin_faceted_slab_for_invariance_test()
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 0.3, 0.2).unit_vector();
  const stk::math::Vector3d u = Cross(normal, stk::math::Vector3d(0., 0., 1.)).unit_vector();
  const stk::math::Vector3d v = Cross(normal, u);
  const double lo = 0.05;
  const double hi = 0.09;
  constexpr double w = 5.;
  const auto point = [&](const double x, const double y, const double z) { return x * normal + y * u + z * v; };
  const std::vector<stk::math::Vector3d> vertices{point(lo, -w, -w),
      point(hi, -w, -w),
      point(hi, w, -w),
      point(lo, w, -w),
      point(lo, -w, w),
      point(hi, -w, w),
      point(hi, w, w),
      point(lo, w, w)};
  std::vector<std::array<unsigned, 3>> facets;
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
    facets = {{{0, 4, 7}},
        {{0, 7, 3}},
        {{1, 2, 6}},
        {{1, 6, 5}},
        {{0, 1, 5}},
        {{0, 5, 4}},
        {{2, 3, 7}},
        {{2, 7, 6}},
        {{0, 3, 2}},
        {{0, 2, 1}},
        {{4, 5, 6}},
        {{4, 6, 7}}};
  return std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, vertices, facets);
}

} // namespace

// The sharp features of a rotated cuboid, cut by rank before the surface, then snapping.
TEST(DecompositionInvariance, snappingRotatedCuboid_independentOfIdsAndProcessorDecomposition)
{
  const Cuboid cuboid = rotated_cuboid_for_invariance_test();
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TETRAHEDRON_4,
      invarianceCubeDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {&cuboid});
}

TEST(DecompositionInvariance, cuttingOnlyRotatedCuboid_independentOfIdsAndProcessorDecomposition)
{
  const Cuboid cuboid = rotated_cuboid_for_invariance_test();
  expect_decomposition_independent_of_ids_and_processor_decomposition(
      stk::topology::TETRAHEDRON_4, invarianceCubeDomain, invarianceMeshSize, SNAP_TO_NODE, {&cuboid});
}

TEST(DecompositionInvariance, snappingRotatedCuboidOnTet10_independentOfIdsAndProcessorDecomposition)
{
  const Cuboid cuboid = rotated_cuboid_for_invariance_test();
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TETRAHEDRON_10,
      invarianceCubeDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {&cuboid});
}

TEST(DecompositionInvariance, snappingRotatedRectangle_independentOfIdsAndProcessorDecomposition)
{
  const Cuboid rectangle = rotated_rectangle_for_invariance_test();
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TRIANGLE_3_2D,
      invarianceSquareDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {&rectangle});
}

TEST(DecompositionInvariance, snappingRotatedRectangleOnTri6_independentOfIdsAndProcessorDecomposition)
{
  const Cuboid rectangle = rotated_rectangle_for_invariance_test();
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TRIANGLE_6_2D,
      invarianceSquareDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {&rectangle});
}

// The multiple crossings of a surface without sharp features, with close crossings: the order in
// which an element's crossings are cut must not depend on the ids.
TEST(DecompositionInvariance, cuttingOnlyThinFacetedSlab_independentOfIdsAndProcessorDecomposition)
{
  const auto slab = tilted_thin_faceted_slab_for_invariance_test();
  ASSERT_FALSE(slab->has_sharp_features());
  expect_decomposition_independent_of_ids_and_processor_decomposition(
      stk::topology::TETRAHEDRON_4, invarianceCubeDomain, invarianceMeshSize, SNAP_TO_NODE, {slab.get()});
}

TEST(DecompositionInvariance, snappingThinFacetedSlab_independentOfIdsAndProcessorDecomposition)
{
  const auto slab = tilted_thin_faceted_slab_for_invariance_test();
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TETRAHEDRON_4,
      invarianceCubeDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {slab.get()});
}

// A smooth surface with nodes snapped onto it, and crossings beside them.
TEST(DecompositionInvariance, snappingSphere_independentOfIdsAndProcessorDecomposition)
{
  const Sphere sphere(stk::math::Vector3d(0.04, -0.03, 0.02), 0.62);
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TETRAHEDRON_4,
      invarianceCubeDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {&sphere});
}

// Two surfaces crossing each other, along a line where four phases meet.
TEST(DecompositionInvariance, snappingTwoCrossingSpheres_independentOfIdsAndProcessorDecomposition)
{
  const Sphere sphere0(stk::math::Vector3d(-0.21, -0.03, 0.02), 0.55);
  const Sphere sphere1(stk::math::Vector3d(0.24, 0.05, -0.04), 0.5);
  expect_decomposition_independent_of_ids_and_processor_decomposition(stk::topology::TETRAHEDRON_4,
      invarianceCubeDomain, invarianceMeshSize, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE, {&sphere0, &sphere1});
}

namespace
{
// Two overlapping spheres as phases, each where its level set (the distance to it) is the smallest,
// in a background phase whose level set is 0 everywhere: the phases meet, three at a time, along the
// circle where the spheres cross, as particles touching in a binder.
void expect_level_set_phases_of_two_spheres_independent_of_ids_and_processor_decomposition(
    const stk::topology elemTopology, const Edge_Degeneracy_Handling edgeDegeneracyHandling)
{
  const auto sphere_distance = [](const stk::math::Vector3d & center, const double radius)
  { return [center, radius](const stk::math::Vector3d & x) { return (x - center).length() - radius; }; };
  const std::vector<std::function<double(const stk::math::Vector3d &)>> levelSets{
      sphere_distance(stk::math::Vector3d(-0.21, -0.03, 0.02), 0.55),
      sphere_distance(stk::math::Vector3d(0.24, 0.05, -0.04), 0.5),
      [](const stk::math::Vector3d &) { return 0.; }};
  expect_decomposition_of_each_variation_the_same(elemTopology, invarianceCubeDomain, invarianceMeshSize,
      levelSets.size(), InvarianceLevelSetPolicy::PER_PHASE, edgeDegeneracyHandling,
      [&levelSets](VariedMeshDecomposition & decomposition)
      { return decomposition.decompose_level_sets_and_gather_signature(levelSets); });
}

} // namespace

// Level sets as phases, as the sphere-and-binder problems are set up, cut only, on tet4 and tet10, and
// snapped.
TEST(DecompositionInvariance, cuttingOnlyLevelSetPhasesOfTwoSpheres_independentOfIdsAndProcessorDecomposition)
{
  expect_level_set_phases_of_two_spheres_independent_of_ids_and_processor_decomposition(
      stk::topology::TETRAHEDRON_4, SNAP_TO_NODE);
}

TEST(DecompositionInvariance, cuttingOnlyLevelSetPhasesOfTwoSpheresOnTet10_independentOfIdsAndProcessorDecomposition)
{
  expect_level_set_phases_of_two_spheres_independent_of_ids_and_processor_decomposition(
      stk::topology::TETRAHEDRON_10, SNAP_TO_NODE);
}

TEST(DecompositionInvariance, snappingLevelSetPhasesOfTwoSpheres_independentOfIdsAndProcessorDecomposition)
{
  expect_level_set_phases_of_two_spheres_independent_of_ids_and_processor_decomposition(
      stk::topology::TETRAHEDRON_4, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
}

TEST(DecompositionInvariance, snappingLevelSetPhasesOfTwoSpheresOnTet10_independentOfIdsAndProcessorDecomposition)
{
  expect_level_set_phases_of_two_spheres_independent_of_ids_and_processor_decomposition(
      stk::topology::TETRAHEDRON_10, SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
}

} // namespace krino
