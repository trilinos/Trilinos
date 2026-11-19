// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SecondOrderMeshSpecs.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace krino
{

// A mesh decomposed for a level set is decomposed again after the level set changes only at the
// nodes the decomposition added (on the interface and, for quadratic elements, the mid-side nodes
// of the subelements), as when an application solves for the level set on the decomposed mesh.  The
// new interface must be where the level set, as it is on the decomposed mesh, is zero: the
// crossings of an edge of the mesh are found along its child nodes, not from the level set at its
// ends, which has not changed.
template <typename MESHSPEC>
class RedecompositionFixture : public DecompositionFixture<MESHSPEC, LSPerInterfacePolicy, 1>
{
public:
  using Base = DecompositionFixture<MESHSPEC, LSPerInterfacePolicy, 1>;
  using LevelSetFunction = std::function<double(const stk::math::Vector3d &)>;

  RedecompositionFixture()
  {
    this->build_mesh_from_spec(this->meshSpec);
    this->setup_ls_fields();
    stk::mesh::get_entities(this->mMesh, stk::topology::NODE_RANK, myMeshNodes);
    std::sort(myMeshNodes.begin(), myMeshNodes.end());
  }

protected:
  // Set the level set at every node, or only at the nodes the decompositions added.
  void set_level_set_at_nodes(const bool onlyAddedNodes, const LevelSetFunction & levelSet)
  {
    stk::mesh::BulkData & mesh = this->mMesh;
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    const FieldRef coordsField = this->get_coordinates_field();
    const FieldRef lsField = this->levelset_fields()[0].isovar;
    auto coordsData = coordsField.template data<double>();
    auto lsData = lsField.template data<double, stk::mesh::ReadWrite>();
    for (auto && bucket : mesh.buckets(stk::topology::NODE_RANK))
    {
      for (auto node : *bucket)
      {
        if (onlyAddedNodes && std::binary_search(myMeshNodes.begin(), myMeshNodes.end(), node))
          continue;
        if (has_field_data(lsField, node))
          lsData.entity_values(node)() =
              levelSet(stk::math::Vector3d(coordsData.entity_values(node), dim));
      }
    }
  }

  void decompose(const int stepCount)
  {
    stk::mesh::BulkData & mesh = this->mMesh;
    const std::unique_ptr<InterfaceGeometry> geometry =
        create_levelset_geometry(mesh.mesh_meta_data().spatial_dimension(),
            this->get_aux_meta().active_part(),
            this->cdfem_support(),
            Phase_Support::get(mesh.mesh_meta_data()),
            this->levelset_fields());
    CDMesh::decompose_mesh(mesh, *geometry, stepCount);
  }

  double phase_volume(const std::string & phaseName)
  {
    const stk::mesh::Selector phaseSelector =
        this->get_aux_meta().get_part("block_1_" + phaseName) & this->get_aux_meta().active_part();
    return compute_volume_or_surface_area(this->mMesh, stk::topology::ELEMENT_RANK, phaseSelector);
  }

  // The right triangle or tetrahedron, with vertices at the origin and at unit distance along each
  // axis, is decomposed, only cutting, for the level set x-0.3, which crosses the edges to the
  // vertex at x=1 at x=0.3.  Then the level set is changed at the nodes the decomposition added to
  //   h(x) = (x-0.18)/0.6      for x <= 0.3,
  //   h(x) = 0.2+(x-0.3)/1.4   for x >= 0.3,
  // which is x-0.3 at the vertices (x=0 and x=1), and linear in every subelement (each is on one
  // side of x=0.3): it is +0.2 at the old interface nodes, so that the new interface is at x=0.18.
  // From the level set at the ends of the edges alone, it would be at x=0.3 again.
  void decompose_cut_only_then_again_after_changing_level_set_at_child_nodes()
  {
    if (stk::parallel_machine_size(this->mComm) > 1) return;

    this->cdfem_support().set_cdfem_edge_degeneracy_handling(SNAP_TO_NODE);

    set_level_set_at_nodes(false, [](const stk::math::Vector3d & x) { return x[0] - 0.3; });
    decompose(0);

    const double totalVolume = phase_volume("P-") + phase_volume("P+");
    expect_negative_phase_volume_where_x_is_less_than(0.3, totalVolume);

    const auto h = [](const stk::math::Vector3d & x)
    { return (x[0] <= 0.3) ? (x[0] - 0.18) / 0.6 : 0.2 + (x[0] - 0.3) / 1.4; };
    set_level_set_at_nodes(true, h);
    decompose(1);

    expect_negative_phase_volume_where_x_is_less_than(0.18, totalVolume);
  }

private:
  std::vector<stk::mesh::Entity> myMeshNodes;

  void expect_negative_phase_volume_where_x_is_less_than(const double x0, const double totalVolume)
  {
    const unsigned dim = this->mMesh.mesh_meta_data().spatial_dimension();
    // The part of the right simplex beyond x0 is the right simplex scaled by 1-x0.
    const double expectedVolume = totalVolume * (1. - std::pow(1. - x0, dim));
    EXPECT_NEAR(expectedVolume, phase_volume("P-"), 1.e-12);
    EXPECT_NEAR(totalVolume - expectedVolume, phase_volume("P+"), 1.e-12);
  }
};

using RightTri6 = HigherOrderSimplexMeshSpec<RightTri, stk::topology::TRIANGLE_6_2D>;

using RedecompositionTri3 = RedecompositionFixture<RightTri>;
using RedecompositionTri6 = RedecompositionFixture<RightTri6>;
using RedecompositionTet4 = RedecompositionFixture<RightTet>;
using RedecompositionTet10 = RedecompositionFixture<RightTet10>;

TEST_F(RedecompositionTri3, cutOnlyAfterLevelSetChangesAtChildNodes_interfaceAlongChildNodes)
{
  decompose_cut_only_then_again_after_changing_level_set_at_child_nodes();
}

TEST_F(RedecompositionTri6, cutOnlyAfterLevelSetChangesAtChildNodes_interfaceAlongChildNodes)
{
  decompose_cut_only_then_again_after_changing_level_set_at_child_nodes();
}

TEST_F(RedecompositionTet4, cutOnlyAfterLevelSetChangesAtChildNodes_interfaceAlongChildNodes)
{
  decompose_cut_only_then_again_after_changing_level_set_at_child_nodes();
}

TEST_F(RedecompositionTet10, cutOnlyAfterLevelSetChangesAtChildNodes_interfaceAlongChildNodes)
{
  decompose_cut_only_then_again_after_changing_level_set_at_child_nodes();
}

// The same, on a patch of elements: the unit square or cube, meshed with cells of size 0.25 split
// into simplices (no nodes inside the cells), so that every vertex of the mesh is on a plane
// x=0.25*k and every element is within a slab between two of them.  The second decomposition
// snaps, as well as cuts, so that vertices move onto the interface; in a patch of elements they are
// inside it or on its flat walls, so that the domain does not change.
class RedecompositionPatchOfElements
{
public:
  explicit RedecompositionPatchOfElements(const stk::topology elemTopology)
      : myBBoxMesh(elemTopology, MPI_COMM_WORLD)
  {
    stk::mesh::MetaData & meta = myBBoxMesh.meta_data();
    myLSFields = LSPerInterfacePolicy::setup_levelsets_on_all_blocks(meta, 1);

    CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
    cdfemSupport.set_coords_field(meta.coordinate_field());
    cdfemSupport.register_cdfem_snap_displacements_field();
    cdfemSupport.register_parent_node_ids_field();
    for (auto && lsField : myLSFields)
      cdfemSupport.add_interpolation_field(lsField.isovar);
    cdfemSupport.finalize_fields();

    myBBoxMesh.set_mesh_structure_type(CUBIC_BOUNDING_BOX_MESH);
    const bool is3D = elemTopology.dimension() == 3;
    populate_bounding_box_mesh_and_activate(
        myBBoxMesh, stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., is3D ? 1. : 0.), 0.25);

    stk::mesh::get_entities(myBBoxMesh.bulk_data(), stk::topology::NODE_RANK, myMeshVertices);
    const std::vector<stk::math::Vector3d> coords = mesh_vertex_coordinates();
    size_t numVertices = 0;
    for (size_t i = 0; i < myMeshVertices.size(); ++i)
      if (is_on_lattice(coords[i])) myMeshVertices[numVertices++] = myMeshVertices[i];
    myMeshVertices.resize(numVertices);
  }

  // Decompose for the level set x-0.4, cutting only, or snapping as well.  Then change the level
  // set at the active nodes (those of the decomposed mesh, where an application would solve for it)
  // to
  //   h(x) = x-0.4                     outside 0.25 <= x <= 0.5,
  //   h(x) = -0.15+(x-0.25)/0.6        for 0.25 <= x <= 0.4,
  //   h(x) = 0.1                       for 0.4 <= x <= 0.5,
  // which is x-0.4 on every plane x=0.25*k, so that only the nodes on the old interface (cut or
  // snapped) and the mid-side nodes of the subelements change, and is linear in every subelement
  // (each is on one side of x=0.4).  Its zero is at x=0.34: decomposed again, snapping and cutting,
  // the interface must be there, and the vertices that snap must snap onto it.  From the level set
  // at the ends of the edges alone, the crossings would be at x=0.4 again.  The phase volumes would
  // still be right, as h is linear in every subelement, but snapping would move vertices onto
  // x=0.4, off the interface.  This needs the edges cut by the first decomposition: one that a
  // snapped vertex ends has no child node.
  void
  decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(const bool snapFirst)
  {
    set_level_set_at_nodes(myBBoxMesh.meta_data().universal_part(),
        [](const stk::math::Vector3d & x) { return x[0] - 0.4; });
    expect_decomposition_gives_negative_phase_volume(0, snapFirst, 0.4);

    const auto h = [](const stk::math::Vector3d & x)
    {
      if (x[0] <= 0.25 || x[0] >= 0.5) return x[0] - 0.4;
      if (x[0] <= 0.4) return -0.15 + (x[0] - 0.25) / 0.6;
      return 0.1;
    };
    set_level_set_at_nodes(AuxMetaData::get(myBBoxMesh.meta_data()).active_part(), h);
    expect_decomposition_gives_negative_phase_volume(1, true, 0.34);
  }

private:
  using LevelSetFunction = std::function<double(const stk::math::Vector3d &)>;

  static bool is_on_lattice(const stk::math::Vector3d & x)
  {
    for (unsigned d = 0; d < 3; ++d)
      if (std::abs(4. * x[d] - std::round(4. * x[d])) > 1.e-9) return false;
    return true;
  }

  void set_level_set_at_nodes(
      const stk::mesh::Selector & nodeSelector, const LevelSetFunction & levelSet)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
    auto coordsData = coordsField.data<double>();
    auto lsData = myLSFields[0].isovar.data<double, stk::mesh::ReadWrite>();
    for (auto && bucket : mesh.get_buckets(stk::topology::NODE_RANK, nodeSelector))
      for (auto node : *bucket)
        lsData.entity_values(node)() =
            levelSet(stk::math::Vector3d(coordsData.entity_values(node), dim));
  }

  std::vector<stk::math::Vector3d> mesh_vertex_coordinates() const
  {
    const stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
    auto coordsData = coordsField.data<double>();
    std::vector<stk::math::Vector3d> coords;
    for (auto node : myMeshVertices)
      coords.emplace_back(coordsData.entity_values(node), dim);
    return coords;
  }

  // The numbers of mesh vertices (owned) that moved, and of those that did not move onto the plane
  // x=x0 of the interface: snapping moves a vertex onto an intersection point of the interface.
  std::array<size_t, 2> num_mesh_vertices_moved_and_moved_off_interface(
      const std::vector<stk::math::Vector3d> & coordsBefore, const double x0) const
  {
    const stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    const std::vector<stk::math::Vector3d> coordsAfter = mesh_vertex_coordinates();
    size_t numMoved = 0;
    size_t numMovedOffInterface = 0;
    for (size_t i = 0; i < myMeshVertices.size(); ++i)
    {
      if (mesh.bucket(myMeshVertices[i]).owned() &&
          (coordsAfter[i] - coordsBefore[i]).length() > 1.e-12)
      {
        ++numMoved;
        if (std::abs(coordsAfter[i][0] - x0) > 1.e-12) ++numMovedOffInterface;
      }
    }
    all_reduce_sum(mesh.parallel(), numMoved, numMovedOffInterface);
    return {numMoved, numMovedOffInterface};
  }

  double phase_volume(const std::string & phaseName) const
  {
    const stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
    const stk::mesh::Selector phaseSelector = auxMeta.get_part("block_1_" + phaseName) &
        auxMeta.active_part() & mesh.mesh_meta_data().locally_owned_part();
    double volume =
        compute_volume_or_surface_area(mesh, stk::topology::ELEMENT_RANK, phaseSelector);
    all_reduce_sum(mesh.parallel(), volume);
    return volume;
  }

  void expect_decomposition_gives_negative_phase_volume(
      const int stepCount, const bool snap, const double x0)
  {
    stk::mesh::BulkData & mesh = myBBoxMesh.bulk_data();
    CDFEM_Support::get(mesh.mesh_meta_data())
        .set_cdfem_edge_degeneracy_handling(
            snap ? SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE : SNAP_TO_NODE);
    const std::vector<stk::math::Vector3d> coordsBefore = mesh_vertex_coordinates();

    const std::unique_ptr<InterfaceGeometry> geometry =
        create_levelset_geometry(mesh.mesh_meta_data().spatial_dimension(),
            AuxMetaData::get(mesh.mesh_meta_data()).active_part(),
            CDFEM_Support::get(mesh.mesh_meta_data()),
            Phase_Support::get(mesh.mesh_meta_data()),
            myLSFields);
    CDMesh::decompose_mesh(mesh, *geometry, stepCount);

    const auto [numMoved, numMovedOffInterface] =
        num_mesh_vertices_moved_and_moved_off_interface(coordsBefore, x0);
    if (snap)
    {
      EXPECT_LT(0u, numMoved) << "decomposition " << stepCount;
    }
    else
    {
      EXPECT_EQ(0u, numMoved) << "decomposition " << stepCount;
    }
    EXPECT_EQ(0u, numMovedOffInterface) << "decomposition " << stepCount;
    EXPECT_NEAR(x0, phase_volume("P-"), 1.e-12) << "decomposition " << stepCount;
    EXPECT_NEAR(1. - x0, phase_volume("P+"), 1.e-12) << "decomposition " << stepCount;
  }

  BoundingBoxMesh myBBoxMesh;
  std::vector<LS_Field> myLSFields;
  std::vector<stk::mesh::Entity> myMeshVertices;
};

// Cut only, then snapping: the edges the interface crosses all have child nodes when snapping.
TEST(RedecompositionPatchOfElementsTri3,
    cutThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterfaceAlongChildNodes)
{
  RedecompositionPatchOfElements(stk::topology::TRIANGLE_3_2D)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(false);
}

TEST(RedecompositionPatchOfElementsTri6,
    cutThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterfaceAlongChildNodes)
{
  RedecompositionPatchOfElements(stk::topology::TRIANGLE_6_2D)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(false);
}

TEST(RedecompositionPatchOfElementsTet4,
    cutThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterfaceAlongChildNodes)
{
  RedecompositionPatchOfElements(stk::topology::TETRAHEDRON_4)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(false);
}

TEST(RedecompositionPatchOfElementsTet10,
    cutThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterfaceAlongChildNodes)
{
  RedecompositionPatchOfElements(stk::topology::TETRAHEDRON_10)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(false);
}

// Snapping both times.
TEST(RedecompositionPatchOfElementsTri3,
    snapThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterface)
{
  RedecompositionPatchOfElements(stk::topology::TRIANGLE_3_2D)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(true);
}

TEST(RedecompositionPatchOfElementsTri6,
    snapThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterface)
{
  RedecompositionPatchOfElements(stk::topology::TRIANGLE_6_2D)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(true);
}

TEST(RedecompositionPatchOfElementsTet4,
    snapThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterface)
{
  RedecompositionPatchOfElements(stk::topology::TETRAHEDRON_4)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(true);
}

TEST(RedecompositionPatchOfElementsTet10,
    snapThenSnapAfterLevelSetChangesAtActiveNodes_snapsOntoInterface)
{
  RedecompositionPatchOfElements(stk::topology::TETRAHEDRON_10)
      .decompose_then_again_with_snapping_after_changing_level_set_at_active_nodes(true);
}

} // namespace krino
