// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_SecondOrderMeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_Unit_BoundingBoxMesh.hpp>
#include <Akri_UnitMeshUtils.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace krino
{

namespace
{

const BoundingBox sidesetQueryBox(stk::math::Vector3d(-3, -3, -3), stk::math::Vector3d(3, 3, 3));

size_t global_sum(const size_t local)
{
  size_t global = 0;
  stk::all_reduce_sum(MPI_COMM_WORLD, &local, &global, 1);
  return global;
}

// Sets the coordinates of every local copy of the nodes at from to to.
void move_nodes_at(stk::mesh::BulkData & mesh,
    stk::mesh::Field<double> & coords,
    const stk::math::Vector3d & from,
    const stk::math::Vector3d & to)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  auto coordsData = coords.data<stk::mesh::ReadWrite>();
  for (auto * bucket : mesh.buckets(stk::topology::NODE_RANK))
    for (auto node : *bucket)
    {
      auto x = coordsData.entity_values(node);
      if ((stk::math::Vector3d(x, dim) - from).length() < 1.e-12)
        for (unsigned d = 0; d < dim; ++d)
          x(stk::mesh::ComponentIdx(d)) = to[d];
    }
}

stk::mesh::Field<double> & coordinates_field(const stk::mesh::MetaData & meta)
{
  return *static_cast<stk::mesh::Field<double> *>(
      &AuxMetaData::get(meta).get_current_coordinates().field());
}

struct TwoTetMeshOptions
{
  // surface_1 touches block_2, which has no elements, instead of block_1.
  bool mapSurfaceToEmptyBlock = false;
  // The sides are also in a second side part, surface_2, touching block_1.
  bool sidesAlsoInSecondSurface = false;
};

// Two tets in block_1 sharing a face, with their exposed boundary sides in a side part "surface_1"
// touching block_1 (or the empty block_2), and optionally also in "surface_2".
struct TwoTetMeshWithBoundarySideset
{
  explicit TwoTetMeshWithBoundarySideset(const TwoTetMeshOptions & options = TwoTetMeshOptions{})
  {
    myMeta = stk::mesh::MeshBuilder().set_spatial_dimension(3).create_meta_data();
    myBlock = &myMeta->declare_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
    stk::io::put_io_part_attribute(*myBlock);
    stk::mesh::Part & emptyBlock =
        myMeta->declare_part_with_topology("block_2", stk::topology::TETRAHEDRON_4);
    stk::io::put_io_part_attribute(emptyBlock);
    mySurface = &myMeta->declare_part_with_topology("surface_1", stk::topology::TRIANGLE_3);
    stk::io::put_io_part_attribute(*mySurface);
    myMeta->set_surface_to_block_mapping(
        mySurface, {options.mapSurfaceToEmptyBlock ? &emptyBlock : myBlock});
    mySecondSurface = &myMeta->declare_part_with_topology("surface_2", stk::topology::TRIANGLE_3);
    stk::io::put_io_part_attribute(*mySecondSurface);
    myMeta->set_surface_to_block_mapping(mySecondSurface, {myBlock});
    myCoords = &myMeta->declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(*myCoords, myMeta->universal_part(), 3, nullptr);
    myMeta->set_coordinate_field(myCoords);
    myMeta->commit();

    myMesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create(myMeta);

    const std::vector<stk::math::Vector3d> nodeLocs{
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}};
    const std::vector<stk::mesh::EntityIdVector> elemNodeIds{{1, 2, 3, 4}, {2, 3, 4, 5}};
    const int me = myMesh->parallel_rank();
    const int numProcs = myMesh->parallel_size();

    myMesh->modification_begin();
    for (size_t e = 0; e < elemNodeIds.size(); ++e)
      if (int(e % numProcs) == me)
        stk::mesh::declare_element(*myMesh, *myBlock, e + 1, elemNodeIds[e]);
    if (numProcs > 1 && me < 2)
      for (stk::mesh::EntityId id : {2u, 3u, 4u}) // the shared face's nodes
        myMesh->add_node_sharing(myMesh->get_entity(stk::topology::NODE_RANK, id), 1 - me);
    auto coordsData = myCoords->data<stk::mesh::ReadWrite>();
    for (stk::mesh::EntityId id = 1; id <= nodeLocs.size(); ++id)
    {
      const stk::mesh::Entity node = myMesh->get_entity(stk::topology::NODE_RANK, id);
      if (!myMesh->is_valid(node)) continue;
      auto coords = coordsData.entity_values(node);
      for (unsigned d = 0; d < 3; ++d)
        coords(stk::mesh::ComponentIdx(d)) = nodeLocs[id - 1][d];
    }
    myMesh->modification_end();

    stk::mesh::PartVector sideParts{mySurface};
    if (options.sidesAlsoInSecondSurface) sideParts.push_back(mySecondSurface);
    stk::mesh::create_exposed_block_boundary_sides(*myMesh, *myBlock, sideParts);
  }

  std::shared_ptr<stk::mesh::MetaData> myMeta;
  std::unique_ptr<stk::mesh::BulkData> myMesh;
  stk::mesh::Part * myBlock{nullptr};
  stk::mesh::Part * mySurface{nullptr};
  stk::mesh::Part * mySecondSurface{nullptr};
  stk::mesh::Field<double> * myCoords{nullptr};
};

// +1 if facet's vertex ids are the node ids of one of the sides in the same cyclic order, -1 if in
// reverse order, 0 if no side has its nodes.
int orientation_relative_to_side(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & sides,
    const FacetMesh & facets,
    const unsigned facet)
{
  std::array<FacetVertexId, 3> ids;
  for (unsigned n = 0; n < 3; ++n)
    ids[n] = facets.vertex_id(facets.facet_vertices(facet)[n]);
  for (auto side : sides)
  {
    const stk::mesh::Entity * sideNodes = mesh.begin_nodes(side);
    for (unsigned r = 0; r < 3; ++r)
    {
      if (ids[0] != mesh.identifier(sideNodes[r])) continue;
      if (ids[1] == mesh.identifier(sideNodes[(r + 1) % 3]) &&
          ids[2] == mesh.identifier(sideNodes[(r + 2) % 3]))
        return 1;
      if (ids[1] == mesh.identifier(sideNodes[(r + 2) % 3]) &&
          ids[2] == mesh.identifier(sideNodes[(r + 1) % 3]))
        return -1;
    }
  }
  return 0;
}

int sign_of(const double x) { return (x > 0.) ? 1 : ((x < 0.) ? -1 : 0); }

// For each point on the surface and its outward unit normal, the point just inside has the sign of
// insidePoint and the point just outside the opposite sign.
void expect_inside_and_outside_signs(const MeshSurface & surface,
    const stk::math::Vector3d & insidePoint,
    const std::vector<std::pair<stk::math::Vector3d, stk::math::Vector3d>> &
        pointsAndOutwardNormals)
{
  const double eps = 0.02;
  const int insideSign = sign_of(surface.point_signed_distance(insidePoint));
  EXPECT_NE(0, insideSign);
  for (auto & [point, normal] : pointsAndOutwardNormals)
  {
    EXPECT_EQ(insideSign, sign_of(surface.point_signed_distance(point - eps * normal)))
        << "inside " << point;
    EXPECT_EQ(-insideSign, sign_of(surface.point_signed_distance(point + eps * normal)))
        << "outside " << point;
  }
}

// The outward unit normal of the face (a,b,c) of the solid with interior point inside.
stk::math::Vector3d outward_unit_normal(const stk::math::Vector3d & a,
    const stk::math::Vector3d & b,
    const stk::math::Vector3d & c,
    const stk::math::Vector3d & inside)
{
  stk::math::Vector3d n = Cross(b - a, c - a).unit_vector();
  return (Dot(n, (a + b + c) / 3. - inside) > 0.) ? n : -n;
}

} // namespace

// The 6 exposed faces of the two tets share the mesh's 5 nodes.  The point inside (closest to the
// faces x=0, y=0 and z=0) and the point outside (closest to the face x=0) have opposite signs, and
// sign = -1 negates them.
TEST(MeshSurface, twoTetBoundarySideset_connectedFacetsAndExactSignedDistances)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;

  TwoTetMeshWithBoundarySideset twoTets;
  const stk::math::Vector3d inside(0.1, 0.2, 0.3);
  const stk::math::Vector3d outside(-0.5, 0.2, 0.2);

  double insideWithPositiveSign = 0.;
  for (const int sign : {1, -1})
  {
    MeshSurface surface(*twoTets.myMeta, *twoTets.myCoords, *twoTets.mySurface, sign);
    surface.prepare_to_compute(0., sidesetQueryBox, 0.);

    EXPECT_EQ(6u, global_sum(surface.num_owned_facets()));
    EXPECT_EQ(5u, surface.global_num_vertices());

    const double insideDistance = surface.point_signed_distance(inside);
    const double outsideDistance = surface.point_signed_distance(outside);
    EXPECT_DOUBLE_EQ(0.1, std::abs(insideDistance));
    EXPECT_DOUBLE_EQ(0.5, std::abs(outsideDistance));
    EXPECT_GT(0., insideDistance * outsideDistance);
    if (sign == 1)
      insideWithPositiveSign = insideDistance;
    else
    {
      EXPECT_DOUBLE_EQ(-insideWithPositiveSign, insideDistance);
    }
  }
}

// The facets are rebuilt from the current coordinates in every prepare_to_compute.
TEST(MeshSurface, nodeMovedBetweenPrepares_distancesFollowTheMovedSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;

  TwoTetMeshWithBoundarySideset twoTets;
  MeshSurface surface(*twoTets.myMeta, *twoTets.myCoords, *twoTets.mySurface, 1);
  const stk::math::Vector3d movedNode(2., 2., 2.);

  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_LT(0.5, std::abs(surface.point_signed_distance(movedNode)));

  move_nodes_at(*twoTets.myMesh, *twoTets.myCoords, stk::math::Vector3d(1., 1., 1.), movedNode);
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_DOUBLE_EQ(0., surface.point_signed_distance(movedNode));
}

// 2x2x2 hexes, with the corner (0,0,0) moved down so that the 3 wall quads around it are not
// planar: each of those is 4 triangles about its centroid, which is a vertex of its own (with an id
// that is no node's), and the other 21 are 2 triangles each.
TEST(MeshSurface, hexWithMovedCorner_nonPlanarQuadsSplitAboutTheirCentroids)
{
  BoundingBoxMeshHex8 hex;
  populate_bounding_box_mesh_and_activate(hex, {0., 0., 0.}, {1., 1., 1.}, 0.5);
  stk::mesh::Field<double> & coords = coordinates_field(hex.meta_data());
  stk::mesh::BulkData & mesh = hex.bulk_data();
  move_nodes_at(mesh, coords, stk::math::Vector3d(0., 0., 0.), stk::math::Vector3d(0., 0., -0.5));

  stk::mesh::PartVector walls;
  for (auto name : {"Xminus", "Xplus", "Yminus", "Yplus", "Zminus", "Zplus"})
    walls.push_back(hex.meta_data().get_part(name));
  MeshSurface surface(hex.meta_data(), coords, stk::mesh::selectUnion(walls), 1);
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);

  EXPECT_EQ(3u * 4u + 21u * 2u, global_sum(surface.num_owned_facets()));
  EXPECT_EQ(26u + 3u, surface.global_num_vertices()); // the 26 boundary nodes and 3 centroids
  const FacetMesh & facets = surface.facet_mesh();
  const unsigned dim = 3;
  auto coordsData = coords.data();
  for (unsigned v = 0; v < facets.num_owned_facet_vertices(); ++v)
  {
    const stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, facets.vertex_id(v));
    if (mesh.is_valid(node))
    {
      EXPECT_EQ(
          stk::math::Vector3d(coordsData.entity_values(node), dim), facets.vertex_coordinates(v))
          << "node " << facets.vertex_id(v);
    }
  }
  // The centroid of the bottom corner quad (z=0 at 3 corners, -0.5 at the moved one) is on the
  // surface.
  EXPECT_NEAR(0., surface.point_signed_distance(stk::math::Vector3d(0.25, 0.25, -0.125)), 1.e-14);
}

// The quadratic sides of a tet10 (TRI_6) and a tri6 (LINE_3) are split at their mid nodes: 4
// triangles per face and 2 segments per edge, using all the element's nodes.
TEST_F(StkMeshTet10Fixture, meshSurfaceOfQuadraticSides_splitAtMidNodes)
{
  if (stk::parallel_machine_size(get_comm()) > 1) return;
  build_mesh_from_spec(RegularTet10());
  mMesh.batch_change_entity_parts(get_aux_meta().exposed_boundary_part(),
      mMesh.mesh_meta_data().side_rank(),
      {&get_aux_meta().active_part()},
      {});

  MeshSurface surface(mMesh.mesh_meta_data(),
      coordinates_field(mMesh.mesh_meta_data()),
      get_aux_meta().exposed_boundary_part(),
      1);
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_EQ(4u * 4u, surface.num_owned_facets());
  EXPECT_EQ(10u, surface.global_num_vertices());
}

TEST_F(StkMeshTri6Fixture, meshSurfaceOfQuadraticSides_splitAtMidNodes)
{
  if (stk::parallel_machine_size(get_comm()) > 1) return;
  build_mesh_from_spec(RegularTri6());
  mMesh.batch_change_entity_parts(get_aux_meta().exposed_boundary_part(),
      mMesh.mesh_meta_data().side_rank(),
      {&get_aux_meta().active_part()},
      {});

  MeshSurface surface(mMesh.mesh_meta_data(),
      coordinates_field(mMesh.mesh_meta_data()),
      get_aux_meta().exposed_boundary_part(),
      1);
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_EQ(3u * 2u, surface.num_owned_facets());
  EXPECT_EQ(6u, surface.global_num_vertices());
}

// Sides not in the active part are left out.
TEST(MeshSurface, inactiveSides_leftOut)
{
  BoundingBoxMeshTet4 tets;
  populate_bounding_box_mesh_and_activate(tets, {0., 0., 0.}, {1., 1., 1.}, 0.5);
  const stk::mesh::MetaData & meta = tets.meta_data();
  stk::mesh::BulkData & mesh = tets.bulk_data();
  stk::mesh::Part & xMinus = *meta.get_part("Xminus");
  const stk::mesh::Selector bothX =
      stk::mesh::selectUnion(stk::mesh::PartVector{&xMinus, meta.get_part("Xplus")});

  MeshSurface surface(meta, coordinates_field(meta), bothX, 1);
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  const size_t numBothX = global_sum(surface.num_owned_facets());

  const size_t numXminus = global_sum(stk::mesh::count_selected_entities(
      xMinus & meta.locally_owned_part(), mesh.buckets(meta.side_rank())));
  mesh.batch_change_entity_parts(
      xMinus, meta.side_rank(), {}, {&AuxMetaData::get(meta).active_part()});
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_LT(0u, numXminus);
  EXPECT_EQ(numBothX - numXminus, global_sum(surface.num_owned_facets()));
}

// Each side's facet keeps the side's node order when the side faces out of the blocks touching its
// sideset, or when the selector is not a union (no orientation), and reverses it otherwise; sign =
// -1 reverses each.  Mapped to block_2, which has no elements, a side's orientation comes from the
// unselected block_1 element it touches, reversed (as for a side seen only from the other side
// without an aura).
TEST(MeshSurface, twoTetSides_orientedByTouchingBlockOrFromTheOtherSideOrNotAtAll)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;

  for (const bool mapSurfaceToEmptyBlock : {false, true})
  {
    TwoTetMeshWithBoundarySideset twoTets({.mapSurfaceToEmptyBlock = mapSurfaceToEmptyBlock});
    const stk::mesh::BulkData & mesh = *twoTets.myMesh;
    std::vector<stk::mesh::Entity> sides;
    stk::mesh::get_selected_entities(*twoTets.mySurface & twoTets.myMeta->locally_owned_part(),
        mesh.buckets(twoTets.myMeta->side_rank()),
        sides);
    for (const bool unionSelector : {true, false})
    {
      const stk::mesh::Selector selector = unionSelector ?
          stk::mesh::Selector(*twoTets.mySurface) :
          (*twoTets.mySurface & twoTets.myMeta->universal_part());
      for (const int sign : {1, -1})
      {
        const int expected = sign * ((unionSelector && mapSurfaceToEmptyBlock) ? -1 : 1);
        MeshSurface surface(*twoTets.myMeta, *twoTets.myCoords, selector, sign);
        surface.prepare_to_compute(0., sidesetQueryBox, 0.);
        const FacetMesh & facets = surface.facet_mesh();
        EXPECT_EQ(sides.size(), facets.num_owned_facets());
        for (unsigned f = 0; f < facets.num_owned_facets(); ++f)
        {
          EXPECT_EQ(expected, orientation_relative_to_side(mesh, sides, facets, f))
              << "empty block " << mapSurfaceToEmptyBlock << ", union " << unionSelector
              << ", sign " << sign;
        }
      }
    }
  }
}

// Every facet of the split quadrilaterals faces out of the block: just inside each sub-triangle of
// the 3 non-planar quads at the moved corner, and of a planar quad, the sign is that of the centre,
// and just outside it is the opposite.  sign = -1 negates them all.
TEST(MeshSurface, hexWithMovedCorner_quadFacetsFaceOutOfTheBlock)
{
  BoundingBoxMeshHex8 hex;
  populate_bounding_box_mesh_and_activate(hex, {0., 0., 0.}, {1., 1., 1.}, 0.5);
  stk::mesh::Field<double> & coords = coordinates_field(hex.meta_data());
  const stk::math::Vector3d corner(0., 0., -0.5);
  move_nodes_at(hex.bulk_data(), coords, stk::math::Vector3d(0., 0., 0.), corner);

  using Quad = std::array<stk::math::Vector3d, 4>;
  const std::vector<std::pair<Quad, stk::math::Vector3d>> cornerQuadsAndOutwardAxes{
      {{{corner, {0.5, 0., 0.}, {0.5, 0.5, 0.}, {0., 0.5, 0.}}}, {0., 0., -1.}},
      {{{corner, {0., 0.5, 0.}, {0., 0.5, 0.5}, {0., 0., 0.5}}}, {-1., 0., 0.}},
      {{{corner, {0.5, 0., 0.}, {0.5, 0., 0.5}, {0., 0., 0.5}}}, {0., -1., 0.}}};
  std::vector<std::pair<stk::math::Vector3d, stk::math::Vector3d>> pointsAndOutwardNormals;
  for (auto & [quad, axis] : cornerQuadsAndOutwardAxes)
  {
    const stk::math::Vector3d centroid = 0.25 * (quad[0] + quad[1] + quad[2] + quad[3]);
    for (unsigned i = 0; i < 4; ++i)
      pointsAndOutwardNormals.emplace_back((quad[i] + quad[(i + 1) % 4] + centroid) / 3., axis);
  }
  // Off both diagonals of a planar top quad, so in each of its triangles whichever way it is split.
  pointsAndOutwardNormals.emplace_back(
      stk::math::Vector3d(0.6, 0.8, 1.), stk::math::Vector3d(0., 0., 1.));
  pointsAndOutwardNormals.emplace_back(
      stk::math::Vector3d(0.8, 0.6, 1.), stk::math::Vector3d(0., 0., 1.));

  stk::mesh::PartVector walls;
  for (auto name : {"Xminus", "Xplus", "Yminus", "Yplus", "Zminus", "Zplus"})
    walls.push_back(hex.meta_data().get_part(name));
  const stk::math::Vector3d centre(0.5, 0.5, 0.5);
  std::array<int, 2> centreSigns{0, 0};
  for (const int sign : {1, -1})
  {
    MeshSurface surface(hex.meta_data(), coords, stk::mesh::selectUnion(walls), sign);
    surface.prepare_to_compute(0., sidesetQueryBox, 0.);
    expect_inside_and_outside_signs(surface, centre, pointsAndOutwardNormals);
    centreSigns[sign == 1 ? 0 : 1] = sign_of(surface.point_signed_distance(centre));
  }
  EXPECT_EQ(-centreSigns[0], centreSigns[1]);
}

// Every facet of the split quadratic sides faces out of the element: just inside the centroid of
// each of the 4 triangles of each tet10 face, the sign is that of the element's centroid, and just
// outside it is the opposite; sign = -1 negates them.
TEST_F(StkMeshTet10Fixture, meshSurfaceOfQuadraticSides_facetsFaceOutOfTheElement)
{
  if (stk::parallel_machine_size(get_comm()) > 1) return;
  build_mesh_from_spec(RegularTet10());
  mMesh.batch_change_entity_parts(get_aux_meta().exposed_boundary_part(),
      mMesh.mesh_meta_data().side_rank(),
      {&get_aux_meta().active_part()},
      {});

  const std::vector<stk::math::Vector3d> v = RegularTet().nodeLocs;
  const stk::math::Vector3d centroid = 0.25 * (v[0] + v[1] + v[2] + v[3]);
  std::vector<std::pair<stk::math::Vector3d, stk::math::Vector3d>> pointsAndOutwardNormals;
  const std::array<std::array<unsigned, 3>, 4> faces{{{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}};
  for (const std::array<unsigned, 3> & face : faces)
  {
    const stk::math::Vector3d & a = v[face[0]];
    const stk::math::Vector3d & b = v[face[1]];
    const stk::math::Vector3d & c = v[face[2]];
    const stk::math::Vector3d ab = 0.5 * (a + b), bc = 0.5 * (b + c), ca = 0.5 * (c + a);
    const stk::math::Vector3d normal = outward_unit_normal(a, b, c, centroid);
    for (auto & tri : {std::array{a, ab, ca},
             std::array{b, bc, ab},
             std::array{c, ca, bc},
             std::array{ab, bc, ca}})
      pointsAndOutwardNormals.emplace_back((tri[0] + tri[1] + tri[2]) / 3., normal);
  }

  std::array<int, 2> centroidSigns{0, 0};
  for (const int sign : {1, -1})
  {
    MeshSurface surface(mMesh.mesh_meta_data(),
        coordinates_field(mMesh.mesh_meta_data()),
        get_aux_meta().exposed_boundary_part(),
        sign);
    surface.prepare_to_compute(0., sidesetQueryBox, 0.);
    expect_inside_and_outside_signs(surface, centroid, pointsAndOutwardNormals);
    centroidSigns[sign == 1 ? 0 : 1] = sign_of(surface.point_signed_distance(centroid));
  }
  EXPECT_EQ(-centroidSigns[0], centroidSigns[1]);
}

// As for the tet10, with the 2 segments of each tri6 edge, probed at the edge's quarter points.
TEST_F(StkMeshTri6Fixture, meshSurfaceOfQuadraticSides_facetsFaceOutOfTheElement)
{
  if (stk::parallel_machine_size(get_comm()) > 1) return;
  build_mesh_from_spec(RegularTri6());
  mMesh.batch_change_entity_parts(get_aux_meta().exposed_boundary_part(),
      mMesh.mesh_meta_data().side_rank(),
      {&get_aux_meta().active_part()},
      {});

  std::vector<stk::math::Vector3d> v;
  for (auto & x : RegularTri().nodeLocs)
    v.emplace_back(x[0], x[1], 0.);
  const stk::math::Vector3d centroid = (v[0] + v[1] + v[2]) / 3.;
  std::vector<std::pair<stk::math::Vector3d, stk::math::Vector3d>> pointsAndOutwardNormals;
  for (unsigned i = 0; i < 3; ++i)
  {
    const stk::math::Vector3d & a = v[i];
    const stk::math::Vector3d & b = v[(i + 1) % 3];
    stk::math::Vector3d normal = stk::math::Vector3d(b[1] - a[1], a[0] - b[0], 0.).unit_vector();
    if (Dot(normal, 0.5 * (a + b) - centroid) < 0.) normal = -normal;
    pointsAndOutwardNormals.emplace_back(0.75 * a + 0.25 * b, normal);
    pointsAndOutwardNormals.emplace_back(0.25 * a + 0.75 * b, normal);
  }

  std::array<int, 2> centroidSigns{0, 0};
  for (const int sign : {1, -1})
  {
    MeshSurface surface(mMesh.mesh_meta_data(),
        coordinates_field(mMesh.mesh_meta_data()),
        get_aux_meta().exposed_boundary_part(),
        sign);
    surface.prepare_to_compute(0., sidesetQueryBox, 0.);
    expect_inside_and_outside_signs(surface, centroid, pointsAndOutwardNormals);
    centroidSigns[sign == 1 ? 0 : 1] = sign_of(surface.point_signed_distance(centroid));
  }
  EXPECT_EQ(-centroidSigns[0], centroidSigns[1]);
}

// With the sides in two sidesets, each is two coincident facets, removed to one by
// set_remove_coincident_facets.  set_check_watertight warns (on processor 0) only once a side is
// taken out of both sidesets, after which the rebuilt surface is open.
TEST(MeshSurface, sidesInTwoSidesets_coincidentFacetsRemovedAndOpenSurfaceWarned)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;

  TwoTetMeshWithBoundarySideset twoTets({.sidesAlsoInSecondSurface = true});
  const stk::mesh::Selector bothSidesets =
      stk::mesh::selectUnion(stk::mesh::PartVector{twoTets.mySurface, twoTets.mySecondSurface});
  {
    MeshSurface surface(*twoTets.myMeta, *twoTets.myCoords, bothSidesets, 1);
    surface.prepare_to_compute(0., sidesetQueryBox, 0.);
    EXPECT_EQ(12u, global_sum(surface.num_owned_facets()));
  }

  MeshSurface surface(*twoTets.myMeta, *twoTets.myCoords, bothSidesets, 1);
  surface.set_remove_coincident_facets(true);
  surface.set_check_watertight(true);
  const bool isProc0 = 0 == stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned numWarningsBefore = stk::get_warning_count();
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_EQ(6u, global_sum(surface.num_owned_facets()));
  if (isProc0)
  {
    EXPECT_EQ(numWarningsBefore, stk::get_warning_count());
  }

  std::vector<stk::mesh::Entity> sidesToRemove;
  stk::mesh::get_selected_entities(*twoTets.mySurface & twoTets.myMeta->locally_owned_part(),
      twoTets.myMesh->buckets(twoTets.myMeta->side_rank()),
      sidesToRemove);
  sidesToRemove.resize((isProc0 && !sidesToRemove.empty()) ? 1 : 0);
  twoTets.myMesh->batch_change_entity_parts(sidesToRemove,
      stk::mesh::PartVector{},
      stk::mesh::PartVector{twoTets.mySurface, twoTets.mySecondSurface});
  surface.prepare_to_compute(0., sidesetQueryBox, 0.);
  EXPECT_EQ(5u, global_sum(surface.num_owned_facets()));
  if (isProc0)
  {
    EXPECT_LT(numWarningsBefore, stk::get_warning_count());
  }
}

} // namespace krino
