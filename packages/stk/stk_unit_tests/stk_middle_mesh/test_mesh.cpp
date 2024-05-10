#include "gtest/gtest.h"

#include <cmath>
#include <fstream>
#include <set>

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "util/meshes.hpp"
#include "stk_middle_mesh/utils.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

namespace {
enum class GeoClassificationE
{
  Vertex = 0,
  Edge,
  Face
};

GeoClassificationE get_geo_classification(const MeshSpec& spec, const utils::Point& pt)
{
  double eps  = 1e-12;
  auto x      = pt.get_x();
  auto y      = pt.get_y();
  bool isXmin = std::abs(x - spec.xmin) < eps;
  bool isXmax = std::abs(x - spec.xmax) < eps;
  bool isYmin = std::abs(y - spec.ymin) < eps;
  bool isYmax = std::abs(y - spec.ymax) < eps;

  if ((isXmin && isYmin) || (isXmin && isYmax) || (isXmax && isYmin) || (isXmax && isYmax))
    return GeoClassificationE::Vertex;
  else if (isXmin || isXmax || isYmin || isYmax)
    return GeoClassificationE::Edge;
  else
    return GeoClassificationE::Face;
}

GeoClassificationE get_edge_geo_classification(const MeshSpec& spec, MeshEntityPtr edge)
{
  MeshEntityPtr verts[MAX_DOWN];
  get_downward(edge, 0, verts);
  GeoClassificationE g1 = get_geo_classification(spec, verts[0]->get_point_orig(0));
  GeoClassificationE g2 = get_geo_classification(spec, verts[1]->get_point_orig(0));

  return static_cast<GeoClassificationE>(std::max(static_cast<int>(g1), static_cast<int>(g2)));
}

MeshEntityPtr get_closest_entity(std::shared_ptr<Mesh> mesh, int dim, const utils::Point& pt, double tol=1e-13)
{
  MeshEntityPtr closestEntity = nullptr;
  double closestDistance    = std::numeric_limits<double>::max();
  for (auto entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      auto ptVert        = mesh::compute_centroid(entity);
      double distSquared = dot(pt - ptVert, pt - ptVert);
      if (distSquared < closestDistance && std::sqrt(distSquared) < tol)
      {
        closestDistance = distSquared;
        closestEntity   = entity;
      }
    }

  return closestEntity;
}  

MeshEntityPtr get_closest_vert(std::shared_ptr<Mesh> mesh, const utils::Point& pt)
{
  return get_closest_entity(mesh, 0, pt, 1);
}


void expect_near(const utils::Point& pt1, const utils::Point& pt2, double tol)
{
  EXPECT_NEAR(pt1.x, pt2.x, tol);
  EXPECT_NEAR(pt1.y, pt2.y, tol);
  EXPECT_NEAR(pt1.z, pt2.z, tol);
}

} // namespace

TEST(Mesh, TopologyDown)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto dx = (spec.xmax - spec.xmin) / spec.numelX;
  auto dy = (spec.ymax - spec.ymin) / spec.numelY;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  // verify counts of mesh entities
  EXPECT_EQ(count_valid(mesh->get_vertices()), (spec.numelX + 1) * (spec.numelY + 1));
  EXPECT_EQ(count_valid(mesh->get_edges()), (spec.numelX + 1) * spec.numelY + (spec.numelY + 1) * spec.numelX);
  EXPECT_EQ(count_valid(mesh->get_elements()), spec.numelX * spec.numelY);

  MeshEntityPtr edges[MAX_DOWN], verts1[MAX_DOWN], verts2[MAX_DOWN];
  auto& elements = mesh->get_elements();
  for (auto el : elements)
  {
    if (!el)
      continue;

    int nedges = get_downward(el, 1, edges);
    EXPECT_EQ(nedges, 4);
    std::set<MeshEntityPtr> vertsFromEdges;
    for (int i = 0; i < nedges; ++i)
    {
      int nverts1 = get_downward(edges[i], 0, verts1);
      EXPECT_EQ(nverts1, 2);
      vertsFromEdges.insert(verts1[0]);
      vertsFromEdges.insert(verts1[1]);
    }

    // check the same vertices from the edges and the vertices
    int nverts2 = get_downward(el, 0, verts2);
    EXPECT_EQ(nverts2, 4);
    for (int i = 0; i < nverts2; ++i)
      EXPECT_EQ(vertsFromEdges.count(verts2[i]), static_cast<unsigned int>(1));
    EXPECT_EQ(vertsFromEdges.size(), static_cast<unsigned int>(4));

    // check the vertices are oriented counterclockwise
    auto p0 = verts2[0]->get_point_orig(0);
    auto p1 = verts2[1]->get_point_orig(0);
    auto p2 = verts2[2]->get_point_orig(0);
    auto p3 = verts2[3]->get_point_orig(0);

    EXPECT_FLOAT_EQ(p0.get_x() + dx, p1.get_x());
    EXPECT_FLOAT_EQ(p0.get_y(), p1.get_y());

    EXPECT_FLOAT_EQ(p0.get_x() + dx, p2.get_x());
    EXPECT_FLOAT_EQ(p0.get_y() + dy, p2.get_y());

    EXPECT_FLOAT_EQ(p0.get_x(), p3.get_x());
    EXPECT_FLOAT_EQ(p0.get_y() + dy, p3.get_y());
  }

  // check edges have correct size
  for (auto& edge : mesh->get_edges())
  {
    get_downward(edge, 0, verts1);
    utils::Point pt1     = verts1[0]->get_point_orig(0);
    utils::Point pt2     = verts1[1]->get_point_orig(0);
    double len           = std::sqrt(std::pow(pt2.get_x() - pt1.get_x(), 2) + std::pow(pt2.get_y() - pt1.get_y(), 2));
    bool isCorrectLength = std::abs(len - dx) < 1e-12 || std::abs(len - dy) < 1e-12;
    EXPECT_TRUE(isCorrectLength);
  }
}

TEST(Mesh, TopologyUp)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  auto& verts = mesh->get_vertices();
  for (auto vert : verts)
  {
    if (!vert)
      continue;

    auto g     = get_geo_classification(spec, vert->get_point_orig(0));
    int nedges = vert->count_up();
    if (g == GeoClassificationE::Vertex)
      EXPECT_EQ(nedges, 2);
    else if (g == GeoClassificationE::Edge)
      EXPECT_EQ(nedges, 3);
    else // Face
      EXPECT_EQ(nedges, 4);

    std::vector<MeshEntityPtr> edges;
    for (int i = 0; i < nedges; ++i)
      edges.push_back(vert->get_up(i));

    EXPECT_TRUE(is_unique(edges));
  }

  auto& edges = mesh->get_edges();
  for (auto edge : edges)
  {
    if (!edge)
      continue;

    auto g     = get_edge_geo_classification(spec, edge);
    int nfaces = edge->count_up();
    if (g == GeoClassificationE::Edge)
      EXPECT_EQ(nfaces, 1);
    else if (g == GeoClassificationE::Face)
      EXPECT_EQ(nfaces, 2);
    else // Vertex
      EXPECT_TRUE(false);

    std::vector<MeshEntityPtr> faces;
    for (int i = 0; i < nfaces; ++i)
      faces.push_back(edge->get_up(i));

    EXPECT_TRUE(is_unique(faces));
  }
}

int count(MeshEntityPtr vals[], const int n, MeshEntityPtr val)
{
  int nfound = 0;
  for (int i = 0; i < n; ++i)
    if (vals[i] == val)
      nfound += 1;

  return nfound;
}

void test_upward(std::shared_ptr<Mesh> mesh, MeshEntityPtr eIn, const int dimOut, std::vector<MeshEntityPtr>& entities)
{
  int dimIn = get_type_dimension(eIn->get_type());
  // check dimension
  for (auto& e : entities)
    EXPECT_EQ(get_type_dimension(e->get_type()), dimOut);

  // test all returned entities have e_in as a downward adjacency,
  // and no other entities in the mesh do
  MeshEntityPtr down[MAX_DOWN];
  for (auto& e : mesh->get_mesh_entities(dimOut))
  {
    int ndown  = get_downward(e, dimIn, down);
    int nfound = count(down, ndown, eIn);

    if (std::find(entities.begin(), entities.end(), e) == entities.end())
      EXPECT_EQ(nfound, 0);
    else
      EXPECT_EQ(nfound, 1);
  }
}

TEST(Mesh, TopologyUpward)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  std::vector<MeshEntityPtr> entities;
  for (auto& v : mesh->get_vertices())
    for (int dim = 1; dim < 3; ++dim)
    {
      get_upward(v, dim, entities);
      test_upward(mesh, v, dim, entities);
    }

  for (auto& edge : mesh->get_edges())
  {
    get_upward(edge, 2, entities);
    test_upward(mesh, edge, 2, entities);
  }
}

void test_bridge(std::shared_ptr<Mesh> mesh, MeshEntityPtr eIn, const int viaDim, const int targetDim,
                 std::vector<MeshEntityPtr>& entities)
{
  // check dimension
  for (auto& e : entities)
    EXPECT_EQ(get_type_dimension(e->get_type()), targetDim);

  // verify uniqueness
  auto it = std::unique(entities.begin(), entities.end());
  EXPECT_EQ(it, entities.end());

  std::vector<MeshEntityPtr> viaEntities(MAX_DOWN);
  int eDim = get_type_dimension(eIn->get_type());

  int nVia;
  if (viaDim < eDim)
    nVia = get_downward(eIn, viaDim, entities.data());
  else
    nVia = get_upward(eIn, viaDim, entities);

  // test that each entity has at least one of entities_via as an adjacency
  std::vector<MeshEntityPtr> viaEntities2(MAX_DOWN);
  for (auto& e : mesh->get_mesh_entities(targetDim))
  {
    int nVia2;
    if (viaDim < targetDim)
      nVia2 = get_downward(e, viaDim, viaEntities2.data());
    else
      nVia2 = get_upward(e, viaDim, viaEntities2);

    viaEntities2.resize(nVia2);

    int nfound = 0;
    for (auto e2 : viaEntities2)
      nfound += count(viaEntities.data(), nVia, e2);

    // test that each entity in entities has at least one of the via_enties
    // as an adjacency, and that none of the other entities do
    if (std::find(entities.begin(), entities.end(), e) == entities.end())
      EXPECT_EQ(nfound, 0);
    else
      EXPECT_GE(nfound, 0);
  }
}

TEST(Mesh, TopologyBridge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 5;
  spec.numelY = 6;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  std::vector<MeshEntityPtr> entities;

  for (auto& v : mesh->get_vertices())
  {
    get_bridge_adjacent(v, 1, 2, entities);
    test_bridge(mesh, v, 1, 2, entities);

    get_bridge_adjacent(v, 2, 1, entities);
    test_bridge(mesh, v, 2, 1, entities);

    get_bridge_adjacent(v, 1, 0, entities);
    test_bridge(mesh, v, 1, 0, entities);

    get_bridge_adjacent(v, 2, 0, entities);
    test_bridge(mesh, v, 2, 0, entities);
  }

  for (auto& e : mesh->get_edges())
  {
    get_bridge_adjacent(e, 0, 1, entities);
    test_bridge(mesh, e, 0, 1, entities);

    get_bridge_adjacent(e, 0, 2, entities);
    test_bridge(mesh, e, 0, 2, entities);

    get_bridge_adjacent(e, 2, 0, entities);
    test_bridge(mesh, e, 2, 0, entities);

    get_bridge_adjacent(e, 2, 1, entities);
    test_bridge(mesh, e, 2, 1, entities);
  }

  for (auto& e : mesh->get_elements())
  {
    get_bridge_adjacent(e, 0, 1, entities);
    test_bridge(mesh, e, 0, 1, entities);

    get_bridge_adjacent(e, 0, 2, entities);
    test_bridge(mesh, e, 0, 2, entities);

    get_bridge_adjacent(e, 1, 0, entities);
    test_bridge(mesh, e, 1, 0, entities);

    get_bridge_adjacent(e, 1, 2, entities);
    test_bridge(mesh, e, 1, 2, entities);
  }
}

TEST(Mesh, EntityCreationErrors)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  MeshEntityPtr down[MAX_DOWN], down2[MAX_DOWN];
  MeshEntityPtr el = mesh->get_elements()[0];
  get_downward(el, 0, down);
  get_downward(el, 1, down2);

  EXPECT_ANY_THROW(mesh->create_quad_from_verts(down[0], down[1], down[2], down[3]));
  EXPECT_ANY_THROW(mesh->create_quad_from_verts(down[1], down[0], down[3], down[2]));

  EXPECT_ANY_THROW(mesh->create_quad_from_verts(down[0], down[0], down[2], down[3]));
  EXPECT_ANY_THROW(mesh->create_edge(down[0], down[1]));
  EXPECT_ANY_THROW(mesh->create_edge(down[1], down[0]));
  EXPECT_ANY_THROW(mesh->create_edge(down[0], down[0]));

  EXPECT_ANY_THROW(mesh->create_quad(down2[0], down2[1], down2[2], down2[3]));
  EXPECT_ANY_THROW(mesh->create_quad(down2[0], down2[0], down2[2], down2[3]));
}

utils::Point apply_rotation(const double theta, const utils::Point& pt)
{
  double x2 = std::cos(theta) * pt.get_x() - std::sin(theta) * pt.get_y();
  double y2 = std::sin(theta) * pt.get_x() + std::cos(theta) * pt.get_y();

  return utils::Point(x2, y2);
}

TEST(Mesh, XiCoords)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto mesh = make_empty_mesh();

  auto v1  = mesh->create_vertex(0, 0);
  auto v2  = mesh->create_vertex(2, 0);
  auto v3  = mesh->create_vertex(2, 2);
  auto v4  = mesh->create_vertex(0, 2);
  auto el  = mesh->create_quad_from_verts(v1, v2, v3, v4);
  auto tri = mesh->create_triangle_from_verts(v1, v2, v3);

  MeshEntityPtr edges[MAX_DOWN];
  get_downward(el, 1, edges);
  utils::Point pt;

  // compute_quad_coords_from_xi
  pt = compute_quad_coords_from_xi(el, utils::Point(0, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = compute_quad_coords_from_xi(el, utils::Point(1, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = compute_quad_coords_from_xi(el, utils::Point(1, 1));
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 2);

  pt = compute_quad_coords_from_xi(el, utils::Point(0, 1));
  EXPECT_FLOAT_EQ(pt.get_x(), 0);
  EXPECT_FLOAT_EQ(pt.get_y(), 2);

  pt = compute_quad_coords_from_xi(el, utils::Point(0.5, 0.5));
  EXPECT_FLOAT_EQ(pt.get_x(), 1);
  EXPECT_FLOAT_EQ(pt.get_y(), 1);

  // compute_tri_coords_from_xi
  pt = compute_coords_from_xi(tri, utils::Point(1.0 / 3.0, 1.0 / 3.0));
  EXPECT_FLOAT_EQ(pt.get_x(), 4.0 / 3.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 2.0 / 3.0);

  pt = compute_tri_coords_from_xi_3d(tri, utils::Point(1.0 / 3.0, 1.0 / 3.0));
  EXPECT_FLOAT_EQ(pt.get_x(), 4.0 / 3.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 2.0 / 3.0);

  pt = compute_coords_from_xi(tri, utils::Point(0, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = compute_tri_coords_from_xi_3d(tri, utils::Point(0, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 0);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = compute_coords_from_xi(tri, utils::Point(1, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = compute_tri_coords_from_xi_3d(tri, utils::Point(1, 0));
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 0);

  pt = compute_coords_from_xi(tri, utils::Point(0, 1));
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 2);

  pt = compute_tri_coords_from_xi_3d(tri, utils::Point(0, 1));
  EXPECT_FLOAT_EQ(pt.get_x(), 2);
  EXPECT_FLOAT_EQ(pt.get_y(), 2);

  // Quad centroid
  pt = compute_centroid(el);
  EXPECT_FLOAT_EQ(pt.get_x(), 1);
  EXPECT_FLOAT_EQ(pt.get_y(), 1);

  // Triangle centroid
  pt = compute_centroid(tri);
  EXPECT_FLOAT_EQ(pt.get_x(), 4.0 / 3.0);
  EXPECT_FLOAT_EQ(pt.get_y(), 2.0 / 3.0);
}

TEST(Mesh, CoordsFromXi3D)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto mesh = make_empty_mesh();

  auto v1  = mesh->create_vertex(0, 0, 0);
  auto v2  = mesh->create_vertex(2, 0, 0);
  auto v3  = mesh->create_vertex(2, 0, 2);
  auto v4  = mesh->create_vertex(0, 0, 2);
  auto el  = mesh->create_quad_from_verts(v1, v2, v3, v4);
  auto tri = mesh->create_triangle_from_verts(v1, v2, v4);

  MeshEntityPtr edges[MAX_DOWN];
  get_downward(el, 1, edges);
  utils::Point pt;

  pt = compute_coords_from_xi_3d(el, {0.25, 0.25});
  EXPECT_NEAR(pt.x, 0.5, 1e-13);
  EXPECT_NEAR(pt.y, 0, 1e-13);
  EXPECT_NEAR(pt.z, 0.5, 1e-13);

  pt = compute_coords_from_xi_3d(tri, {0.25, 0.25});
  EXPECT_NEAR(pt.x, 0.5, 1e-13);
  EXPECT_NEAR(pt.y, 0, 1e-13);
  EXPECT_NEAR(pt.z, 0.5, 1e-13);
}

TEST(Mesh, Centroid3D)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto mesh = make_empty_mesh();

  auto v1  = mesh->create_vertex(0, 0, 0);
  auto v2  = mesh->create_vertex(2, 0, 0);
  auto v3  = mesh->create_vertex(2, 0, 2);
  auto v4  = mesh->create_vertex(0, 0, 2);
  auto el  = mesh->create_quad_from_verts(v1, v2, v3, v4);
  auto tri = mesh->create_triangle_from_verts(v1, v2, v3);

  expect_near(compute_quad_centroid_3d(el), utils::Point(1, 0, 1), 1e-13);
  expect_near(compute_tri_centroid_3d(tri), utils::Point(4.0 / 3.0, 0, 2.0 / 3.0), 1e-13);
}

TEST(Mesh, ErrorChecking)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 10;
  spec.numelY = 10;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh1 = create_mesh(spec, func);

  EXPECT_EQ(check_angles(mesh1, 2, 178), 0);
  EXPECT_NO_THROW(check_topology(mesh1));
  EXPECT_NO_THROW(check_coordinate_field(mesh1));
}

TEST(Mesh, SetDownOrientation)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 1;
  spec.numelY = 1;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;  
  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  mesh::MeshEntityPtr el1 = get_closest_entity(mesh, 2, {0.5, 0.5});

  el1->set_down_orientation(0, mesh::EntityOrientation::Standard);
  EXPECT_EQ(el1->get_down_orientation(0), mesh::EntityOrientation::Standard);

  el1->set_down_orientation(0, mesh::EntityOrientation::Reversed);
  EXPECT_EQ(el1->get_down_orientation(0), mesh::EntityOrientation::Reversed);
}

TEST(Mesh, ReplaceDown)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 1;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;  
  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  mesh::MeshEntityPtr el1 = get_closest_entity(mesh, 2, {0.5, 0.25});

  mesh::MeshEntityPtr v1 = mesh->create_vertex(0,   0.5, 0);
  mesh::MeshEntityPtr v2 = mesh->create_vertex(0.5, 0.5, 0);

  mesh::MeshEntityPtr newEdge = mesh->create_edge(v1, v2);

  el1->replace_down(2, newEdge, mesh::EntityOrientation::Reversed);
  EXPECT_EQ(el1->get_down(2), newEdge);
  EXPECT_EQ(el1->get_down_orientation(2), mesh::EntityOrientation::Reversed);
}

TEST(Mesh, EntityOrientationReverse)
{
  EXPECT_EQ(reverse(mesh::EntityOrientation::Standard), mesh::EntityOrientation::Reversed);
  EXPECT_EQ(reverse(mesh::EntityOrientation::Reversed), mesh::EntityOrientation::Standard);
}

TEST(Mesh, ReverseEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 1;
  spec.numelY = 2;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;  
  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  mesh::MeshEntityPtr edge = get_closest_entity(mesh, 1, {0.5, 0.5});
  mesh::MeshEntityPtr el1  = get_closest_entity(mesh, 2,  {0.5, 0.25});
  mesh::MeshEntityPtr el2  = get_closest_entity(mesh, 2,  {0.5, 0.75});
  mesh::MeshEntityPtr v0   = edge->get_down(0);
  mesh::MeshEntityPtr v1   = edge->get_down(1);
  mesh::EntityOrientation orient1 = el1->get_down_orientation(2);
  mesh::EntityOrientation orient2 = el2->get_down_orientation(0);

  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts1, verts2;
  mesh::get_downward(el1, 0, verts1.data());
  mesh::get_downward(el2, 0, verts2.data());

  mesh::reverse_edge(edge);

  EXPECT_EQ(edge->get_down(0), v1);
  EXPECT_EQ(edge->get_down(1), v0);
  EXPECT_EQ(el1->get_down(2), edge);
  EXPECT_EQ(el2->get_down(0), edge);
  EXPECT_EQ(el1->get_down_orientation(2), reverse(orient1));
  EXPECT_EQ(el2->get_down_orientation(0), reverse(orient2));


  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts1After, verts2After;
  mesh::get_downward(el1, 0, verts1After.data());
  mesh::get_downward(el2, 0, verts2After.data());

  for (int i=0; i < 4; ++i)
  {
    EXPECT_EQ(verts1[i], verts1After[i]);
    EXPECT_EQ(verts2[i], verts2After[i]);
  }
}

TEST(Mesh, DeleteFace)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  int nverts = 16;
  int nedges = 24;
  int nelem  = 9;

  auto func = [&](const utils::Point& pt) { return pt; };

  // delete corner element
  {
    std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

    // find element
    std::vector<MeshEntityPtr> elements;
    for (auto& vert : mesh->get_vertices())
      if (vert && get_upward(vert, 2, elements) == 1)
        break;

    auto el = elements.front();
    mesh->delete_face(el);
    EXPECT_EQ(count_valid(mesh->get_vertices()), nverts - 1);
    EXPECT_EQ(count_valid(mesh->get_edges()), nedges - 2);
    EXPECT_EQ(count_valid(mesh->get_elements()), nelem - 1);
    EXPECT_NO_THROW(check_topology(mesh));
  }

  // delete element on boundary of domain but not in corner
  {
    std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

    // find element
    std::vector<MeshEntityPtr> elements;
    MeshEntityPtr el = nullptr;
    for (auto& edge : mesh->get_edges())
      if (edge && edge->count_up() == 1)
      {
        auto v1  = edge->get_down(0);
        auto v2  = edge->get_down(1);
        int nup1 = get_upward(v1, 2, elements);
        int nup2 = get_upward(v2, 2, elements);
        if (nup1 > 1 && nup2 > 1)
        {
          el = edge->get_up(0);
          break;
        }
      }

    assert(el);

    mesh->delete_face(el);
    EXPECT_EQ(count_valid(mesh->get_vertices()), nverts);
    EXPECT_EQ(count_valid(mesh->get_edges()), nedges - 1);
    EXPECT_EQ(count_valid(mesh->get_elements()), nelem - 1);
    EXPECT_NO_THROW(check_topology(mesh));
  }

  // delete interior element
  {
    std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

    // find element
    std::vector<MeshEntityPtr> elements;
    MeshEntityPtr el = nullptr;
    for (auto& elI : mesh->get_elements())
      if (elI)
      {
        int minUp = std::numeric_limits<int>::max();
        for (int i = 0; i < elI->count_down(); ++i)
          minUp = std::min(minUp, elI->get_down(i)->count_up());

        if (minUp > 1)
        {
          el = elI;
          break;
        }
      }

    assert(el);

    mesh->delete_face(el);
    EXPECT_EQ(count_valid(mesh->get_vertices()), nverts);
    EXPECT_EQ(count_valid(mesh->get_edges()), nedges);
    EXPECT_EQ(count_valid(mesh->get_elements()), nelem - 1);
    EXPECT_NO_THROW(check_topology(mesh));
  }
}

TEST(Mesh, SplitEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 3;
  spec.numelY = 3;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  int nverts = 16;
  int nedges = 24 + 9;
  int nelem  = 18;

  auto func = [&](const utils::Point& pt) { return pt; };

  {
    for (int i = 0; i < nedges; ++i)
    {
      std::shared_ptr<Mesh> mesh = create_mesh(spec, func, MPI_COMM_WORLD, true);
      auto edge                  = mesh->get_edges()[i];
      int nUp                    = edge->count_up();
      if (edge)
      {
        mesh->split_edge(edge, 0.5);
        EXPECT_EQ(count_valid(mesh->get_vertices()), nverts + 1);
        if (nUp == 2)
        {
          EXPECT_EQ(count_valid(mesh->get_edges()), nedges + 3);
          EXPECT_EQ(count_valid(mesh->get_elements()), nelem + 2);
        } else // n_up == 1
        {
          EXPECT_EQ(count_valid(mesh->get_edges()), nedges + 2);
          EXPECT_EQ(count_valid(mesh->get_elements()), nelem + 1);
        }

        EXPECT_NO_THROW(check_topology(mesh));
        check_topology(mesh);
      }
    }
  }
}

TEST(Mesh, RemoteSharedEntityAccessors)
{
  auto mesh = make_empty_mesh();
  auto vert = mesh->create_vertex(0, 0, 0);

  RemoteSharedEntity remote{1, 42};

  EXPECT_EQ(vert->count_remote_shared_entities(), 0);
  vert->add_remote_shared_entity(remote);
  EXPECT_EQ(vert->count_remote_shared_entities(), 1);

  auto remoteRetrieved = vert->get_remote_shared_entity(0);
  EXPECT_EQ(remote.remoteRank, remoteRetrieved.remoteRank);
  EXPECT_EQ(remote.remoteId, remoteRetrieved.remoteId);
}

TEST(Mesh, RemoteSharedEntityDeletion)
{
  auto mesh = make_empty_mesh();
  auto vert = mesh->create_vertex(0, 0, 0);

  RemoteSharedEntity remote1{1, 42}, remote2{2, 7}, remote3{3, 8};
  vert->add_remote_shared_entity(remote1);
  vert->add_remote_shared_entity(remote2);
  vert->add_remote_shared_entity(remote3);

  vert->delete_remote_shared_entity(1);

  auto remoteRetrieved1 = vert->get_remote_shared_entity(0);
  auto remoteRetrieved3 = vert->get_remote_shared_entity(1);

  EXPECT_EQ(remote1.remoteRank, remoteRetrieved1.remoteRank);
  EXPECT_EQ(remote1.remoteId, remoteRetrieved1.remoteId);

  EXPECT_EQ(remote3.remoteRank, remoteRetrieved3.remoteRank);
  EXPECT_EQ(remote3.remoteId, remoteRetrieved3.remoteId);
}

TEST(Mesh, ConvertXiToRange)
{
  EXPECT_NEAR(convert_xi_coords_to_range(-1, 1, 0.75),  0.5, 1e-13);
  EXPECT_NEAR(convert_xi_coords_to_range(-1, 1, 0.25), -0.5, 1e-13);
  EXPECT_NEAR(convert_xi_coords_to_range(-1, 1, 0),    -1.0, 1e-13);

  expect_near(convert_xi_coords_to_range(-1, 1, {0.75, 0.25}), utils::Point(0.5, -0.5), 1e-13);
}

TEST(Mesh, ConvertXiFromRange)
{
  EXPECT_NEAR(convert_xi_coords_from_range(-1, 1, 0.5),  0.75, 1e-13);
  EXPECT_NEAR(convert_xi_coords_from_range(-1, 1, -0.5), 0.25, 1e-13);
  EXPECT_NEAR(convert_xi_coords_from_range(-1, 1, 0),    0.5,  1e-13);

  expect_near(convert_xi_coords_from_range(-1, 1, {0.5, -0.5}), utils::Point(0.75, 0.25), 1e-13);
}

namespace {

template <typename T>
class SymmetricCommunication
{
  public:
    explicit SymmetricCommunication(MPI_Comm comm, int tag)
      : m_comm(comm)
      , m_tag(tag)
      , m_sendBufs(utils::impl::comm_size(comm))
      , m_recvBufs(utils::impl::comm_size(comm))
      , m_sendReqs(utils::impl::comm_size(comm))
      , m_recvReqs(utils::impl::comm_size(comm))
    {}

    std::vector<T>& get_send_buf(int rank) { return m_sendBufs[rank]; }

    const std::vector<T>& get_recv_buf(int rank) const { return m_recvBufs[rank]; };

    void start_communication()
    {
      std::fill(m_sendReqs.begin(), m_sendReqs.end(), MPI_REQUEST_NULL);
      std::fill(m_recvReqs.begin(), m_recvReqs.end(), MPI_REQUEST_NULL);

      for (size_t i = 0; i < m_recvBufs.size(); ++i)
      {
        if (m_sendBufs[i].size() > 0)
        {
          m_recvBufs[i].resize(m_sendBufs[i].size());
          MPI_Irecv(m_recvBufs[i].data(), m_recvBufs[i].size() * sizeof(T), MPI_BYTE, i, m_tag, m_comm,
                    &(m_recvReqs[i]));
          MPI_Isend(m_sendBufs[i].data(), m_sendBufs[i].size() * sizeof(T), MPI_BYTE, i, m_tag, m_comm,
                    &(m_sendReqs[i]));
        }
      }
    }

    void finish_communication()
    {
      MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(m_recvReqs.size(), m_recvReqs.data(), MPI_STATUSES_IGNORE);
    }

  private:
    MPI_Comm m_comm;
    int m_tag;
    std::vector<std::vector<T>> m_sendBufs;
    std::vector<std::vector<T>> m_recvBufs;
    std::vector<MPI_Request> m_sendReqs;
    std::vector<MPI_Request> m_recvReqs;
};

struct RemoteInfo
{
    int srcRank;
    int srcId;
    int destRank;
    int destId;
    utils::Point pt;
};

std::map<int, std::vector<RemoteInfo>> get_vert_shared_entity_info(std::shared_ptr<Mesh> mesh)
{
  SymmetricCommunication<RemoteInfo> comm(mesh->get_comm(), 110);
  int myrank = utils::impl::comm_rank(mesh->get_comm());

  for (auto vert : mesh->get_vertices())
    if (vert)
      for (int i = 0; i < vert->count_remote_shared_entities(); ++i)
      {
        auto remoteInfo = vert->get_remote_shared_entity(i);
        RemoteInfo info{myrank, vert->get_id(), remoteInfo.remoteRank, remoteInfo.remoteId, vert->get_point_orig(0)};
        comm.get_send_buf(remoteInfo.remoteRank).push_back(info);
      }

  comm.start_communication();
  comm.finish_communication();

  std::map<int, std::vector<RemoteInfo>> receivedInfo;
  for (int rank = 0; rank < utils::impl::comm_size(mesh->get_comm()); ++rank)
    for (auto& info : comm.get_recv_buf(rank))
      receivedInfo[info.destId].push_back(info);

  return receivedInfo;
}

struct RemoteEdgeInfo
{
    int srcRank;
    int srcId;
    int destRank;
    int destId;
    int destVert1Id;
    int destVert2Id;
};

std::map<int, std::vector<RemoteEdgeInfo>> get_edge_shared_entity_info(std::shared_ptr<Mesh> mesh)
{
  SymmetricCommunication<RemoteEdgeInfo> comm(mesh->get_comm(), 110);
  int myrank = utils::impl::comm_rank(mesh->get_comm());

  for (auto& edge : mesh->get_edges())
    if (edge)
    {
      auto v1 = edge->get_down(0);
      auto v2 = edge->get_down(1);
      for (int i = 0; i < edge->count_remote_shared_entities(); ++i)
      {
        RemoteSharedEntity remote = edge->get_remote_shared_entity(i);
        int v1DestId = -1, v2DestId = -1;
        for (int j = 0; j < v1->count_remote_shared_entities(); ++j)
          if (v1->get_remote_shared_entity(j).remoteRank == remote.remoteRank)
          {
            v1DestId = v1->get_remote_shared_entity(j).remoteId;
            break;
          }

        for (int j = 0; j < v2->count_remote_shared_entities(); ++j)
          if (v2->get_remote_shared_entity(j).remoteRank == remote.remoteRank)
          {
            v2DestId = v2->get_remote_shared_entity(j).remoteId;
            break;
          }

        comm.get_send_buf(remote.remoteRank)
            .push_back({myrank, edge->get_id(), remote.remoteRank, remote.remoteId, v1DestId, v2DestId});
      }
    }

  comm.start_communication();
  comm.finish_communication();

  std::map<int, std::vector<RemoteEdgeInfo>> receivedInfo;
  for (int rank = 0; rank < utils::impl::comm_size(mesh->get_comm()); ++rank)
    for (auto& info : comm.get_recv_buf(rank))
      receivedInfo[info.destId].push_back(info);

  return receivedInfo;
}

void test_all_shared_entites_present(MeshEntityPtr vert, std::vector<RemoteInfo> remoteInfo)
{
  std::vector<RemoteSharedEntity> remotesLocal;
  for (int i = 0; i < vert->count_remote_shared_entities(); ++i)
    remotesLocal.push_back(vert->get_remote_shared_entity(i));

  auto compareLocal = [](const RemoteSharedEntity& lhs, const RemoteSharedEntity& rhs) {
    return lhs.remoteRank < rhs.remoteRank;
  };
  std::sort(remotesLocal.begin(), remotesLocal.end(), compareLocal);

  auto compareRemote = [](const RemoteInfo& lhs, const RemoteInfo& rhs) { return lhs.srcRank < rhs.srcRank; };
  std::sort(remoteInfo.begin(), remoteInfo.end(), compareRemote);

  for (int i = 0; i < vert->count_remote_shared_entities(); ++i)
  {
    EXPECT_EQ(remotesLocal[i].remoteId, remoteInfo[i].srcId);
    EXPECT_EQ(remotesLocal[i].remoteRank, remoteInfo[i].srcRank);
  }
}

void test_vert_shared_entities(std::shared_ptr<Mesh> mesh, bool checkCoords = true)
{
  // verifies that:
  //  1. all remote entities have the same coordinates
  //  2. Remotes are symmetric (ie if proc 0 has remote on proc 1, then proc 1 has remote on proc 0)
  auto& verts                                        = mesh->get_vertices();
  std::map<int, std::vector<RemoteInfo>> remoteInfos = get_vert_shared_entity_info(mesh);
  for (auto& p : remoteInfos)
  {
    MeshEntityPtr vert = verts[p.first];

    std::vector<RemoteInfo>& remoteInfoForVert = p.second;
    // EXPECT_EQ(vert->count_remote_shared_entities(), remote_info_for_vert.size());
    test_all_shared_entites_present(vert, remoteInfoForVert);

    if (checkCoords)
    {
      utils::Point vertCoords = vert->get_point_orig(0);
      for (int i = 0; i < vert->count_remote_shared_entities(); ++i)
      {
        utils::Point remoteCoords = remoteInfoForVert[i].pt;
        double dist               = std::sqrt(dot(vertCoords - remoteCoords, vertCoords - remoteCoords));
        EXPECT_NEAR(dist, 0, 1e-12);
      }
    }
  }
}

void test_edge_shared_entities(std::shared_ptr<Mesh> mesh)
{
  std::map<int, std::vector<RemoteEdgeInfo>> receivedInfo = get_edge_shared_entity_info(mesh);

  for (auto& p : receivedInfo)
  {
    int edgeId                                 = p.first;
    std::vector<RemoteEdgeInfo>& receivedEdges = p.second;

    EXPECT_EQ(receivedEdges.size(), 1u);

    MeshEntityPtr edge           = mesh->get_edges()[edgeId];
    RemoteEdgeInfo& receivedEdge = receivedEdges[0];
    int v1Id                     = edge->get_down(0)->get_id();
    int v2Id                     = edge->get_down(1)->get_id();

    // Edges are not required to have same orientation, need to figure it out here
    int remoteV1Id = v1Id == receivedEdge.destVert1Id ? receivedEdge.destVert1Id : receivedEdge.destVert2Id;
    int remoteV2Id = v1Id == receivedEdge.destVert1Id ? receivedEdge.destVert2Id : receivedEdge.destVert1Id;
    EXPECT_EQ(v1Id, remoteV1Id);
    EXPECT_EQ(v2Id, remoteV2Id);
  }

  for (auto edge : mesh->get_edges())
    if (edge && edge->count_remote_shared_entities() > 0)
    {
      EXPECT_EQ(receivedInfo.count(edge->get_id()), 1u);
    }
}

void check_coordinates(std::map<int, utils::Point>& expectedCoordsPerRank, std::vector<RemoteInfo>& remoteInfos)
{
  for (auto& remoteInfo : remoteInfos)
  {
    utils::Point expectedCoords = expectedCoordsPerRank[remoteInfo.srcRank];
    double dist                 = std::sqrt(dot(expectedCoords - remoteInfo.pt, expectedCoords - remoteInfo.pt));
    EXPECT_NEAR(dist, 0.0, 1e-12);
  }
}

} // namespace

TEST(Mesh, ParallelCounts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX      = 2;
  spec.numelY      = 1;
  spec.xmin        = 0;
  spec.xmax        = 1;
  spec.ymin        = 0;
  spec.ymax        = 1;
  double xMidpoint = 0.5 * (spec.xmin + spec.xmax);

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  EXPECT_EQ(count_valid(mesh->get_vertices()), 4);
  EXPECT_EQ(count_valid(mesh->get_edges()), 4);
  EXPECT_EQ(count_valid(mesh->get_elements()), 1);

  auto lowerVert = get_closest_vert(mesh, {xMidpoint, spec.ymin});
  auto upperVert = get_closest_vert(mesh, {xMidpoint, spec.ymax});

  int otherRank = 1 - utils::impl::comm_rank(mesh->get_comm());
  EXPECT_EQ(lowerVert->count_remote_shared_entities(), 1);
  EXPECT_EQ(lowerVert->get_remote_shared_entity(0).remoteRank, otherRank);

  EXPECT_EQ(upperVert->count_remote_shared_entities(), 1);
  EXPECT_EQ(lowerVert->get_remote_shared_entity(0).remoteRank, otherRank);

  test_vert_shared_entities(mesh);
  test_edge_shared_entities(mesh);
}

TEST(Mesh, RemoteSharedEntities)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 8)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 8;
  spec.numelY = 8;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);
  check_topology(mesh);

  test_vert_shared_entities(mesh);
  test_edge_shared_entities(mesh);

  if (utils::impl::comm_size(MPI_COMM_WORLD) == 4 || utils::impl::comm_size(MPI_COMM_WORLD) == 8)
  {
    auto vert = get_closest_vert(mesh, {(spec.xmin + spec.xmax) / 2, (spec.ymin + spec.ymax) / 2, 0});
    EXPECT_EQ(vert->count_remote_shared_entities(), 3);
  }
}

TEST(Mesh, RemoteSharedEntitiesPeriodicSquare)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX    = 4;
  spec.numelY    = 4;
  spec.xmin      = 0;
  spec.xmax      = 1;
  spec.ymin      = 0;
  spec.ymax      = 1;
  spec.xPeriodic = true;
  spec.yPeriodic = true;
  double xmid    = 0.5 * (spec.xmax + spec.xmin);
  double ymid    = 0.5 * (spec.ymax + spec.ymin);

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);
  check_topology(mesh);

  std::map<int, std::vector<RemoteInfo>> remoteInfos = get_vert_shared_entity_info(mesh);

  if (utils::impl::comm_rank(MPI_COMM_WORLD) == 0)
  {
    auto vertLL = get_closest_vert(mesh, {spec.xmin, spec.ymin});
    auto vertLR = get_closest_vert(mesh, {xmid, spec.ymin});
    auto vertTL = get_closest_vert(mesh, {spec.xmin, ymid});
    auto vertTR = get_closest_vert(mesh, {xmid, ymid});

    EXPECT_EQ(vertLL->count_remote_shared_entities(), 3);
    EXPECT_EQ(vertLR->count_remote_shared_entities(), 3);
    EXPECT_EQ(vertTL->count_remote_shared_entities(), 3);
    EXPECT_EQ(vertTR->count_remote_shared_entities(), 3);

    std::map<int, utils::Point> expectedCoordsLl = {{1, utils::Point(spec.xmax, spec.ymin)},
                                                    {2, utils::Point(spec.xmin, spec.ymax)},
                                                    {3, utils::Point(spec.xmax, spec.ymax)}};

    std::map<int, utils::Point> expectedCoordsLr = {
        {1, utils::Point(xmid, spec.ymin)}, {2, utils::Point(xmid, spec.ymax)}, {3, utils::Point(xmid, spec.ymax)}};

    std::map<int, utils::Point> expectedCoordsTl = {
        {1, utils::Point(spec.xmax, ymid)}, {2, utils::Point(spec.xmin, ymid)}, {3, utils::Point(spec.xmax, ymid)}};

    std::map<int, utils::Point> expectedCoordsTr = {
        {1, utils::Point(xmid, ymid)}, {2, utils::Point(xmid, ymid)}, {3, utils::Point(xmid, ymid)}};

    check_coordinates(expectedCoordsLl, remoteInfos[vertLL->get_id()]);
    check_coordinates(expectedCoordsLr, remoteInfos[vertLR->get_id()]);
    check_coordinates(expectedCoordsTl, remoteInfos[vertTL->get_id()]);
    check_coordinates(expectedCoordsTr, remoteInfos[vertTR->get_id()]);
  }
}

TEST(Mesh, RemoteSharedEntitiesPeriodicAnnulus)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 9)
    GTEST_SKIP();

  std::shared_ptr<Mesh> mesh = make_annulus_mesh(9, 9, 0.5, 1.5, 0);
  test_vert_shared_entities(mesh, true);
  test_edge_shared_entities(mesh);
}

TEST(Mesh, MPIStuff)
{
  int npts = 2;
  std::array<int, 2> gridShape{0, 0};
  MPI_Dims_create(npts, 2, gridShape.data());

  EXPECT_EQ(npts, gridShape[0] * gridShape[1]);
}

TEST(Mesh, ParallelEdgeOrientation)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto mesh = make_empty_mesh(MPI_COMM_WORLD);
  if (utils::impl::comm_rank(mesh->get_comm()) == 0)
  {
    MeshEntityPtr v1 = mesh->create_vertex(0, 0);
    MeshEntityPtr v2 = mesh->create_vertex(0, 1);
    MeshEntityPtr v3 = mesh->create_vertex(-0.5, 0.5);
    MeshEntityPtr e1 = mesh->create_edge(v1, v2);
    MeshEntityPtr e2 = mesh->create_edge(v2, v3);
    MeshEntityPtr e3 = mesh->create_edge(v3, v1);
    mesh->create_triangle(e1, e2, e3);

    v1->add_remote_shared_entity({1, 0});
    v2->add_remote_shared_entity({1, 2});
    e1->add_remote_shared_entity({1, 2});
  } else
  {
    MeshEntityPtr v1 = mesh->create_vertex(0, 0);
    MeshEntityPtr v2 = mesh->create_vertex(0.5, 0.5);
    MeshEntityPtr v3 = mesh->create_vertex(0, 1);
    MeshEntityPtr e1 = mesh->create_edge(v1, v2);
    MeshEntityPtr e2 = mesh->create_edge(v2, v3);
    MeshEntityPtr e3 = mesh->create_edge(v3, v1);
    mesh->create_triangle(e1, e2, e3);
    v1->add_remote_shared_entity({0, 0});
    v3->add_remote_shared_entity({0, 1});
    e3->add_remote_shared_entity({0, 0});
  }

  if (utils::impl::comm_rank(mesh->get_comm()) > 0)
    EXPECT_ANY_THROW(check_topology(mesh));
  else
    EXPECT_NO_THROW(check_topology(mesh));
}

TEST(Mesh, GetOwnerWithRemotes)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);
  auto vert = mesh->create_vertex(utils::Point(0, 0, 0));

  vert->add_remote_shared_entity({1, 7});
  vert->add_remote_shared_entity({2, 8});

  int owner = std::min(utils::impl::comm_rank(MPI_COMM_WORLD), 1);
  EXPECT_EQ(mesh::get_owner(mesh, vert), owner);
}

TEST(Mesh, GetOwnerWithoutRemotes)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);
  auto vert = mesh->create_vertex(utils::Point(0, 0, 0));
  EXPECT_EQ(mesh::get_owner(mesh, vert), utils::impl::comm_rank(MPI_COMM_WORLD));
}

TEST(Mesh, GetRemoteSharedEntity)
{
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);
  auto vert = mesh->create_vertex(utils::Point(0, 0, 0));

  vert->add_remote_shared_entity({1, 7});
  vert->add_remote_shared_entity({2, 8});

  EXPECT_EQ(mesh::get_remote_shared_entity(vert, 1).remoteRank, 1);
  EXPECT_EQ(mesh::get_remote_shared_entity(vert, 1).remoteId, 7);

  EXPECT_EQ(mesh::get_remote_shared_entity(vert, 2).remoteRank, 2);
  EXPECT_EQ(mesh::get_remote_shared_entity(vert, 2).remoteId, 8);

  EXPECT_ANY_THROW(mesh::get_remote_shared_entity(vert, 0));
  EXPECT_ANY_THROW(mesh::get_remote_shared_entity(vert, 3));
}


TEST(Mesh, ErrorRemotesNonSymmetric)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto mesh = make_empty_mesh(MPI_COMM_WORLD);
  if (utils::impl::comm_rank(MPI_COMM_WORLD) == 0)
  {
    auto v1 = mesh->create_vertex(0, 0, 0);
    auto v2 = mesh->create_vertex(0, 1, 0);
    auto v3 = mesh->create_vertex(-1, 0.5);
    mesh->create_triangle_from_verts(v1, v2, v3);

    v1->add_remote_shared_entity({1, 0});
  } else
  {
    auto v1 = mesh->create_vertex(0, 0,   0);
    auto v2 = mesh->create_vertex(1, 0.5, 0);
    auto v3 = mesh->create_vertex(0, 1,   0);
    mesh->create_triangle_from_verts(v1, v2, v3);

    v3->add_remote_shared_entity({0, 2});
  }

  EXPECT_ANY_THROW(check_topology(mesh));
}

TEST(Mesh, ErrorRemotesNotUnique)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

#ifndef NDEBUG
    GTEST_SKIP();
#endif
  auto mesh = make_empty_mesh(MPI_COMM_WORLD);
  auto v1 = mesh->create_vertex(0, 0, 0);
  auto v2 = mesh->create_vertex(0, 1, 0);
  auto v3 = mesh->create_vertex(-1, 0.5);
  mesh->create_triangle_from_verts(v1, v2, v3);

  v1->add_remote_shared_entity({1, 0});
  v1->add_remote_shared_entity({1, 0});

  EXPECT_ANY_THROW(check_topology(mesh));
}

TEST(Mesh, AnnulusRemotes)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
    GTEST_SKIP();

  // Note: 2x2 doesnt work, find out why
  std::shared_ptr<mesh::Mesh> mesh1 = make_annulus_mesh(4, 4, 0.5, 1.5, 0);
  mesh::check_topology(mesh1);
}

TEST(Mesh, getLocalId)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshSpec spec;
  spec.numelX = 1;
  spec.numelY = 1;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<Mesh> mesh = create_mesh(spec, func);

  mesh::MeshEntityPtr el = mesh->get_elements()[0];
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  mesh::get_downward(el, 0, verts.data());

  for (int i=0; i < 4; ++i)
  {
    EXPECT_EQ(mesh::get_local_id(el, verts[i]), i);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
