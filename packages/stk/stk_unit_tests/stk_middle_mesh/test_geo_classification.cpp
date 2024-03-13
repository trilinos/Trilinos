
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_geo_classifier.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;
using namespace utils::impl;

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

} // namespace

TEST(GeoClassification, Square)
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

  MeshGeoClassifier classifier;
  classifier.classify(mesh);

  // test the dimensions are correct
  for (auto& v : mesh->get_vertices())
  {
    int dim = static_cast<int>(get_geo_classification(spec, v->get_point_orig(0)));
    EXPECT_EQ(mesh->get_geo_class(v).dim, dim);
  }

  // check vert ids are unique
  std::vector<MeshEntityPtr> verts;
  for (auto& v : mesh->get_vertices())
    if (mesh->get_geo_class(v).dim == 0)
      verts.push_back(v);

  for (auto& v1 : verts)
    for (auto& v2 : verts)
      if (v1 != v2)
      {
        EXPECT_FALSE(mesh->get_geo_class(v1).id == mesh->get_geo_class(v2).id);
      }

  // test edges
  for (auto& edge : mesh->get_edges())
  {
    auto g = mesh->get_geo_class(edge);
    if (edge->count_up() == 1)
    {
      auto pt1    = edge->get_down(0)->get_point_orig(0);
      auto pt2    = edge->get_down(1)->get_point_orig(0);
      auto pt     = 0.5 * (pt1 + pt2);
      double eps  = 1e-13;
      bool isXmin = std::abs(pt.x - spec.xmin) < eps;
      bool isXmax = std::abs(pt.x - spec.xmax) < eps;
      bool isYmin = std::abs(pt.y - spec.ymin) < eps;
      bool isYmax = std::abs(pt.y - spec.ymax) < eps;

      EXPECT_EQ(g.dim, 1);
      if (isYmin)
        EXPECT_EQ(g.id, 0);
      else if (isXmax)
        EXPECT_EQ(g.id, 1);
      else if (isYmax)
        EXPECT_EQ(g.id, 2);
      else if (isXmin)
        EXPECT_EQ(g.id, 3);
      else
        assert(false);
    } else
    {
      EXPECT_EQ(g, GeoClassification(2, 0));
    }
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
