#include "stk_middle_mesh/adjacency_search_closest_points.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {
void check_vectors_equal(const std::vector<mesh::MeshEntityPtr>& vector1,
                         const std::vector<mesh::MeshEntityPtr>& vector2)
{
  std::vector<mesh::MeshEntityPtr> vector1Sorted(vector1), vector2Sorted(vector2);
  std::sort(vector1Sorted.begin(), vector1Sorted.end(), mesh::is_less);
  std::sort(vector2Sorted.begin(), vector2Sorted.end(), mesh::is_less);

  EXPECT_EQ(vector1.size(), vector2.size());
  for (size_t i = 0; i < vector1.size(); ++i)
    EXPECT_EQ(vector1Sorted[i], vector2Sorted[i]);
}

mesh::MeshEntityPtr get_closest_vert(std::shared_ptr<mesh::Mesh> mesh, const utils::Point& pt)
{
  double distMin                = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr entityMin = nullptr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point disp = vert->get_point_orig(0) - pt;
      double dist       = dot(disp, disp);
      if (dist < distMin)
      {
        distMin   = dist;
        entityMin = vert;
      }
    }

  return entityMin;
}
} // namespace

// Test cases:
//    1. uniform mesh, test radiating outward
//    2. graded mesh, test that points returned really are the closest
//    3. test filtering

TEST(AdjacencySearchClosestPoints, UniformMesh)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.ymin   = 0;
  spec.xmax   = 1;
  spec.ymax   = 1;

  auto mesh = mesh::impl::create_mesh(spec, [](const utils::Point& pt) { return pt; });

  mesh::impl::AdjacencySearchClosestPoints searcher(mesh);
  std::vector<mesh::MeshEntityPtr> closestVerts;

  mesh::MeshEntityPtr startVert = get_closest_vert(mesh, {0.5, 0.5});

  {
    searcher.search(startVert, 4, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {0.5, 0.25}), get_closest_vert(mesh, {0.75, 0.5}), get_closest_vert(mesh, {0.5, 0.75}),
        get_closest_vert(mesh, {0.25, 0.5})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 8, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {0.5, 0.25}),  get_closest_vert(mesh, {0.75, 0.5}),
        get_closest_vert(mesh, {0.5, 0.75}),  get_closest_vert(mesh, {0.25, 0.5}),
        get_closest_vert(mesh, {0.25, 0.25}), get_closest_vert(mesh, {0.75, 0.25}),
        get_closest_vert(mesh, {0.75, 0.75}), get_closest_vert(mesh, {0.25, 0.75}),
    };

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 25, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts;
    double deltaX = (spec.xmax - spec.xmin) / spec.numelX;
    double deltaY = (spec.ymax - spec.ymin) / spec.numelX;
    for (int i = 0; i < spec.numelX + 1; ++i)
      for (int j = 0; j < spec.numelY + 1; ++j)
      {
        mesh::MeshEntityPtr vert = get_closest_vert(mesh, {spec.xmin + deltaX * i, spec.ymin + deltaY * j});
        if (vert != startVert)
          expectedVerts.push_back(vert);
      }

    check_vectors_equal(closestVerts, expectedVerts);
  }
}

TEST(AdjacencySearchClosestPoints, GradedMesh)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.ymin   = 0;
  spec.xmax   = 1;
  spec.ymax   = 1;

  auto mesh =
      mesh::impl::create_mesh(spec, [](const utils::Point& pt) { return utils::Point(pt.x, pt.y * pt.y * pt.y, 0); });

  mesh::impl::AdjacencySearchClosestPoints searcher(mesh);
  std::vector<mesh::MeshEntityPtr> closestVerts;

  double x0 = 0, y0 = 0;
  double x1 = 0.25, y1 = 0.015625;
  double x2 = 0.50, y2 = 0.125;
  double y3 = 0.421875;

  mesh::MeshEntityPtr startVert = get_closest_vert(mesh, {x1, y2});

  {
    searcher.search(startVert, 2, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {get_closest_vert(mesh, {x1, y1}),
                                                      get_closest_vert(mesh, {x1, y0})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 4, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {x1, y1}), get_closest_vert(mesh, {x1, y0}), get_closest_vert(mesh, {x0, y2}),
        get_closest_vert(mesh, {x2, y2})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 6, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {x1, y1}), get_closest_vert(mesh, {x1, y0}), get_closest_vert(mesh, {x0, y2}),
        get_closest_vert(mesh, {x2, y2}), get_closest_vert(mesh, {x0, y1}), get_closest_vert(mesh, {x2, y1})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 8, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {x1, y1}), get_closest_vert(mesh, {x1, y0}), get_closest_vert(mesh, {x0, y2}),
        get_closest_vert(mesh, {x2, y2}), get_closest_vert(mesh, {x0, y1}), get_closest_vert(mesh, {x2, y1}),
        get_closest_vert(mesh, {x0, y0}), get_closest_vert(mesh, {x2, y0})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 9, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {x1, y1}), get_closest_vert(mesh, {x1, y0}), get_closest_vert(mesh, {x0, y2}),
        get_closest_vert(mesh, {x2, y2}), get_closest_vert(mesh, {x0, y1}), get_closest_vert(mesh, {x2, y1}),
        get_closest_vert(mesh, {x0, y0}), get_closest_vert(mesh, {x2, y0}), get_closest_vert(mesh, {x1, y3})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 11, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {x1, y1}), get_closest_vert(mesh, {x1, y0}), get_closest_vert(mesh, {x0, y2}),
        get_closest_vert(mesh, {x2, y2}), get_closest_vert(mesh, {x0, y1}), get_closest_vert(mesh, {x2, y1}),
        get_closest_vert(mesh, {x0, y0}), get_closest_vert(mesh, {x2, y0}), get_closest_vert(mesh, {x1, y3}),
        get_closest_vert(mesh, {x0, y3}), get_closest_vert(mesh, {x2, y3})};

    check_vectors_equal(closestVerts, expectedVerts);
  }
}

TEST(AdjacencySearchClosestPoints, UniformMeshFilter)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 4;
  spec.numelY = 4;
  spec.xmin   = 0;
  spec.ymin   = 0;
  spec.xmax   = 1;
  spec.ymax   = 1;

  auto mesh = mesh::impl::create_mesh(spec, [](const utils::Point& pt) { return pt; });

  auto filterFieldPtr = mesh::create_field<mesh::impl::AdjacencySearchClosestPoints::Bool>(
      mesh, mesh::impl::FieldShape(1, 0, 0), 1, true);
  auto& filterField                                      = *filterFieldPtr;
  filterField(get_closest_vert(mesh, {0.25, 0.5}), 0, 0) = false;

  mesh::impl::AdjacencySearchClosestPoints searcher(mesh, filterFieldPtr);
  std::vector<mesh::MeshEntityPtr> closestVerts;

  mesh::MeshEntityPtr startVert = get_closest_vert(mesh, {0.5, 0.5});

  {
    searcher.search(startVert, 3, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {0.5, 0.25}), get_closest_vert(mesh, {0.75, 0.5}), get_closest_vert(mesh, {0.5, 0.75})};

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 7, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts = {
        get_closest_vert(mesh, {0.5, 0.25}),  get_closest_vert(mesh, {0.75, 0.5}),
        get_closest_vert(mesh, {0.5, 0.75}),  get_closest_vert(mesh, {0.25, 0.25}),
        get_closest_vert(mesh, {0.75, 0.25}), get_closest_vert(mesh, {0.75, 0.75}),
        get_closest_vert(mesh, {0.25, 0.75}),
    };

    check_vectors_equal(closestVerts, expectedVerts);
  }

  {
    searcher.search(startVert, 25, closestVerts);
    std::vector<mesh::MeshEntityPtr> expectedVerts;
    double deltaX = (spec.xmax - spec.xmin) / spec.numelX;
    double deltaY = (spec.ymax - spec.ymin) / spec.numelX;
    for (int i = 0; i < spec.numelX + 1; ++i)
      for (int j = 0; j < spec.numelY + 1; ++j)
      {
        mesh::MeshEntityPtr vert = get_closest_vert(mesh, {spec.xmin + deltaX * i, spec.ymin + deltaY * j});
        if (vert != startVert && filterField(vert, 0, 0))
          expectedVerts.push_back(vert);
      }

    check_vectors_equal(closestVerts, expectedVerts);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
