#include "gtest/gtest.h"

#include "stk_middle_mesh/adjacency_search.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/predicates/intersection.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

std::shared_ptr<mesh::Mesh> create_mesh(const std::vector<utils::Point>& vertPts, const std::vector<std::array<int, 4>>& elVertIds)
{
  auto mesh = mesh::make_empty_mesh();

  std::vector<mesh::MeshEntityPtr> verts;
  for (auto& pt : vertPts)
    verts.push_back(mesh->create_vertex(pt));

  for (auto& vertIds : elVertIds)
  {
    if (vertIds[3] == -1)
    {
      mesh->create_triangle_from_verts(verts[vertIds[0]], verts[vertIds[1]], verts[vertIds[2]]);
    } else
    {
      mesh->create_quad_from_verts(verts[vertIds[0]], verts[vertIds[1]], verts[vertIds[2]], verts[vertIds[3]]);
    }
  }

  return mesh;
}

bool any_contains(const std::vector<mesh::MeshEntityPtr>& entities, mesh::MeshEntityPtr el1)
{
  predicates::impl::FaceContainer container;
  utils::Point pt = compute_centroid(el1);
  for (auto e : entities)
    if (container.contains(e, pt))
      return true;

  return false;
}

// test that "every" point of el1 is contained in one of entities
// In practice, we don't test every point, just a bunch of them
void test_contains_inverse(const std::vector<mesh::MeshEntityPtr>& entities, mesh::MeshEntityPtr el1)
{
  int npts   = 5; // number of test points
  double dxi = 1.0 / npts;

  predicates::impl::FaceContainer container(1e-12);
  for (int i = 0; i < npts; ++i)
    for (int j = 0; j < npts; ++j)
    {
      double xi1      = i * dxi;
      double xi2      = j * dxi;
      utils::Point pt = compute_coords_from_xi(el1, utils::Point(xi1, xi2));

      bool found = false;
      for (auto e : entities)
      {
        if (container.contains(e, pt))
        {
          found = true;
          break;
        }
      }

      EXPECT_TRUE(found);
    }
}

}

TEST(AdjacencySearch, Identical)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // mesh1 and mesh2 are identical
  mesh::impl::MeshSpec spec;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec, func);
  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    // the algorithm could produce either the single element that overlaps
    // el1, or that element plus all adjacencies, or something in between
    EXPECT_GE(entities.size(), static_cast<unsigned int>(1));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(9));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

TEST(AdjacencySearch, Subdivision)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 10;
  spec2.numelY = 10;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);
  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    // the algorithm could produce either the single element that overlaps
    // el1, or that element plus all adjacencies, or something in between
    EXPECT_GE(entities.size(), static_cast<unsigned int>(4));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(16));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

TEST(AdjacencySearch, ContainedElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // mesh1 and mesh2 are identical
  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 15;
  spec2.numelY = 15;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    EXPECT_GE(entities.size(), static_cast<unsigned int>(9));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(9 + 16));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

TEST(AdjacencySearch, ContainedElements)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // elements on mesh2 are contained within mesh1
  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 5;
  spec.numelY = 5;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 20;
  spec2.numelY = 20;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    EXPECT_GE(entities.size(), static_cast<unsigned int>(16));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(16 + 20));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

TEST(AdjacencySearch, TwoToFive)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // in this case, there is no mesh2 edge coincident with any internal edges,
  // and there are contained elements
  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 10;
  spec.numelY = 10;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 25;
  spec2.numelY = 25;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    EXPECT_GE(entities.size(), static_cast<unsigned int>(9));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(9 + 16));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

TEST(AdjacencySearch, NoVerts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // for some elements on mesh1, there are no mesh2 vertices contained
  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 12;
  spec.numelY = 12;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 18;
  spec2.numelY = 4;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    EXPECT_GE(entities.size(), static_cast<unsigned int>(2));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(11));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

TEST(AdjacencySearch, CutCorner)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec;
  spec.numelX = 1;
  spec.numelY = 1;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);

  auto mesh2 = mesh::make_empty_mesh();
  auto v1    = mesh2->create_vertex(-0.25, 0.25);
  auto v2    = mesh2->create_vertex(0.5, 1.25);
  auto v3    = mesh2->create_vertex(-0.25, 1.25);
  auto v4    = mesh2->create_vertex(0, 0);
  auto v5    = mesh2->create_vertex(1, 1);
  auto v6    = mesh2->create_vertex(1, 0);
  mesh2->create_triangle_from_verts(v1, v2, v3);
  mesh2->create_quad_from_verts(v4, v5, v2, v1);
  mesh2->create_triangle_from_verts(v4, v6, v5);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);
  std::vector<mesh::MeshEntityPtr> entities;
  auto el1 = search.get_next(entities);
  // std::cout << "mesh2 elements:" << std::endl;
  // for (auto& el : entities)
  //   std::cout << el << std::endl;
  EXPECT_EQ(entities.size(), static_cast<unsigned int>(3));
  EXPECT_TRUE(any_contains(entities, el1));
  test_contains_inverse(entities, el1);

  // TODO: this avoids a memory leak of the Projection created
  //        inside the adjacency search, figure out why
  while ((el1 = search.get_next(entities)))
  {}
}

TEST(AdjacencySearch, Contains)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  // some elements on mesh1 are completely contained in elements on mesh2
  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 36;
  spec.numelY = 36;
  spec.xmin   = 0;
  spec.xmax   = 1;
  spec.ymin   = 0;
  spec.ymax   = 1;

  spec2.numelX = 12;
  spec2.numelY = 12;
  spec2.xmin   = 0;
  spec2.xmax   = 1;
  spec2.ymin   = 0;
  spec2.ymax   = 1;

  auto func = [&](const utils::Point& pt) { return pt; };

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ((el1 = search.get_next(entities)))
  {
    EXPECT_GE(entities.size(), static_cast<unsigned int>(1));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(4));
    EXPECT_TRUE(any_contains(entities, el1));
    test_contains_inverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX * spec.numelY);
}

/*
TEST(AdjacencySearch, Contains2)
{
  // some elements on mesh1 are completely contained in elements on mesh2
  mesh::impl::MeshSpec spec, spec2;
  spec.numelX = 200;
  spec.numelY = 200;
  spec.xmin = 0;
  spec.xmax = 1;
  spec.ymin = 0;
  spec.ymax = 1;

  spec2.numelX = 400;
  spec2.numelY = 400;
  spec2.xmin = 0;
  spec2.xmax = 1;
  spec2.ymin = 0;
  spec2.ymax = 1;

  auto func = [&](const utils::Point& pt) { return pt;};

  std::shared_ptr<mesh::Mesh> mesh1 = create_mesh(spec, func);
  std::shared_ptr<mesh::Mesh> mesh2 = create_mesh(spec2, func);

  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  int nentities = 0;
  while ( (el1 = search.get_next(entities)) )
  {
    EXPECT_GE(entities.size(), static_cast<unsigned int>(1));
    EXPECT_LE(entities.size(), static_cast<unsigned int>(100));
    //EXPECT_TRUE(anyContains(entities, el1));
    //testContainsInverse(entities, el1);
    nentities++;
  }
  EXPECT_EQ(nentities, spec.numelX*spec.numelY);
}
*/

TEST(AdjacencySearch, UnreachableCutCorner)
{
  auto mesh1 = create_mesh( {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}},
                            {{0, 1, 2, 3}});

  auto mesh2 = create_mesh( {{0, 0, 0}, {0.25, 0, 0}, {0.75, 0, 0}, {1.25, 0, 0},
                             {0, 0.75, 0}, {0.25, 0.75, 0}, {0.75, 0.75, 0}, {1.25, 0.75, 0},
                             {0, 1, 0}, {0.25, 1, 0}, {0.5, 1.25, 0}, {0.9, 1.0, 0},
                             {1.1, 1.25, 0}, {1.5, 0.9, 0}},
                            {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6},
                             {4, 5, 9, 8}, {5, 6, 10, 9}, {6, 7, 11, 10},
                             {7, 13, 12, 11}});

                               
  mesh::impl::AdjacencySearch search(mesh1, mesh2);

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> entities;
  while ((el1 = search.get_next(entities)))
  {
    EXPECT_EQ(entities.size(), 7u);

    auto pred = [](const mesh::MeshEntityPtr& e) { return e->get_id() == 6; };
    EXPECT_NE(std::find_if(entities.begin(), entities.end(), pred), entities.end());
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
