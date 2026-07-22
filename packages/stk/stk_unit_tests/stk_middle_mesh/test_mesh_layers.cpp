#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_layers.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace mesh;
using namespace mesh::impl;

namespace {

MeshEntityPtr find_closest_el(const std::vector<MeshEntityPtr>& entities, const utils::Point& pt, double tol)
{
  MeshEntityPtr entity = nullptr;

  double min_dist = tol;
  for (auto& e : entities)
  {
    auto disp = e->get_point_orig(0) - pt;
    double dist = std::sqrt(dot(disp, disp));
    if (dist < min_dist)
    {
      min_dist = dist;
      entity = e;
    }
  }

  return entity;
}

class AlwaysTrue
{
  public:
    bool operator()(MeshEntityPtr) { return true; }
};

class MeshLayersSerialNoFilterTester : public ::testing::Test
{
  protected:
    void SetUp() override
    {
      if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
        GTEST_SKIP();

      MeshSpec spec;
      spec.numelX = 3;
      spec.numelY = 3;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;

      auto func = [&](const utils::Point& pt) { return pt; };

      mesh = create_mesh(spec, func);   

      for (auto& v : mesh->get_vertices())
        if (v)
        {
          auto pt = v->get_point_orig(0);
          if (std::abs(pt.x) < 1e-13)
            roots.push_back(v);
        }         
    }

    std::shared_ptr<Mesh> mesh;
    std::vector<MeshEntityPtr> roots;
    AlwaysTrue filter;
};


class ExcludedVertFilter
{
  public:
    ExcludedVertFilter(MeshEntityPtr vert = nullptr) :
      m_entity(vert)
    {}

    bool operator()(MeshEntityPtr v) { return v != m_entity; }

  private:
    MeshEntityPtr m_entity;
};

class MeshLayersSerialFilterTester : public ::testing::Test
{
  protected:
    void SetUp() override
    {
      if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
        GTEST_SKIP();

      MeshSpec spec;
      spec.numelX = 3;
      spec.numelY = 3;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;

      auto func = [&](const utils::Point& pt) { return pt; };

      mesh = create_mesh(spec, func);   

      std::vector<MeshEntityPtr> output;
      MeshEntityPtr vExclude = nullptr;
      for (auto& v : mesh->get_vertices())
        if (v)
        {
          auto pt = v->get_point_orig(0);
          if (std::abs(pt.x) < 1e-13)
            roots.push_back(v);

          if (std::abs(pt.x - 1.0 / 3.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13)
            vExclude = v;
        }

      filter = ExcludedVertFilter(vExclude);     
    }

    std::shared_ptr<Mesh> mesh;
    std::vector<MeshEntityPtr> roots;
    ExcludedVertFilter filter;
};


class MeshLayersParallelNoFilterTester : public ::testing::Test
{
  protected:
    void SetUp() override
    {
      if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
        GTEST_SKIP();

      MeshSpec spec;
      spec.numelX = 4;
      spec.numelY = 4;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;

      auto func = [&](const utils::Point& pt) { return pt; };

      mesh = create_mesh(spec, func);

      MeshLayers layers(mesh);

      myRank = utils::impl::comm_rank(MPI_COMM_WORLD);
      if (myRank == 0)
        for (auto& v : mesh->get_vertices())
          if (v)
          {
            auto pt = v->get_point_orig(0);
            if (std::abs(pt.y) < 1e-13)
              roots.push_back(v);
          }
    }

    std::shared_ptr<Mesh> mesh;
    std::vector<MeshEntityPtr> roots;
    AlwaysTrue filter;
    int myRank = utils::impl::comm_rank(MPI_COMM_WORLD);
};

class MeshLayersParallelFilterTester : public ::testing::Test
{
  protected:
    void SetUp() override
    {
      if (utils::impl::comm_size(MPI_COMM_WORLD) > 4)
        GTEST_SKIP();

      MeshSpec spec;
      spec.numelX = 4;
      spec.numelY = 4;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;

      auto func = [&](const utils::Point& pt) { return pt; };

      mesh = create_mesh(spec, func);

      MeshLayers layers(mesh);

      myRank = utils::impl::comm_rank(MPI_COMM_WORLD);
      if (myRank == 0)
        for (auto& v : mesh->get_vertices())
          if (v)
          {
            auto pt = v->get_point_orig(0);
            if (std::abs(pt.y) < 1e-13)
              roots.push_back(v);
          }

      filter = ExcludedVertFilter(find_closest_el(mesh->get_vertices(), {0.5, 0.5}, 1e-10));
    }

    std::shared_ptr<Mesh> mesh;
    std::vector<MeshEntityPtr> roots;
    ExcludedVertFilter filter;
    int myRank = utils::impl::comm_rank(MPI_COMM_WORLD);
};

}

TEST_F(MeshLayersSerialNoFilterTester, OneLayer)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 1, output);

  EXPECT_EQ(output.size(), static_cast<unsigned int>(8));
  for (auto v : output)
  {
    double xCoord  = v->get_point_orig(0).x;
    bool isCorrect = std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13;
    EXPECT_TRUE(isCorrect);
  }
}


TEST_F(MeshLayersSerialNoFilterTester, TwoLayers)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 2, output);

  EXPECT_EQ(output.size(), static_cast<unsigned int>(12));
  for (auto v : output)
  {
    double xCoord = v->get_point_orig(0).x;
    bool isCorrect =
        std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13 || std::abs(xCoord - 2.0 / 3.0) < 1e-13;
    EXPECT_TRUE(isCorrect);
  }
}


TEST_F(MeshLayersSerialNoFilterTester, ThreeLayers)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 3, output);
  EXPECT_EQ(output.size(), static_cast<unsigned int>(16));
  EXPECT_TRUE(is_unique(output));
}

TEST_F(MeshLayersSerialFilterTester, OneLayer)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 1, output);

  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(7));
  for (auto v : output)
  {
    double xCoord  = v->get_point_orig(0).x;
    double yCoord  = v->get_point_orig(0).y;
    bool isCorrect = std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13;

    if (std::abs(xCoord - 1.0 / 3.0) < 1e-13 && std::abs(yCoord - 1.0 / 3.0) < 1e-13)
      EXPECT_FALSE(true);
    else
      EXPECT_TRUE(isCorrect);
  }
}

TEST_F(MeshLayersSerialFilterTester, TwoLayers)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 2, output);

  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(10));
  for (auto v : output)
  {
    double xCoord = v->get_point_orig(0).x;
    double yCoord = v->get_point_orig(0).y;
    bool isCorrect =
        std::abs(xCoord) < 1e-13 || std::abs(xCoord - 1.0 / 3.0) < 1e-13 || std::abs(xCoord - 2.0 / 3.0) < 1e-13;

    bool isExcluded = (std::abs(xCoord - 1.0 / 3.0) < 1e-13 && std::abs(yCoord - 1.0 / 3.0) < 1e-13) ||
                      (std::abs(xCoord - 2.0 / 3.0) < 1e-13 && std::abs(yCoord - 1.0 / 3.0) < 1e-13);

    if (isExcluded)
      EXPECT_FALSE(true);
    else
      EXPECT_TRUE(isCorrect);
  }
}

TEST_F(MeshLayersSerialFilterTester, ThreeLayers)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 3, output);

  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(14));

  std::set<MeshEntityPtr> verts;
  for (auto v : mesh->get_vertices())
  {
    auto pt         = v->get_point_orig(0);
    bool isExcluded = (std::abs(pt.x - 1.0 / 3.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13) ||
                      (std::abs(pt.x - 1.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13);
    if (!isExcluded)
      verts.insert(v);
  }

  for (auto v : output)
  {
    EXPECT_TRUE(verts.count(v) == 1);
  }
}

TEST_F(MeshLayersSerialFilterTester, FourLayers)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 4, output);
  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(15));

  std::set<MeshEntityPtr> verts;
  for (auto v : mesh->get_vertices())
  {
    auto pt         = v->get_point_orig(0);
    bool isExcluded = (std::abs(pt.x - 1.0 / 3.0) < 1e-13 && std::abs(pt.y - 1.0 / 3.0) < 1e-13);
    if (!isExcluded)
      verts.insert(v);
  }

  for (auto v : output)
  {
    EXPECT_TRUE(verts.count(v) == 1);
  }
}


TEST_F(MeshLayersSerialFilterTester, AllLayers)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_all_layers(filter, roots, output);
  EXPECT_TRUE(is_unique(output));
  EXPECT_EQ(output.size(), static_cast<unsigned int>(15));
}


TEST_F(MeshLayersParallelNoFilterTester, OneLayerNP2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 1, output);

  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25} };
  } else
  {
    expectedPts = { {0.5, 0},    {0.75, 0}, {0.5, 0.25}};  
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));  
}

TEST_F(MeshLayersParallelNoFilterTester, TwoLayersNP2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 2, output);

  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25},
                    {0, 0.50}, {0.25, 0.50}, {0.5, 0.50} };
  } else
  {
    expectedPts = { {0.5, 0}, {0.75, 0},    {0.5, 0.25},
                    {1, 0},   {0.75, 0.25}, {0.5, 0.5}};
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));  
}

TEST_F(MeshLayersParallelNoFilterTester, ThreeLayersNP2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 3, output);

  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25},
                    {0, 0.50}, {0.25, 0.50}, {0.5, 0.50},
                    {0, 0.75}, {0.25, 0.75}, {0.5, 0.75} };
  } else
  {
    expectedPts = { {0.5, 0},  {0.75, 0},    {0.5, 0.25},
                    {1, 0},    {0.75, 0.25}, {0.5, 0.5},
                    {1, 0.25}, {0.75, 0.50}, {0.5, 0.75}}; 
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));   
}

TEST_F(MeshLayersParallelNoFilterTester, AllLayersNP2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_all_layers(filter, roots, output);


  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25},
                    {0, 0.50}, {0.25, 0.50}, {0.5, 0.50},
                    {0, 0.75}, {0.25, 0.75}, {0.5, 0.75},
                    {0, 1.00}, {0.25, 1.00}, {0.5, 1.00} };
  } else
  {
    expectedPts = { {0.5, 0},    {0.75, 0},    {1.00, 0},
                    {0.5, 0.25}, {0.75, 0.25}, {1.00, 0.25},
                    {0.5, 0.50}, {0.75, 0.50}, {1.00, 0.50},
                    {0.5, 0.75}, {0.75, 0.75}, {1.00, 0.75},
                    {0.5, 1.00}, {0.75, 1.00}, {1.00, 1.00}};  
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));   
}

TEST_F(MeshLayersParallelNoFilterTester, TwoLayersNP4)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 2, output);

  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25},
                    {0, 0.50}, {0.25, 0.50}, {0.5, 0.50} };
  } else if (myRank == 1)
  {
    expectedPts = { {0.5, 0}, {0.75, 0},    {0.5, 0.25},
                    {1, 0},   {0.75, 0.25}, {0.5, 0.5}};  
  } else if (myRank == 2)
  {
    expectedPts = { {0, 0.5}, {0.25, 0.5}, {0.5, 0.5}};
  } else if (myRank == 3)
  {
    expectedPts = { {0.5, 0.5}};
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));  
}

TEST_F(MeshLayersParallelNoFilterTester, ThreeLayersNP4)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 3, output);

  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25},
                    {0, 0.50}, {0.25, 0.50}, {0.5, 0.50} };
  } else if (myRank == 1)
  {
    expectedPts = { {0.5, 0}, {0.75, 0},    {0.5, 0.25},
                    {1, 0},   {0.75, 0.25}, {0.5, 0.5},
                    {0.75, 0.5}, {1.0, 0.25}};  
  } else if (myRank == 2)
  {
    expectedPts = { {0, 0.5},  {0.25, 0.5},  {0.5, 0.5},
                    {0, 0.75}, {0.25, 0.75}, {0.5, 0.75}};
  } else if (myRank == 3)
  {
    expectedPts = { {0.5, 0.5}, {0.75, 0.5}, {0.5, 0.75}};
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));  
}

TEST_F(MeshLayersParallelFilterTester, ThreeLayersNP4)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  MeshLayers layers(mesh);
  std::vector<MeshEntityPtr> output;
  layers.get_layers(filter, roots, 3, output);

  std::vector<utils::Point> expectedPts;
  if (myRank == 0)
  {
    expectedPts = { {0, 0},    {0.25, 0},    {0.5, 0},
                    {0, 0.25}, {0.25, 0.25}, {0.5, 0.25},
                    {0, 0.50}, {0.25, 0.50} };
  } else if (myRank == 1)
  {
    expectedPts = { {0.5, 0}, {0.75, 0},    {0.5, 0.25},
                    {1, 0},   {0.75, 0.25},
                    {0.75, 0.5}, {1.0, 0.25}};  
  } else if (myRank == 2)
  {
    expectedPts = { {0, 0.5},  {0.25, 0.5},
                    {0, 0.75}, {0.25, 0.75} };
  } else if (myRank == 3)
  {
    expectedPts = {{0.75, 0.5}};
  }

  EXPECT_EQ(output.size(), expectedPts.size());
  for (auto& pt : expectedPts)
    EXPECT_TRUE(find_closest_el(output, pt, 1e-10));  
}


} // namespace impl
} // namespace middle_mesh
} // namespace stk
