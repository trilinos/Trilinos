#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/mesh_projection_calculator.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

using namespace nonconformal4::impl;
using namespace predicates::impl;

class MeshProjectionCalculatorTester : public ::testing::Test
{
  protected:
    MeshProjectionCalculatorTester()
      : m_mesh1(mesh::make_empty_mesh())
      , m_mesh2(mesh::make_empty_mesh())
      , m_meshIn(mesh::make_empty_mesh())
    {}

    void create_mesh1(const std::vector<utils::Point>& verts, const std::vector<std::array<int, 4>>& elementVerts)
    {
      create_mesh(m_mesh1, verts, elementVerts);
    }

    void create_mesh2(const std::vector<utils::Point>& verts, const std::vector<std::array<int, 4>>& elementVerts)
    {
      create_mesh(m_mesh2, verts, elementVerts);
    }

    void run()
    {
      double eps  = 1e-13;
      mClassifier = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(
          m_mesh2, PointClassifierNormalWrapperTolerances(eps));
      mRelationalData = std::make_shared<nonconformal4::impl::MeshRelationalData>(m_mesh1, m_mesh2, m_meshIn);
      mProjector      = std::make_shared<nonconformal4::impl::MeshProjectionCalculator>(
          m_mesh1, m_mesh2, mRelationalData, mClassifier, middle_mesh::impl::EdgeTracerTolerances(eps));
      mProjector->project();
    }

    std::set<nonconformal4::impl::FakeVert> get_all_fake_verts()
    {
      std::set<nonconformal4::impl::FakeVert> fakeVerts;

      auto& verts1ToFakeVerts = *(mRelationalData->verts1ToFakeVerts);
      for (auto& vert : m_mesh1->get_vertices())
        if (vert)
          fakeVerts.insert(verts1ToFakeVerts(vert, 0, 0));

      auto& verts2ToFakeVerts = *(mRelationalData->verts2ToFakeVerts);
      for (auto& vert : m_mesh2->get_vertices())
        if (vert)
          fakeVerts.insert(verts2ToFakeVerts(vert, 0, 0));

      auto& edges2ToFakeVertsIn = *(mRelationalData->edges2ToFakeVertsIn);
      for (auto& edge : m_mesh2->get_edges())
        if (edge)
        {
          for (int i = 0; i < edges2ToFakeVertsIn.get_num_comp(edge, 0); ++i)
            fakeVerts.insert(edges2ToFakeVertsIn(edge, 0, i).vert);
        }

      return fakeVerts;
    }

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> mClassifier;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> mRelationalData;
    std::shared_ptr<nonconformal4::impl::MeshProjectionCalculator> mProjector;

  private:
    void create_mesh(std::shared_ptr<mesh::Mesh> mesh, const std::vector<utils::Point>& vertCoords,
                     const std::vector<std::array<int, 4>>& elementVerts)
    {
      for (auto& vert : vertCoords)
        mesh->create_vertex(vert);

      auto& verts = mesh->get_vertices();
      for (auto& elVertIds : elementVerts)
      {
        if (elVertIds[3] == -1)
        {
          auto v0 = verts[elVertIds[0]];
          auto v1 = verts[elVertIds[1]];
          auto v2 = verts[elVertIds[2]];
          mesh->create_triangle_from_verts(v0, v1, v2);
        } else
        {
          auto v0 = verts[elVertIds[0]];
          auto v1 = verts[elVertIds[1]];
          auto v2 = verts[elVertIds[2]];
          auto v3 = verts[elVertIds[3]];
          mesh->create_quad_from_verts(v0, v1, v2, v3);
        }
      }
    }
};

void expect_near(const utils::Point& pt1, const utils::Point& pt2, double tol)
{
  EXPECT_NEAR(pt1.x, pt2.x, tol);
  EXPECT_NEAR(pt1.y, pt2.y, tol);
  EXPECT_NEAR(pt1.z, pt2.z, tol);
}

mesh::MeshEntityPtr find_closest_edge(std::shared_ptr<mesh::Mesh> mesh, const utils::Point& midpoint)
{
  double minDist              = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEdge = nullptr;
  for (auto& edge : mesh->get_edges())
    if (edge)
    {
      utils::Point edgeMidpoint = compute_centroid_3d(edge);
      utils::Point disp         = edgeMidpoint - midpoint;
      double dist               = std::sqrt(dot(disp, disp));

      if (dist < minDist)
      {
        minDist = dist;
        minEdge = edge;
      }
    }

  return minEdge;
}

void test_classification(std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalData,
                         std::shared_ptr<mesh::Mesh> mesh2,
                         const std::vector<predicates::impl::PointRecord>& expectedRecords)
{
  int idx                  = 0;
  auto& verts2ClassOnMesh1 = *(meshRelationalData->verts2ClassOnMesh1);
  for (auto& vert2 : mesh2->get_vertices())
    if (vert2)
    {
      const predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert2, 0, 0);
      EXPECT_EQ(record.type, expectedRecords[idx].type);
      EXPECT_EQ(get_entity(record), get_entity(expectedRecords[idx]));
      idx++;
    }
}

void test_edge_fake_verts(std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalData,
                          std::shared_ptr<mesh::Mesh> mesh2,
                          std::map<mesh::MeshEntityPtr, std::vector<utils::Point>>& edges2FakeVertsExpected)
{
  auto cmp = [](const utils::Point& lhs, const utils::Point& rhs) {
    if (lhs.x == rhs.x)
      return lhs.y < rhs.y;
    else
      return lhs.x < rhs.x;
  };

  auto& edges2ToFakeVertsIn = *(meshRelationalData->edges2ToFakeVertsIn);
  for (auto& edge : mesh2->get_edges())
    if (edge)
    {
      ASSERT_EQ(edges2ToFakeVertsIn.get_num_comp(edge, 0), (int)edges2FakeVertsExpected[edge].size());

      std::vector<utils::Point> edge2FakeVertsExpected = edges2FakeVertsExpected[edge];
      std::vector<utils::Point> edge2FakeVerts;
      for (int i = 0; i < edges2ToFakeVertsIn.get_num_comp(edge, 0); ++i)
        edge2FakeVerts.push_back(edges2ToFakeVertsIn(edge, 0, i).vert.pt);

      std::sort(edge2FakeVerts.begin(), edge2FakeVerts.end(), cmp);
      std::sort(edge2FakeVertsExpected.begin(), edge2FakeVertsExpected.end(), cmp);

      for (size_t i = 0; i < edge2FakeVerts.size(); ++i)
      {
        utils::Point disp = edge2FakeVerts[i] - edge2FakeVertsExpected[i];
        double dist       = std::sqrt(dot(disp, disp));
        EXPECT_NEAR(dist, 0.0, 1e-13);
      }
    }
}

void test_edge_splits(std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalData,
                      std::shared_ptr<mesh::Mesh> mesh1,
                      std::map<mesh::MeshEntityPtr, std::vector<double>>& mesh1EdgesToSplitExpected)
{
  auto& mesh1EdgesToSplit = *(meshRelationalData->mesh1EdgesToSplit);
  for (auto& edge : mesh1->get_edges())
    if (edge)
    {
      ASSERT_EQ(mesh1EdgesToSplit.get_num_comp(edge, 0), (int)mesh1EdgesToSplitExpected[edge].size());

      std::vector<double> splitXiExpected = mesh1EdgesToSplitExpected[edge];
      std::vector<double> splitXi;
      for (int i = 0; i < mesh1EdgesToSplit.get_num_comp(edge, 0); ++i)
        splitXi.push_back(mesh1EdgesToSplit(edge, 0, i).xi);

      std::sort(splitXi.begin(), splitXi.end());
      std::sort(splitXiExpected.begin(), splitXiExpected.end());

      for (size_t i = 0; i < splitXi.size(); ++i)
      {
        EXPECT_NEAR(splitXi[i], splitXiExpected[i], 1e-13);
      }
    }
}

} // namespace

TEST_F(MeshProjectionCalculatorTester, IdenticalQuads)
{
  std::vector<utils::Point> verts = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(1, 1), utils::Point(0, 1)};
  std::vector<std::array<int, 4>> elVerts = {{0, 1, 2, 3}};
  create_mesh1(verts, elVerts);
  create_mesh2(verts, elVerts);
  run();

  auto& verts1ToFakeVerts = *(mRelationalData->verts1ToFakeVerts);
  for (auto& vert : m_mesh1->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& verts2ToFakeVerts   = *(mRelationalData->verts2ToFakeVerts);
  auto& verts2ClassOnMesh1  = *(mRelationalData->verts2ClassOnMesh1);
  auto& edges2ToFakeVertsIn = *(mRelationalData->edges2ToFakeVertsIn);
  for (auto& vert : m_mesh2->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);

    predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert, 0, 0);
    EXPECT_EQ(record.type, PointClassification::Vert);
  }

  for (auto& edge2 : m_mesh2->get_edges())
  {
    EXPECT_EQ(edges2ToFakeVertsIn.get_num_comp(edge2, 0), 2);
    std::array<nonconformal4::impl::FakeVert, 2> edge2FakeVerts = {edges2ToFakeVertsIn(edge2, 0, 0).vert,
                                                                   edges2ToFakeVertsIn(edge2, 0, 1).vert};
    std::sort(edge2FakeVerts.begin(), edge2FakeVerts.end());

    std::array<mesh::MeshEntityPtr, 2> edge2Verts                       = {edge2->get_down(0), edge2->get_down(1)};
    std::array<nonconformal4::impl::FakeVert, 2> edge2FakeVertsExpected = {verts2ToFakeVerts(edge2Verts[0], 0, 0),
                                                                           verts2ToFakeVerts(edge2Verts[1], 0, 0)};
    std::sort(edge2FakeVertsExpected.begin(), edge2FakeVertsExpected.end());

    EXPECT_EQ(edge2FakeVerts[0].id, edge2FakeVertsExpected[0].id);
    EXPECT_EQ(edge2FakeVerts[1].id, edge2FakeVertsExpected[1].id);
  }

  auto& mesh1EdgesToSplit = *(mRelationalData->mesh1EdgesToSplit);
  for (auto& edge1 : m_mesh1->get_edges())
    EXPECT_EQ(mesh1EdgesToSplit.get_num_comp(edge1, 0), 0);

  EXPECT_EQ(get_all_fake_verts().size(), 4u);
}

TEST_F(MeshProjectionCalculatorTester, IdenticalTris)
{
  std::vector<utils::Point> verts = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1)};
  std::vector<std::array<int, 4>> elVerts = {{0, 1, 2, -1}};
  create_mesh1(verts, elVerts);
  create_mesh2(verts, elVerts);
  run();

  auto& verts1ToFakeVerts = *(mRelationalData->verts1ToFakeVerts);
  for (auto& vert : m_mesh1->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& verts2ToFakeVerts   = *(mRelationalData->verts2ToFakeVerts);
  auto& verts2ClassOnMesh1  = *(mRelationalData->verts2ClassOnMesh1);
  auto& edges2ToFakeVertsIn = *(mRelationalData->edges2ToFakeVertsIn);
  for (auto& vert : m_mesh2->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);

    predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert, 0, 0);
    EXPECT_EQ(record.type, PointClassification::Vert);
  }

  for (auto& edge2 : m_mesh2->get_edges())
  {
    EXPECT_EQ(edges2ToFakeVertsIn.get_num_comp(edge2, 0), 2);
    std::array<nonconformal4::impl::FakeVert, 2> edge2FakeVerts = {edges2ToFakeVertsIn(edge2, 0, 0).vert,
                                                                   edges2ToFakeVertsIn(edge2, 0, 1).vert};
    std::sort(edge2FakeVerts.begin(), edge2FakeVerts.end());

    std::array<mesh::MeshEntityPtr, 2> edge2Verts                       = {edge2->get_down(0), edge2->get_down(1)};
    std::array<nonconformal4::impl::FakeVert, 2> edge2FakeVertsExpected = {verts2ToFakeVerts(edge2Verts[0], 0, 0),
                                                                           verts2ToFakeVerts(edge2Verts[1], 0, 0)};
    std::sort(edge2FakeVertsExpected.begin(), edge2FakeVertsExpected.end());

    EXPECT_EQ(edge2FakeVerts[0].id, edge2FakeVertsExpected[0].id);
    EXPECT_EQ(edge2FakeVerts[1].id, edge2FakeVertsExpected[1].id);
  }

  auto& mesh1EdgesToSplit = *(mRelationalData->mesh1EdgesToSplit);
  for (auto& edge1 : m_mesh1->get_edges())
    EXPECT_EQ(mesh1EdgesToSplit.get_num_comp(edge1, 0), 0);

  EXPECT_EQ(get_all_fake_verts().size(), 3u);
}


TEST_F(MeshProjectionCalculatorTester, EdgeIntersection)
{
  std::vector<utils::Point> verts1         = {utils::Point(0, 0), utils::Point(1, 0),   utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> elVerts1 = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {3, 2, 6, 7}, {4, 0, 3, 7}};
  std::vector<utils::Point> verts2         = {utils::Point(-1, -1),  utils::Point(2, -1), utils::Point(2, 0.5),
                                              utils::Point(-1, 0.5), utils::Point(2, 2),  utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> elVerts2 = {{0, 1, 2, 3}, {3, 2, 4, 5}};
  create_mesh1(verts1, elVerts1);
  create_mesh2(verts2, elVerts2);

  run();

  auto& verts1ToFakeVerts = *(mRelationalData->verts1ToFakeVerts);
  for (auto& vert : m_mesh1->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& verts2ToFakeVerts = *(mRelationalData->verts2ToFakeVerts);
  for (auto& vert : m_mesh2->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& els                                                              = m_mesh1->get_elements();
  std::vector<predicates::impl::PointRecord> expectedVertClassifications = {
      predicates::impl::PointRecord(PointClassification::Vert, 0, els[1]),
      predicates::impl::PointRecord(PointClassification::Vert, 1, els[1]),
      predicates::impl::PointRecord(PointClassification::Edge, 0, els[2]),
      predicates::impl::PointRecord(PointClassification::Edge, 3, els[4]),
      predicates::impl::PointRecord(PointClassification::Vert, 2, els[3]),
      predicates::impl::PointRecord(PointClassification::Vert, 3, els[3])};
  test_classification(mRelationalData, m_mesh2, expectedVertClassifications);

  std::map<mesh::MeshEntityPtr, std::vector<utils::Point>> edge2FakeVertsExpected;
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, -1))]  = {utils::Point(-1, -1),
                                                                                utils::Point(2, -1)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(2, -0.25))] = {utils::Point(2, -1),
                                                                                utils::Point(2, 0.5)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 0.5))] = {
      utils::Point(-1, 0.5), utils::Point(0, 0.5), utils::Point(1, 0.5), utils::Point(2, 0.5)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(-1, -0.25))] = {utils::Point(-1, -1),
                                                                                 utils::Point(-1, 0.5)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(2, 1.25))]   = {utils::Point(2, 0.5),
                                                                                 utils::Point(2, 2)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 2))] = {utils::Point(-1, 2), utils::Point(2, 2)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(-1, 1.25))] = {utils::Point(-1, 0.5),
                                                                                utils::Point(-1, 2)};
  test_edge_fake_verts(mRelationalData, m_mesh2, edge2FakeVertsExpected);

  std::map<mesh::MeshEntityPtr, std::vector<double>> edge1SplitsExpected;
  edge1SplitsExpected[find_closest_edge(m_mesh1, utils::Point(0, 0.5))]  = {0.5};
  edge1SplitsExpected[find_closest_edge(m_mesh1, utils::Point(1, 0.5))]  = {0.5};
  edge1SplitsExpected[find_closest_edge(m_mesh1, utils::Point(-1, 0.5))] = {0.5};
  edge1SplitsExpected[find_closest_edge(m_mesh1, utils::Point(2, 0.5))]  = {0.5};
  test_edge_splits(mRelationalData, m_mesh1, edge1SplitsExpected);

  EXPECT_EQ(get_all_fake_verts().size(), 12u);
}

TEST_F(MeshProjectionCalculatorTester, EdgeIntersectionTri)
{
  std::vector<utils::Point> verts          = {utils::Point(0, 0), utils::Point(1, 0),   utils::Point(1, 1),
                                              utils::Point(0, 1)};
  std::vector<std::array<int, 4>> elVerts1 = {{0, 1, 2, -1}, {0, 2, 3, -1}};
  std::vector<std::array<int, 4>> elVerts2 = {{0, 1, 3, -1}, {1, 2, 3, -1}};
  create_mesh1(verts, elVerts1);
  create_mesh2(verts, elVerts2);

  run();

  auto& verts1ToFakeVerts = *(mRelationalData->verts1ToFakeVerts);
  for (auto& vert : m_mesh1->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& verts2ToFakeVerts = *(mRelationalData->verts2ToFakeVerts);
  for (auto& vert : m_mesh2->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& els                                                              = m_mesh1->get_elements();
  std::vector<predicates::impl::PointRecord> expectedVertClassifications = {
      predicates::impl::PointRecord(PointClassification::Vert, 0, els[0]),
      predicates::impl::PointRecord(PointClassification::Vert, 1, els[0]),
      predicates::impl::PointRecord(PointClassification::Vert, 1, els[1]),
      predicates::impl::PointRecord(PointClassification::Vert, 2, els[1])};
  test_classification(mRelationalData, m_mesh2, expectedVertClassifications);

  std::map<mesh::MeshEntityPtr, std::vector<utils::Point>> edge2FakeVertsExpected;
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 0))]   = {utils::Point(0, 0),
                                                                                utils::Point(1, 0)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(1,   0.5))] = {utils::Point(1, 0),
                                                                                utils::Point(1, 1)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 1.0))] = {utils::Point(1, 1),
                                                                                utils::Point(0, 1)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0,   0.5))] = {utils::Point(0, 1),
                                                                                utils::Point(0, 0)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 0.5))] = {utils::Point(1, 0),
                                                                                utils::Point(0.5, 0.5),
                                                                                utils::Point(0, 1)};
  test_edge_fake_verts(mRelationalData, m_mesh2, edge2FakeVertsExpected);

  std::map<mesh::MeshEntityPtr, std::vector<double>> edge1SplitsExpected;
  edge1SplitsExpected[find_closest_edge(m_mesh1, utils::Point(0.5, 0.5))]  = {0.5};

  test_edge_splits(mRelationalData, m_mesh1, edge1SplitsExpected);

  EXPECT_EQ(get_all_fake_verts().size(), 5u);
}

TEST_F(MeshProjectionCalculatorTester, EdgeCoincident)
{
  std::vector<utils::Point> verts1         = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(1, 1),
                                              utils::Point(0, 1),    utils::Point(0.5, 1.5), utils::Point(0.25, 2),
                                              utils::Point(-1, 1.5), utils::Point(-0.5, 0.5)};
  std::vector<std::array<int, 4>> elVerts1 = {{0, 1, 2, 3},  {2, 4, 3, -1}, {3, 4, 5, -1},
                                              {3, 5, 6, -1}, {3, 6, 7, -1}, {0, 3, 7, -1}};
  std::vector<utils::Point> verts2 = {verts1[0], verts1[1], verts1[2], verts1[4], verts1[7], verts1[5], verts1[6]};
  std::vector<std::array<int, 4>> elVerts2 = {{0, 1, 2, -1}, {0, 2, 3, 4}, {3, 5, 6, 4}};
  create_mesh1(verts1, elVerts1);
  create_mesh2(verts2, elVerts2);

  run();

  auto& verts1ToFakeVerts = *(mRelationalData->verts1ToFakeVerts);
  for (auto& vert : m_mesh1->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& verts2ToFakeVerts = *(mRelationalData->verts2ToFakeVerts);
  for (auto& vert : m_mesh2->get_vertices())
  {
    nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
    expect_near(fv.pt, vert->get_point_orig(0), 1e-13);
  }

  auto& els                                                              = m_mesh1->get_elements();
  std::vector<predicates::impl::PointRecord> expectedVertClassifications = {
      predicates::impl::PointRecord(PointClassification::Vert, 0, els[0]),
      predicates::impl::PointRecord(PointClassification::Vert, 1, els[0]),
      predicates::impl::PointRecord(PointClassification::Vert, 2, els[0]),
      predicates::impl::PointRecord(PointClassification::Vert, 1, els[1]),
      predicates::impl::PointRecord(PointClassification::Vert, 2, els[5]),
      predicates::impl::PointRecord(PointClassification::Vert, 2, els[2]),
      predicates::impl::PointRecord(PointClassification::Vert, 1, els[4])};
  test_classification(mRelationalData, m_mesh2, expectedVertClassifications);

  std::map<mesh::MeshEntityPtr, std::vector<utils::Point>> edge2FakeVertsExpected;
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 0))]   = {utils::Point(0, 0), utils::Point(1, 0)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(1, 0.5))]   = {utils::Point(1, 0), utils::Point(1, 1)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.5, 0.5))] = {utils::Point(0, 0), utils::Point(1, 1)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.75, 1.25))] = {utils::Point(1, 1),
                                                                                  utils::Point(0.5, 1.5)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0, 1))] = {utils::Point(-0.5, 0.5), utils::Point(0, 1),
                                                                            utils::Point(0.5, 1.5)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(-0.25, 0.25))]  = {utils::Point(0, 0),
                                                                                    utils::Point(-0.5, 0.5)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(0.375, 1.75))]  = {utils::Point(0.5, 1.5),
                                                                                    utils::Point(0.25, 2)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(-0.375, 1.75))] = {utils::Point(-1, 1.5),
                                                                                    utils::Point(0.25, 2)};
  edge2FakeVertsExpected[find_closest_edge(m_mesh2, utils::Point(-0.75, 1))]     = {utils::Point(-0.5, 0.5),
                                                                                    utils::Point(-1, 1.5)};
  test_edge_fake_verts(mRelationalData, m_mesh2, edge2FakeVertsExpected);

  std::map<mesh::MeshEntityPtr, std::vector<double>> edge1SplitsExpected;
  test_edge_splits(mRelationalData, m_mesh1, edge1SplitsExpected);

  EXPECT_EQ(get_all_fake_verts().size(), 8u);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
