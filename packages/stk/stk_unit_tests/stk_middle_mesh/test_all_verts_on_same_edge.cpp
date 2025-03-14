#include "gtest/gtest.h"

#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/element_mesh_classifier.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/predicates/point_classifier_normal_wrapper.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

using namespace predicates::impl;

class AllVertsOnSameEdgeTester : public ::testing::Test
{
  protected:
    void setup(const std::vector<utils::Point>& pts)
    {
      mesh::impl::MeshSpec spec;
      spec.numelX = 4;
      spec.numelY = 4;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;

      auto func                         = [&](const utils::Point& pt) { return utils::Point(pt.x, pt.y); };
      mesh2 = mesh::impl::create_mesh(spec, func);

      mesh1 = mesh::make_empty_mesh();
      for (auto& pt : pts)
      {
        mesh::MeshEntityPtr vert = mesh1->create_vertex(pt);
        verts.insert(vert);
      }
      vertClassifications =
          mesh::create_field<predicates::impl::PointRecord>(mesh1, mesh::impl::FieldShape(1, 0, 0), 1);

      m_pointClassifier = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(mesh2);
      for (auto& vert : mesh1->get_vertices())
      {
        for (auto mesh2El : mesh2->get_elements())
        {
          if (mesh2El)
          {
            predicates::impl::PointRecord record = m_pointClassifier->classify_reverse(mesh2El, vert->get_point_orig(0));
            if (record.type != PointClassification::Exterior)
            {
              (*vertClassifications)(vert, 0, 0) = record;
              continue;
            }
          }
        }
      }
    }

    std::shared_ptr<mesh::Mesh> mesh1;
    std::shared_ptr<mesh::Mesh> mesh2;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifier;
    mesh::FieldPtr<predicates::impl::PointRecord> vertClassifications;
    std::set<mesh::MeshEntityPtr, mesh::MeshEntityCompare> verts;
};

} // namespace

TEST_F(AllVertsOnSameEdgeTester, ZeroVerts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, OneVertInterior)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.1, 0.1)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, OneVertOnVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.25, 0.25)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, OneVertOnEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.1)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, TwoVertsOnSameEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.1), utils::Point(0.0, 0.2)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, TwoVertsOnEdgeAndVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.1), utils::Point(0.0, 0.25)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, TwoVertsOnDifferentEdges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.1), utils::Point(0.25, 0.1)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, TwoVertsOnEdgeAndDifferentVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.1), utils::Point(0.25, 0.25)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, ThreeVertsOnEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.05), utils::Point(0.0, 0.1), utils::Point(0.0, 0.2)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, ThreeVertsOnDifferentEdges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.05), utils::Point(0.0, 0.1), utils::Point(0.25, 0.2)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, ThreeVertsOnEdgeAndVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.05), utils::Point(0.0, 0.1), utils::Point(0.0, 0.25)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, ThreeVertsOnEdgeAndDifferentVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.05), utils::Point(0.0, 0.1), utils::Point(0.25, 0.25)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, ThreeVertsOnDifferentVerts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.0), utils::Point(0.0, 0.25), utils::Point(0.25, 0.25)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, FourVertsOnEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.05), utils::Point(0.0, 0.1), utils::Point(0.0, 0.15), utils::Point(0.0, 0.20)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, FourVertsOnEdgeAndVerts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.00), utils::Point(0.0, 0.1), utils::Point(0.0, 0.15), utils::Point(0.0, 0.25)});
  EXPECT_TRUE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, FourVertsOnEdgeAndDifferentEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  setup({utils::Point(0.0, 0.00), utils::Point(0.0, 0.1), utils::Point(0.0, 0.15), utils::Point(0.25, 0.2)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

TEST_F(AllVertsOnSameEdgeTester, FourVertsOnEdgeAndDifferentVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();
  setup({utils::Point(0.0, 0.0), utils::Point(0.0, 0.1), utils::Point(0.0, 0.15), utils::Point(0.25, 0.25)});
  EXPECT_FALSE(nonconformal4::impl::are_all_verts_on_same_edge(verts, vertClassifications));
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
