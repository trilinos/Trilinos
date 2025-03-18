#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/edge_tracer.hpp"
#include "stk_middle_mesh/field.hpp"
#include "gtest/gtest.h"
#include <limits>

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace stk::middle_mesh::predicates::impl;

namespace testing {


class EdgeTracerTesterBase
{
  protected:
    EdgeTracerTesterBase()
    {}

    // forward to private members of class
    std::vector<mesh::MeshEntityPtr> get_elements(mesh::MeshEntityPtr entity) { return tracer->get_elements(entity); }
    std::vector<mesh::MeshEntityPtr> get_excluded_elements(const predicates::impl::PointRecord& prevIntersection,
                                                           const predicates::impl::PointRecord& currentIntersection)
    {
      return tracer->get_excluded_elements(prevIntersection, currentIntersection);
    }

    std::vector<mesh::MeshEntityPtr> get_included_entities(std::vector<mesh::MeshEntityPtr>& allEntities,
                                                           const std::vector<mesh::MeshEntityPtr>& excludedEntities)
    {
      return tracer->get_included_entities(allEntities, excludedEntities);
    }

    bool at_end_of_line(const std::vector<mesh::MeshEntityPtr>& includedElements,
                        const predicates::impl::PointRecord& endVert)
    {
      return tracer->at_end_of_line(includedElements, endVert);
    }

    std::vector<mesh::MeshEntityPtr> get_edges(const std::vector<mesh::MeshEntityPtr>& elements)
    {
      return tracer->get_edges(elements);
    }

    std::vector<mesh::MeshEntityPtr> get_excluded_edges(const predicates::impl::PointRecord& record)
    {
      return tracer->get_excluded_edges(record);
    }

    int choose_intersection_with_minimum_angle(const std::vector<double>& alphas, const std::vector<double>& betas,
                                               const std::vector<int>& intersectionsInRange,
                                               const std::vector<mesh::MeshEntityPtr>& edges,
                                               const std::vector<int>& intersectionIdxToEdgeIdx,
                                               const utils::Point& edgeStartPt, const utils::Point& edgeEndPt)
    {
      return tracer->choose_intersection_with_minimum_angle(alphas, betas, intersectionsInRange, edges,
                                                           intersectionIdxToEdgeIdx, edgeStartPt, edgeEndPt);
    }

    nonconformal4::impl::EdgeIntersection
    create_edge_intersection(mesh::MeshEntityPtr edge, const std::vector<mesh::MeshEntityPtr>& includedElements,
                             double alpha, double beta)
    {
      return tracer->create_edge_intersection(edge, includedElements, alpha, beta);
    }

    void trace_edge(const utils::Point& ptStart, const utils::Point& normalStart,
                    const predicates::impl::PointRecord& recordStart, const utils::Point& ptEnd,
                    const utils::Point& normalEnd, const predicates::impl::PointRecord& recordEnd,
                    std::vector<nonconformal4::impl::EdgeIntersection>& intersections)
    {
      tracer->trace_edge(ptStart, normalStart, recordStart, ptEnd, normalEnd, recordEnd, intersections);
    }

    void redistribute_points_near_endpoints(std::vector<nonconformal4::impl::EdgeIntersection>& intersections,
                                            double eps)
    {
      tracer->redistribute_points_near_endpoints(intersections, eps);
    }

    std::shared_ptr<nonconformal4::impl::EdgeTracer> tracer;

    mesh::FieldPtr<utils::Point> create_normal_field(std::shared_ptr<mesh::Mesh> mesh)
    {
      auto field = mesh::create_field<utils::Point>(mesh, mesh::impl::FieldShape(1, 0, 0), 1);
      for (auto& vert : mesh->get_vertices())
        if (vert)
          (*field)(vert, 0, 0) = utils::Point(0, 0, 1);

      return field;
    }
};

} // namespace testing

namespace {
class EdgeTracerTesterQuad : public ::testing::Test,
                             public testing::EdgeTracerTesterBase
{
  protected:
      EdgeTracerTesterQuad()
      : mesh(create_test_mesh(false))
      , mesh1(create_test_mesh(false))
    {
      tracer = std::make_shared<nonconformal4::impl::EdgeTracer>(mesh, create_normal_field(mesh), middle_mesh::impl::EdgeTracerTolerances(1e-13));
    }

    std::shared_ptr<mesh::Mesh> create_test_mesh(bool createTris=false)
    {
      mesh::impl::MeshSpec spec;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;
      spec.numelX = 3;
      spec.numelY = 3;

      auto meshTest = create_mesh(
          spec, [](const utils::Point& pt) { return pt; }, MPI_COMM_SELF, createTris);
      interiorVert = meshTest->get_vertices()[5];
      assert(interiorVert->count_up() == 4);

      interiorEdge = meshTest->get_edges()[6];
      assert(interiorEdge->count_up() == 2);

      interiorEl = meshTest->get_elements()[4];

      return meshTest;
    }

    std::shared_ptr<mesh::Mesh> mesh;
    std::shared_ptr<mesh::Mesh> mesh1;
    mesh::MeshEntityPtr interiorVert;
    mesh::MeshEntityPtr interiorEdge;
    mesh::MeshEntityPtr interiorEl;
};

class EdgeTracerTesterTri : public ::testing::Test,
                            public testing::EdgeTracerTesterBase
{
  protected:
      EdgeTracerTesterTri()
      : mesh(create_test_mesh(true))
      , mesh1(create_test_mesh(true))
    {
      tracer = std::make_shared<nonconformal4::impl::EdgeTracer>(mesh, create_normal_field(mesh), middle_mesh::impl::EdgeTracerTolerances(1e-13));
    }

    std::shared_ptr<mesh::Mesh> create_test_mesh(bool createTris=false)
    {
      mesh::impl::MeshSpec spec;
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;
      spec.numelX = 3;
      spec.numelY = 3;

      auto meshTest = create_mesh(
          spec, [](const utils::Point& pt) { return pt; }, MPI_COMM_SELF, createTris);

      return meshTest;
    }

    std::shared_ptr<mesh::Mesh> mesh;
    std::shared_ptr<mesh::Mesh> mesh1;
    mesh::MeshEntityPtr interiorVert;
    mesh::MeshEntityPtr interiorEdge;
    mesh::MeshEntityPtr interiorEl;
};


template <typename T>
void expect_equal(std::vector<T> vec1, std::vector<T> vec2)
{
  std::sort(vec1.begin(), vec1.end());
  std::sort(vec2.begin(), vec2.end());
  EXPECT_EQ(vec1.size(), vec2.size());
  for (size_t i = 0; i < vec1.size(); ++i)
    EXPECT_EQ(vec1[i], vec2[i]);
}

mesh::MeshEntityPtr get_closest_element(std::shared_ptr<mesh::Mesh> mesh, const utils::Point& centroid)
{
  double minDist            = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEl = nullptr;
  for (auto& el : mesh->get_elements())
    if (el)
    {
      utils::Point disp = centroid - compute_centroid(el);
      double dist       = std::sqrt(dot(disp, disp));
      if (dist < minDist)
      {
        minDist = dist;
        minEl   = el;
      }
    }

  return minEl;
}

} // namespace

TEST_F(EdgeTracerTesterQuad, get_entity)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el = mesh->get_elements()[0];
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts, edges;
  get_downward(el, 1, edges.data());
  get_downward(el, 0, verts.data());

  for (size_t i = 0; i < verts.size(); ++i)
  {
    predicates::impl::PointRecord record(PointClassification::Vert, i, el);
    EXPECT_EQ(get_entity(record), verts[i]);
  }

  for (size_t i = 0; i < edges.size(); ++i)
  {
    predicates::impl::PointRecord record(PointClassification::Edge, i, el);
    EXPECT_EQ(get_entity(record), edges[i]);
  }
}

TEST_F(EdgeTracerTesterQuad, get_elements)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr vert = interiorVert;
  std::vector<mesh::MeshEntityPtr> els, elsExpected;
  get_upward(vert, 2, elsExpected);
  els = get_elements(vert);
  expect_equal(els, elsExpected);

  mesh::MeshEntityPtr edge = interiorEdge;
  get_upward(edge, 2, elsExpected);
  els = get_elements(edge);
  expect_equal(els, elsExpected);
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_FirstPoint)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el = mesh->get_elements()[0];
  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(),
                                   predicates::impl::PointRecord(PointClassification::Interior, -1, el));
  expect_equal(elements, {});

  elements = get_excluded_elements(predicates::impl::PointRecord(),
                                   predicates::impl::PointRecord(PointClassification::Edge, 1, el));
  expect_equal(elements, {});

  elements = get_excluded_elements(predicates::impl::PointRecord(),
                                   predicates::impl::PointRecord(PointClassification::Vert, 1, el));
  expect_equal(elements, {});
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_SecondPointFromInterior)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el = mesh->get_elements()[0];
  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Interior, -1, el),
                                   predicates::impl::PointRecord(PointClassification::Vert, 1, el));
  expect_equal(elements, {el});

  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Interior, -1, el),
                                   predicates::impl::PointRecord(PointClassification::Edge, 1, el));
  expect_equal(elements, {el});
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_TwoEdges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Edge, 0, interiorEl),
                                   predicates::impl::PointRecord(PointClassification::Edge, 2, interiorEl));
  expect_equal(elements, {interiorEl});
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_VertEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Vert, 0, interiorEl),
                                   predicates::impl::PointRecord(PointClassification::Edge, 2, interiorEl));
  expect_equal(elements, {interiorEl});
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_EdgeVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Edge, 0, interiorEl),
                                   predicates::impl::PointRecord(PointClassification::Vert, 2, interiorEl));
  expect_equal(elements, {interiorEl});
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_VertVertDiagonal)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Vert, 0, interiorEl),
                                   predicates::impl::PointRecord(PointClassification::Vert, 2, interiorEl));
  expect_equal(elements, {interiorEl});
}

TEST_F(EdgeTracerTesterQuad, get_excluded_elements_VertVertSharedEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> elements;
  elements = get_excluded_elements(predicates::impl::PointRecord(PointClassification::Vert, 0, interiorEl),
                                   predicates::impl::PointRecord(PointClassification::Vert, 1, interiorEl));
  mesh::MeshEntityPtr edge = interiorEl->get_down(0);
  expect_equal(elements, {edge->get_up(0), edge->get_up(1)});
}

TEST_F(EdgeTracerTesterQuad, get_included_entities)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> allEntities      = mesh->get_vertices();
  std::vector<mesh::MeshEntityPtr> excludedEntities = {mesh->get_vertices()[2], mesh->get_vertices()[7]};
  std::sort(excludedEntities.begin(), excludedEntities.end(), mesh::MeshEntityCompare());
  std::vector<mesh::MeshEntityPtr> includedEntities = get_included_entities(allEntities, excludedEntities);

  EXPECT_EQ(includedEntities.size(), allEntities.size() - excludedEntities.size());
  EXPECT_TRUE(std::find(includedEntities.begin(), includedEntities.end(), excludedEntities[0]) ==
              includedEntities.end());
  EXPECT_TRUE(std::find(includedEntities.begin(), includedEntities.end(), excludedEntities[1]) ==
              includedEntities.end());
}

TEST_F(EdgeTracerTesterQuad, at_end_of_line_Interior)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> includedEls = {interiorEl};
  EXPECT_TRUE(
      at_end_of_line(includedEls, predicates::impl::PointRecord(PointClassification::Interior, -1, interiorEl)));
}

TEST_F(EdgeTracerTesterQuad, at_end_of_line_EdgeSameElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> includedEls = {interiorEl};
  EXPECT_TRUE(at_end_of_line(includedEls, predicates::impl::PointRecord(PointClassification::Edge, 1, interiorEl)));
}

TEST_F(EdgeTracerTesterQuad, at_end_of_line_VertSameElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> includedEls = {interiorEl};
  EXPECT_TRUE(at_end_of_line(includedEls, predicates::impl::PointRecord(PointClassification::Vert, 1, interiorEl)));
}

TEST_F(EdgeTracerTesterQuad, at_end_of_line_EdgeDifferentElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> includedEls = {interiorEl};
  int edgeLocalId                              = 1;
  mesh::MeshEntityPtr edge                     = interiorEl->get_down(edgeLocalId);
  mesh::MeshEntityPtr otherEl                  = edge->get_up(0) == interiorEl ? edge->get_up(1) : edge->get_up(0);
  int otherElEdgeLocalId                       = -1;
  for (int i = 0; i < otherEl->count_down(); ++i)
    if (otherEl->get_down(i) == edge)
    {
      otherElEdgeLocalId = i;
      break;
    }
  assert(otherElEdgeLocalId != -1);

  EXPECT_TRUE(at_end_of_line(includedEls,
                             predicates::impl::PointRecord(PointClassification::Edge, otherElEdgeLocalId, otherEl)));
}

TEST_F(EdgeTracerTesterQuad, at_end_of_line_VertDifferentElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> includedEls = {interiorEl};
  int vertLocalId                              = 1;

  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elVerts;
  int nverts               = get_downward(interiorEl, 0, elVerts.data());
  mesh::MeshEntityPtr vert = elVerts[vertLocalId];

  std::vector<mesh::MeshEntityPtr> vertEls;
  get_upward(vert, 2, vertEls);
  mesh::MeshEntityPtr otherEl = vertEls[0] == interiorEl ? vertEls[1] : vertEls[0];

  nverts                 = get_downward(otherEl, 0, elVerts.data());
  int otherElVertLocalId = -1;
  for (int i = 0; i < nverts; ++i)
    if (elVerts[i] == vert)
    {
      otherElVertLocalId = i;
      break;
    }
  assert(otherElVertLocalId != -1);

  EXPECT_TRUE(at_end_of_line(includedEls,
                             predicates::impl::PointRecord(PointClassification::Vert, otherElVertLocalId, otherEl)));
}

TEST_F(EdgeTracerTesterQuad, get_edges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el1 = interiorEdge->get_up(0);
  mesh::MeshEntityPtr el2 = interiorEdge->get_up(1);
  std::vector<mesh::MeshEntityPtr> els{el1, el2}, edges;
  for (int i = 0; i < el1->count_down(); ++i)
    edges.push_back(el1->get_down(i));

  for (int i = 0; i < el2->count_down(); ++i)
  {
    mesh::MeshEntityPtr edge = el2->get_down(i);
    if (std::find(edges.begin(), edges.end(), edge) == edges.end())
      edges.push_back(edge);
  }

  std::vector<mesh::MeshEntityPtr> edges2 = get_edges(els);
  expect_equal(edges, edges2);
}

TEST_F(EdgeTracerTesterQuad, get_excluded_edges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<mesh::MeshEntityPtr> edges;
  edges = get_excluded_edges(predicates::impl::PointRecord(PointClassification::Interior, -1, interiorEl));
  expect_equal(edges, {});

  edges = get_excluded_edges(predicates::impl::PointRecord(PointClassification::Edge, 0, interiorEl));
  expect_equal(edges, {interiorEl->get_down(0)});

  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  get_downward(interiorEl, 0, verts.data());
  mesh::MeshEntityPtr vert = verts[0];
  std::vector<mesh::MeshEntityPtr> edgesExpected;
  for (int i = 0; i < vert->count_up(); ++i)
    edgesExpected.push_back(vert->get_up(i));

  edges = get_excluded_edges(predicates::impl::PointRecord(PointClassification::Vert, 0, interiorEl));
  expect_equal(edges, edgesExpected);
}

TEST_F(EdgeTracerTesterQuad, getBestIntersectionWhenMultipleInRange)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  auto meshTest                 = mesh::make_empty_mesh();
  mesh::MeshEntityPtr v1        = meshTest->create_vertex(0, 0, 0);
  mesh::MeshEntityPtr v2        = meshTest->create_vertex(1, 0, 0);
  mesh::MeshEntityPtr v3        = meshTest->create_vertex(1, 1, 0);
  mesh::MeshEntityPtr edge1     = meshTest->create_edge(v1, v2);
  mesh::MeshEntityPtr edge2     = meshTest->create_edge(v2, v3);

  utils::Point inputEdgeStart = utils::Point(0.5, 0.5);
  utils::Point inputEdgeEnd   = utils::Point(0.5, -0.5);

  std::vector<mesh::MeshEntityPtr> edges    = {edge1, edge2};
  std::vector<double> alphas                = {0.25, 0.5, -1, 0.75};
  std::vector<double> betas                 = {0.3, 0.5, -12, 0.8};
  std::vector<int> intersectionsInRange     = {0, 1, 3};
  std::vector<int> intersectionIdxToEdgeIdx = {0, 0, 1, 1};

  int bestIntersectionIdx = choose_intersection_with_minimum_angle(
      alphas, betas, intersectionsInRange, edges, intersectionIdxToEdgeIdx, inputEdgeStart, inputEdgeEnd);
  EXPECT_EQ(bestIntersectionIdx, 1);
}

TEST_F(EdgeTracerTesterQuad, create_edge_intersection)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el1             = interiorEl;
  mesh::MeshEntityPtr commonEdge      = interiorEl->get_down(0);
  mesh::MeshEntityPtr el2             = commonEdge->get_up(0) == el1 ? commonEdge->get_up(1) : commonEdge->get_up(0);
  mesh::MeshEntityPtr intersectedEdge = el2->get_down(1);
  double alpha = 0.5, beta = 0.1;

  nonconformal4::impl::EdgeIntersection edgeIntersection =
      create_edge_intersection(intersectedEdge, {el1, el2}, alpha, beta);
  EXPECT_EQ(edgeIntersection.record.type, PointClassification::Edge);
  EXPECT_EQ(edgeIntersection.record.id, 1);
  EXPECT_EQ(edgeIntersection.record.el, el2);
  EXPECT_EQ(edgeIntersection.alpha, alpha);
  EXPECT_EQ(edgeIntersection.beta, beta);

  beta             = 0;
  edgeIntersection = create_edge_intersection(intersectedEdge, {el1, el2}, alpha, beta);
  EXPECT_EQ(edgeIntersection.record.type, PointClassification::Vert);
  EXPECT_EQ(edgeIntersection.record.id, 1);
  EXPECT_EQ(edgeIntersection.record.el, el2);
  EXPECT_EQ(edgeIntersection.alpha, alpha);
  EXPECT_EQ(edgeIntersection.beta, -1);
}

TEST_F(EdgeTracerTesterQuad, get_entity_id_Verts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el = interiorEl;
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  int nverts = get_downward(el, 0, verts.data());

  for (int i = 0; i < nverts; ++i)
    EXPECT_EQ(predicates::impl::get_entity_id(el, verts[i]), i);
}

TEST_F(EdgeTracerTesterQuad, get_entity_id_Edges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  mesh::MeshEntityPtr el = interiorEl;
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> edges;
  int nedges = get_downward(el, 0, edges.data());

  for (int i = 0; i < nedges; ++i)
    EXPECT_EQ(predicates::impl::get_entity_id(el, edges[i]), i);
}

TEST_F(EdgeTracerTesterQuad, trace_edgeThroughEdges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point ptStart(0.25, 0.25), ptEnd(0.75, 0.25), normal(0, 0, 1);
  mesh::MeshEntityPtr elStart  = get_closest_element(mesh1, utils::Point(1.0 / 6, 1.0 / 6));
  mesh::MeshEntityPtr elMiddle = get_closest_element(mesh1, utils::Point(1.0 / 2, 1.0 / 6));
  mesh::MeshEntityPtr elEnd    = get_closest_element(mesh1, utils::Point(5.0 / 6, 1.0 / 6));
  predicates::impl::PointRecord recordStart(PointClassification::Interior, -1, elStart);
  predicates::impl::PointRecord recordEnd(PointClassification::Interior, -1, elEnd);
  std::vector<nonconformal4::impl::EdgeIntersection> intersections;

  trace_edge(ptStart, normal, recordStart, ptEnd, normal, recordEnd, intersections);

  EXPECT_EQ(intersections.size(), 2u);
  EXPECT_EQ(intersections[0].record.el, elStart);
  EXPECT_EQ(intersections[0].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[0].record.id, 1);
  EXPECT_NEAR(intersections[0].alpha, 1.0 / 6, 1e-13);
  EXPECT_NEAR(intersections[0].beta, 0.75, 1e-13);

  EXPECT_EQ(intersections[1].record.el, elMiddle);
  EXPECT_EQ(intersections[1].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[1].record.id, 1);
  EXPECT_NEAR(intersections[1].alpha, 10.0 / 12, 1e-13);
  EXPECT_NEAR(intersections[1].beta, 0.75, 1e-13);
}

TEST_F(EdgeTracerTesterTri, trace_edgeThroughEdges)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point ptStart(0.25, 0.10), ptEnd(0.75, 0.1), normal(0, 0, 1);
  mesh::MeshEntityPtr el1 = get_closest_element(mesh1, utils::Point(2.0 / 9, 1.0 / 9));
  mesh::MeshEntityPtr el2 = get_closest_element(mesh1, utils::Point(4.0 / 9, 2.0 / 9));
  mesh::MeshEntityPtr el3 = get_closest_element(mesh1, utils::Point(5.0 / 9, 1.0 / 9));
  mesh::MeshEntityPtr el4 = get_closest_element(mesh1, utils::Point(7.0 / 9, 2.0 / 9));
  predicates::impl::PointRecord recordStart(PointClassification::Interior, -1, el1);
  predicates::impl::PointRecord recordEnd(PointClassification::Interior, -1, el4);
  std::vector<nonconformal4::impl::EdgeIntersection> intersections;

  trace_edge(ptStart, normal, recordStart, ptEnd, normal, recordEnd, intersections);

  EXPECT_EQ(intersections.size(), 3u);
  EXPECT_EQ(intersections[0].record.el, el1);
  EXPECT_EQ(intersections[0].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[0].record.id, 1);
  EXPECT_NEAR(intersections[0].alpha, 1.0 / 6, 1e-13);
  EXPECT_NEAR(intersections[0].beta, 0.3, 1e-13);

  EXPECT_EQ(intersections[1].record.el, el2);
  EXPECT_EQ(intersections[1].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[1].record.id, 0);
  EXPECT_NEAR(intersections[1].alpha, 2*(0.1 + 1.0/3 - 1.0/4), 1e-13);
  EXPECT_NEAR(intersections[1].beta, 0.3, 1e-13);

  EXPECT_EQ(intersections[2].record.el, el3);
  EXPECT_EQ(intersections[2].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[2].record.id, 1);
  EXPECT_NEAR(intersections[2].alpha, 5.0/6, 1e-13);
  EXPECT_NEAR(intersections[2].beta, 0.3, 1e-13);  
}

TEST_F(EdgeTracerTesterQuad, trace_edgeThroughEdgesStartingFromEdge)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  // both endpoints are on edges
  utils::Point ptStart(0, 0.25), ptEnd(1, 0.25), normal(0, 0, 1);
  mesh::MeshEntityPtr elStart  = get_closest_element(mesh1, utils::Point(1.0 / 6, 1.0 / 6));
  mesh::MeshEntityPtr elMiddle = get_closest_element(mesh1, utils::Point(1.0 / 2, 1.0 / 6));
  mesh::MeshEntityPtr elEnd    = get_closest_element(mesh1, utils::Point(5.0 / 6, 1.0 / 6));
  predicates::impl::PointRecord recordStart(PointClassification::Edge, 3, elStart);
  predicates::impl::PointRecord recordEnd(PointClassification::Edge, 1, elEnd);
  std::vector<nonconformal4::impl::EdgeIntersection> intersections;

  trace_edge(ptStart, normal, recordStart, ptEnd, normal, recordEnd, intersections);

  EXPECT_EQ(intersections.size(), 2u);
  EXPECT_EQ(intersections[0].record.el, elStart);
  EXPECT_EQ(intersections[0].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[0].record.id, 1);
  EXPECT_NEAR(intersections[0].alpha, 1.0 / 3, 1e-13);
  EXPECT_NEAR(intersections[0].beta, 0.75, 1e-13);

  EXPECT_EQ(intersections[1].record.el, elMiddle);
  EXPECT_EQ(intersections[1].record.type, PointClassification::Edge);
  EXPECT_EQ(intersections[1].record.id, 1);
  EXPECT_NEAR(intersections[1].alpha, 2.0 / 3, 1e-13);
  EXPECT_NEAR(intersections[1].beta, 0.75, 1e-13);
}

TEST_F(EdgeTracerTesterQuad, trace_edgeThroughDiagonalVertices)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point ptStart(0.25, 0.25), ptEnd(0.75, 0.75), normal(0, 0, 1);
  mesh::MeshEntityPtr elStart  = get_closest_element(mesh1, utils::Point(1.0 / 6, 1.0 / 6));
  mesh::MeshEntityPtr elMiddle = get_closest_element(mesh1, utils::Point(1.0 / 2, 1.0 / 2));
  mesh::MeshEntityPtr elEnd    = get_closest_element(mesh1, utils::Point(5.0 / 6, 5.0 / 6));
  predicates::impl::PointRecord recordStart(PointClassification::Interior, -1, elStart);
  predicates::impl::PointRecord recordEnd(PointClassification::Interior, -1, elEnd);
  std::vector<nonconformal4::impl::EdgeIntersection> intersections;

  trace_edge(ptStart, normal, recordStart, ptEnd, normal, recordEnd, intersections);

  EXPECT_EQ(intersections.size(), 2u);
  EXPECT_EQ(intersections[0].record.el, elStart);
  EXPECT_EQ(intersections[0].record.type, PointClassification::Vert);
  EXPECT_EQ(intersections[0].record.id, 2);
  EXPECT_NEAR(intersections[0].alpha, 1.0 / 6, 1e-13);
  EXPECT_NEAR(intersections[0].beta, -1, 1e-13);

  EXPECT_EQ(intersections[1].record.el, elMiddle);
  EXPECT_EQ(intersections[1].record.type, PointClassification::Vert);
  EXPECT_EQ(intersections[1].record.id, 2);
  EXPECT_NEAR(intersections[1].alpha, 10.0 / 12, 1e-13);
  EXPECT_NEAR(intersections[1].beta, -1, 1e-13);
}

TEST_F(EdgeTracerTesterQuad, trace_edgeThroughHorizontalVertices)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point ptStart(0.25, 1.0 / 3.0), ptEnd(0.75, 1.0 / 3.0), normal(0, 0, 1);
  mesh::MeshEntityPtr elStart  = get_closest_element(mesh1, utils::Point(1.0 / 6, 1.0 / 6));
  mesh::MeshEntityPtr elMiddle = get_closest_element(mesh1, utils::Point(1.0 / 2, 1.0 / 6));
  mesh::MeshEntityPtr elEnd    = get_closest_element(mesh1, utils::Point(5.0 / 6, 1.0 / 6));
  predicates::impl::PointRecord recordStart(PointClassification::Edge, 2, elStart);
  predicates::impl::PointRecord recordEnd(PointClassification::Edge, 2, elEnd);
  std::vector<nonconformal4::impl::EdgeIntersection> intersections;

  trace_edge(ptStart, normal, recordStart, ptEnd, normal, recordEnd, intersections);

  EXPECT_EQ(intersections.size(), 2u);
  EXPECT_EQ(intersections[0].record.el, elStart);
  EXPECT_EQ(intersections[0].record.type, PointClassification::Vert);
  EXPECT_EQ(intersections[0].record.id, 2);
  EXPECT_NEAR(intersections[0].alpha, 1.0 / 6, 1e-13);
  EXPECT_NEAR(intersections[0].beta, -1, 1e-13);

  EXPECT_EQ(intersections[1].record.el, elMiddle);
  EXPECT_EQ(intersections[1].record.type, PointClassification::Vert);
  EXPECT_EQ(intersections[1].record.id, 2);
  EXPECT_NEAR(intersections[1].alpha, 10.0 / 12, 1e-13);
  EXPECT_NEAR(intersections[1].beta, -1, 1e-13);
}

TEST_F(EdgeTracerTesterQuad, trace_edgeBetweenHorizontalVertices)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  utils::Point ptStart(0, 1.0 / 3.0), ptEnd(1, 1.0 / 3.0), normal(0, 0, 1);
  mesh::MeshEntityPtr elStart  = get_closest_element(mesh1, utils::Point(1.0 / 6, 1.0 / 6));
  mesh::MeshEntityPtr elMiddle = get_closest_element(mesh1, utils::Point(1.0 / 2, 1.0 / 6));
  mesh::MeshEntityPtr elEnd    = get_closest_element(mesh1, utils::Point(5.0 / 6, 1.0 / 6));
  predicates::impl::PointRecord recordStart(PointClassification::Vert, 3, elStart);
  predicates::impl::PointRecord recordEnd(PointClassification::Vert, 2, elEnd);
  std::vector<nonconformal4::impl::EdgeIntersection> intersections;

  trace_edge(ptStart, normal, recordStart, ptEnd, normal, recordEnd, intersections);

  EXPECT_EQ(intersections.size(), 2u);
  // EXPECT_EQ(intersections[0].record.el, el_start);
  EXPECT_EQ(intersections[0].record.type, PointClassification::Vert);
  // EXPECT_EQ(intersections[0].record.id, 2);
  EXPECT_EQ(get_entity(intersections[0].record),
            get_entity(predicates::impl::PointRecord(PointClassification::Vert, 2, elStart)));
  EXPECT_NEAR(intersections[0].alpha, 1.0 / 3, 1e-13);
  EXPECT_NEAR(intersections[0].beta, -1, 1e-13);

  // EXPECT_EQ(intersections[1].record.el, el_middle);
  EXPECT_EQ(intersections[1].record.type, PointClassification::Vert);
  // EXPECT_EQ(intersections[1].record.id, 2);
  EXPECT_EQ(get_entity(intersections[1].record),
            get_entity(predicates::impl::PointRecord(PointClassification::Vert, 2, elMiddle)));

  EXPECT_NEAR(intersections[1].alpha, 2.0 / 3, 1e-13);
  EXPECT_NEAR(intersections[1].beta, -1, 1e-13);
}

// Test cases:
//   1. If no points within m_eps, do nothing
//   2. Single point withing range
//   3. 3 points within range
TEST_F(EdgeTracerTesterQuad, redistribute_points_near_endpointNoChange)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<nonconformal4::impl::EdgeIntersection> intersections;
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 0.5, 0.5});

  redistribute_points_near_endpoints(intersections, 0.1);

  EXPECT_EQ(intersections[0].alpha, 0.5);
  EXPECT_EQ(intersections[0].beta, 0.5);
}

TEST_F(EdgeTracerTesterQuad, redistribute_points_near_endpointSinglePointLeft)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<nonconformal4::impl::EdgeIntersection> intersections;
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 0, 0.5});
  double eps = 0.1;

  redistribute_points_near_endpoints(intersections, eps);

  EXPECT_EQ(intersections[0].alpha, eps / 2);
  EXPECT_EQ(intersections[0].beta, 0.5);
}

TEST_F(EdgeTracerTesterQuad, redistribute_points_near_endpointSinglePointRight)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<nonconformal4::impl::EdgeIntersection> intersections;
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 1, 0.5});
  double eps = 0.1;

  redistribute_points_near_endpoints(intersections, eps);

  EXPECT_EQ(intersections[0].alpha, 1 - eps / 2);
  EXPECT_EQ(intersections[0].beta, 0.5);
}

TEST_F(EdgeTracerTesterQuad, redistribute_points_near_endpointThreePointsLeft)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<nonconformal4::impl::EdgeIntersection> intersections;
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 0, 0.5});
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 0, 0.5});
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 0, 0.5});

  double eps = 0.1;

  redistribute_points_near_endpoints(intersections, eps);

  EXPECT_EQ(intersections[0].alpha, eps / 4);
  EXPECT_EQ(intersections[0].beta, 0.5);

  EXPECT_EQ(intersections[1].alpha, eps / 2);
  EXPECT_EQ(intersections[1].beta, 0.5);

  EXPECT_EQ(intersections[2].alpha, 3 * eps / 4);
  EXPECT_EQ(intersections[2].beta, 0.5);
}

TEST_F(EdgeTracerTesterQuad, redistribute_points_near_endpointThreePointsRight)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<nonconformal4::impl::EdgeIntersection> intersections;
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 1, 0.5});
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 1, 0.5});
  intersections.push_back(nonconformal4::impl::EdgeIntersection{predicates::impl::PointRecord(), 1, 0.5});

  double eps = 0.1;

  redistribute_points_near_endpoints(intersections, eps);

  EXPECT_EQ(intersections[0].alpha, 1 - 3 * eps / 4);
  EXPECT_EQ(intersections[0].beta, 0.5);

  EXPECT_EQ(intersections[1].alpha, 1 - eps / 2);
  EXPECT_EQ(intersections[1].beta, 0.5);

  EXPECT_EQ(intersections[2].alpha, 1 - eps / 4);
  EXPECT_EQ(intersections[2].beta, 0.5);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
