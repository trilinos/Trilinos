
#include "gtest/gtest.h"
#include "stk_middle_mesh/mesh_relational_data_scatter.hpp"
#include "stk_middle_mesh/predicates/point_classifier_normal_wrapper.hpp"
#include "stk_middle_mesh/mesh_scatter.hpp"
#include "stk_middle_mesh/mesh_scatter_spec.hpp"
#include "stk_middle_mesh/bounding_box_search.hpp"
#include "stk_middle_mesh/middle_grid_constraint_generator.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_projection_calculator.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/middle_grid_triangulator.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "util/parallel_search_test_util.hpp"
#include "util/nonconformal_interface_helpers.hpp"
#include "stk_middle_mesh/variable_size_field.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/mesh_relational_data.hpp"
#include "stk_middle_mesh/edge_tracer_opts.hpp"
#include "stk_middle_mesh/bounding_box_search.hpp"

namespace {

using namespace stk::middle_mesh;

class MiddleGridTriangulatorParallelTester : public ::testing::Test
{
  protected:

    using BoundingBoxSearch = search::ElementToElementBoundingBoxSearch;

    void setup_spmd(const mesh::impl::MeshSpec& spec1, const mesh::impl::MeshSpec& /*spec2*/, bool createTriangles=false)
    {
      auto f = [](const utils::Point& pt) { return pt; };

      auto mesh1 = mesh::impl::create_mesh(spec1, f, MPI_COMM_WORLD, createTriangles);
      auto mesh2 = mesh::impl::create_mesh(spec1, f, MPI_COMM_WORLD, createTriangles);

      setup(mesh1, mesh2);
    }

    void setup_mpmd(const mesh::impl::MeshSpec& spec1, const mesh::impl::MeshSpec& spec2,
                    const SplitCommTestUtil& splitter, bool createTriangles=false)
    {
      auto f = [](const utils::Point& pt) { return pt; };

      std::shared_ptr<mesh::Mesh> mesh1, mesh2;

      if (splitter.get_color() == SplitCommColor::SEND)
      {
        mesh1 = mesh::impl::create_mesh(spec1, f, splitter.get_comm(), createTriangles);
      }

      if (splitter.get_color() == SplitCommColor::RECV)
      {
        mesh2 = mesh::impl::create_mesh(spec2, f, splitter.get_comm(), createTriangles);
      }

      setup(mesh1, mesh2);      
    }

    void setup_mpmd(const std::vector<utils::Point>& vertCoords1,
                    const std::vector<std::array<int, 4>>& elementVerts1,
                    const std::vector<utils::Point>& vertCoords2,
                    const std::vector<std::array<int, 4>>& elementVerts2)
    {
      assert(utils::impl::comm_size(MPI_COMM_WORLD) == 2);
      std::shared_ptr<mesh::Mesh> mesh1, mesh2;
      if (utils::impl::comm_rank(MPI_COMM_WORLD) == 0)
      {
        mesh1 = mesh::make_empty_mesh(MPI_COMM_SELF);
        create_mesh(mesh1, vertCoords1, elementVerts1);
      } else
      {
        mesh2 = mesh::make_empty_mesh(MPI_COMM_SELF);
        create_mesh(mesh2, vertCoords2, elementVerts2);
      }

      setup(mesh1, mesh2);
    }

    void setup(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2)
    {
      m_mesh1 = mesh1;
      m_mesh2 = mesh2;

      auto scatterSpec = create_scatter_spec(mesh1, mesh2);
      scatter_mesh_1to2(scatterSpec);
      if (m_mesh2)
      {
        m_meshRelationalDataOnMesh2 = std::make_shared<nonconformal4::impl::MeshRelationalData>(m_mesh1ScatteredToMesh2, m_mesh2);
      }

      do_mesh_projections();

      scatterSpec = create_scatter_spec(mesh2, mesh1);
      scatter_mesh_2to1(scatterSpec);
      scatter_mesh_relational_data_2to1();

      create_middle_mesh_verts_and_edges();
      create_middle_mesh_triangles();
    }

    void run_tests()
    {
      if (m_mesh1)
      {
        m_meshInElementsToMesh1                      = m_meshRelationalDataOnMesh1->meshInElementsToMesh1Elements;
        m_meshInElementsToMesh2ScatteredToMesh1      = m_meshRelationalDataOnMesh1->meshInElementsToMesh2Elements;
        m_mesh1InverseClassification                 = nonconformal4::impl::invert_classification_field(
                                                        m_meshInOnMesh1Procs, m_mesh1, m_meshInElementsToMesh1);
        m_mesh2ScatteredToMesh1InverseClassification = nonconformal4::impl::invert_classification_field(
                                                        m_meshInOnMesh1Procs, m_mesh2ScatteredToMesh1,
                                                        m_meshInElementsToMesh2ScatteredToMesh1);

        stk::middle_mesh::test_util::test_areas_positive(m_meshInOnMesh1Procs);
        stk::middle_mesh::test_util::test_every_element_classified(m_meshInOnMesh1Procs, m_meshInElementsToMesh1);
        stk::middle_mesh::test_util::test_every_element_classified(m_meshInOnMesh1Procs, m_meshInElementsToMesh2ScatteredToMesh1);
        stk::middle_mesh::test_util::test_every_element_classified_inverse(m_mesh1, m_mesh1InverseClassification);
        stk::middle_mesh::test_util::test_area_per_element(m_mesh1, m_mesh1InverseClassification);

        mesh::check_topology(m_meshInOnMesh1Procs);
      }

      stk::middle_mesh::test_util::test_total_areas_same(m_mesh1, m_mesh2, m_meshInOnMesh1Procs);
      test_mesh2_area_per_element(m_mesh2ScatteredToMesh1InverseClassification);

    }


    std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec(
        std::shared_ptr<mesh::Mesh> sendMesh, std::shared_ptr<mesh::Mesh> recvMesh)
    {
      std::shared_ptr<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBox> sendMeshSearchAdaptor, recvMeshSearchAdaptor;
      if (sendMesh)
      {
        sendMeshSearchAdaptor = std::make_shared<mesh::impl::SearchMeshElementBoundingBox>(sendMesh, MPI_COMM_WORLD);
      }

      if (recvMesh)
      {
        recvMeshSearchAdaptor = std::make_shared<mesh::impl::SearchMeshElementBoundingBox>(recvMesh, MPI_COMM_WORLD);
      }

      auto search = std::make_shared<search::ElementToElementBoundingBoxSearch>(sendMeshSearchAdaptor, recvMeshSearchAdaptor, 
                                                                                    "mesh1To2Search", MPI_COMM_WORLD);

      search->coarse_search();

      auto scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, sendMesh);
      if (sendMesh)
      {

        const std::vector<BoundingBoxSearch::BoundingBoxB>&unpairedEntities = search->get_unpaired_recv_entities();
        BoundingBoxSearch::EntityProcRelationVec& mesh2To1Relations = search->get_range_to_domain();

        if (unpairedEntities.size() > 0)
          throw std::runtime_error("coarse search could not find a destination for some mesh1 entities");

        for (const BoundingBoxSearch::EntityProcRelation& searchRelation : mesh2To1Relations)
        {
          int destProc = searchRelation.first.proc();
          int mesh1ElementLocalId = searchRelation.second.id();

          mesh::MeshEntityPtr mesh1El = sendMesh->get_elements()[mesh1ElementLocalId];
          scatterSpec->add_destination(mesh1El, destProc);
        }
      }

      return scatterSpec;

    }

    void scatter_mesh_1to2(std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec)
    {
      MPI_Comm scatteredMeshComm = m_mesh2 ? m_mesh2->get_comm() : MPI_COMM_NULL;
      mesh::impl::MeshScatter scatterer(scatterSpec, m_mesh1, scatteredMeshComm, true);
      m_mesh1ScatteredToMesh2   = scatterer.scatter();
      m_mesh1EntityOrigins      = scatterer.get_entity_origins();
      m_mesh1EntityDestinations = scatterer.get_entity_destinations();
    }

    void scatter_mesh_2to1(std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec)
    {
      MPI_Comm scatteredMeshComm = m_mesh1 ? m_mesh1->get_comm() : MPI_COMM_NULL;
      mesh::impl::MeshScatter scatterer(scatterSpec, m_mesh2, scatteredMeshComm, true);
      m_mesh2ScatteredToMesh1   = scatterer.scatter();
      m_mesh2EntityOrigins      = scatterer.get_entity_origins();
      m_mesh2EntityDestinations = scatterer.get_entity_destinations();
    }


    void do_mesh_projections()
    {
      if (m_mesh2)
      {
        m_pointClassifierOnMesh2 = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(m_mesh2);
        nonconformal4::impl::MeshProjectionCalculator meshProjection(m_mesh1ScatteredToMesh2, m_mesh2, m_meshRelationalDataOnMesh2, m_pointClassifierOnMesh2,
                                                                     stk::middle_mesh::impl::EdgeTracerTolerances());
        meshProjection.project();
      }
    }

    void scatter_mesh_relational_data_2to1()
    {
      if (m_mesh2ScatteredToMesh1)
      {
        m_pointClassifierOnMesh1 = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(m_mesh2ScatteredToMesh1);
      }

      nonconformal4::impl::MeshRelationalDataScatterInput scatterInput;
      scatterInput.mesh1 = m_mesh1;
      scatterInput.mesh2 = m_mesh2;
      scatterInput.mesh1ScatteredToMesh2 = m_mesh1ScatteredToMesh2;
      scatterInput.mesh2ScatteredToMesh1 = m_mesh2ScatteredToMesh1;
      scatterInput.mesh1EntityOrigins = m_mesh1EntityOrigins;
      scatterInput.mesh1EntityDestinations = m_mesh1EntityDestinations;
      scatterInput.mesh2EntityOrigins = m_mesh2EntityOrigins;
      scatterInput.mesh2EntityDestinations = m_mesh2EntityDestinations;

      nonconformal4::impl::MeshRelationalDataScatter scatterer(scatterInput, m_meshRelationalDataOnMesh2,
                                                               m_pointClassifierOnMesh2,
                                                               m_pointClassifierOnMesh1,
                                                               MPI_COMM_WORLD);
      m_meshRelationalDataOnMesh1 = scatterer.scatter();
      m_meshInOnMesh1Procs        = scatterer.get_middle_mesh();
    }

    void create_middle_mesh_verts_and_edges()
    {
      if (m_mesh1)
      {
        nonconformal4::impl::MiddleGridConstraintGenerator generator(
          m_mesh1, m_mesh2ScatteredToMesh1, m_meshInOnMesh1Procs, m_meshRelationalDataOnMesh1, m_pointClassifierOnMesh1);

        generator.generate();
      }
    }

    void create_middle_mesh_triangles()
    {
      if (m_mesh1)
      {
        nonconformal4::impl::MiddleGridTriangulator triangulator(m_mesh1, m_mesh2ScatteredToMesh1,
                               m_meshInOnMesh1Procs, m_meshRelationalDataOnMesh1, m_pointClassifierOnMesh1);

        triangulator.triangulate();
      }
    }

    void test_mesh2_area_per_element(mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> mesh2ScatteredToMesh1InverseClassificationPtr)
    {
      stk::DataExchangeUnknownPatternNonBlockingCommBuffer exchanger(MPI_COMM_WORLD);

      if (m_mesh2ScatteredToMesh1)
      {
        auto& mesh2EntityOrigins = *m_mesh2EntityOrigins;
        auto& invClassification = *mesh2ScatteredToMesh1InverseClassificationPtr;

        mesh::impl::ElementOperations2D elemOps;
        for (int phase=0; phase < 2; ++phase)
        {
          for (auto& el2 : m_mesh2ScatteredToMesh1->get_elements())
            if (el2)
            {
              double area = 0;
              for (auto& elIn : invClassification(el2, 0))
              {
                area += elemOps.compute_area(elIn);
              }

              mesh::RemoteSharedEntity remote = mesh2EntityOrigins(el2, 0, 0);
              exchanger.get_send_buf(remote.remoteRank).pack(remote.remoteId);
              exchanger.get_send_buf(remote.remoteRank).pack(area);
            }

          if (phase == 0)
            exchanger.allocate_send_buffers();
        }
      }

      exchanger.start_nonblocking();
      exchanger.post_nonblocking_receives();

      mesh::FieldPtr<double> elAreasPtr;
      if (m_mesh2)
      {
        elAreasPtr = mesh::create_field<double>(m_mesh2, mesh::FieldShape(0, 0, 1), 1, 0);
      }

      auto unpacker = [&](int /*rank*/, stk::CommBuffer& buf)
      {
        auto& elAreas = *elAreasPtr;
        while (buf.remaining() > 0)
        {
          int elId;
          double area;
          buf.unpack(elId);
          buf.unpack(area);
          mesh::MeshEntityPtr el = m_mesh2->get_elements()[elId];
          elAreas(el, 0, 0) += area;
        }
      };

      exchanger.complete_receives(unpacker);

      if (m_mesh2)
      {
        auto& elAreas = *elAreasPtr;
        mesh::impl::ElementOperations2D elemOps;
        for (auto& el : m_mesh2->get_elements())
          if (el)
          {
            double areaEx = elemOps.compute_area(el);
            EXPECT_NEAR(areaEx, elAreas(el, 0, 0), 1e-13);
          }
      }
    }

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

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;

    std::shared_ptr<mesh::Mesh> m_mesh1ScatteredToMesh2;
    std::shared_ptr<mesh::Mesh> m_mesh2ScatteredToMesh1;

    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh1EntityOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh1EntityDestinations;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityDestinations;

    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierOnMesh2;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierOnMesh1;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> m_meshRelationalDataOnMesh2;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> m_meshRelationalDataOnMesh1;
    std::shared_ptr<mesh::Mesh> m_meshInOnMesh1Procs;

    mesh::FieldPtr<mesh::MeshEntityPtr> m_meshInElementsToMesh1;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_meshInElementsToMesh2ScatteredToMesh1;    
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_mesh1InverseClassification;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_mesh2ScatteredToMesh1InverseClassification;
};


}  // namespace


TEST_F(MiddleGridTriangulatorParallelTester, SingleEdgeOverlap)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // The top edge of el1 overlaps with the top edge of el2.  No other
  // edges itnersect
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0),         utils::Point(1, 1),
                                              utils::Point(0, 1),    utils::Point(-0.25, -0.25), utils::Point(1.25, -0.25),
                                              utils::Point(1.25, 1), utils::Point(-0.25, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {0, 3, 7, 4}};

  std::vector<utils::Point> mesh2Verts = {utils::Point(-0.25, -0.25), utils::Point(1.25, -0.25), utils::Point(1.25, 1),
                                          utils::Point(-0.25, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  
  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 8);
}

TEST_F(MiddleGridTriangulatorParallelTester, SingleEdgeOverlapTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // The left edge of el1 overlaps with the top edge of el2.  No other
  // edges itnersect
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0), utils::Point(0, 1),    
                                              utils::Point(0.25, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 3, 2, -1}, {0, 1, 3, -1}, {3, 1, 2, -1}};

  std::vector<utils::Point> mesh2Verts = {mesh1Verts[0], mesh1Verts[1], mesh1Verts[2]};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 3);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapBisection)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // One edge of mesh2 cuts el1 in half.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0),   utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}, {3, 2, 6, 7}, {4, 0, 3, 7}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1),  utils::Point(2, -1), utils::Point(2, 0.5),
                                              utils::Point(-1, 0.5), utils::Point(2, 2),  utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {3, 2, 4, 5}};
  
  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 16);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapBisectionTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }  

  // One edge of mesh2 cuts el1 in half.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0), utils::Point(1, 0), utils::Point(0, 1), utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {0, 3, 2, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 2);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapSharedCorner)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // The top left corner of el1 overlaps the top left corner of el2, but el2
  // is larger than el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, -1), utils::Point(1, -1), utils::Point(2, -1),
                                              utils::Point(0, 0),  utils::Point(1, 0),  utils::Point(2, 0),
                                              utils::Point(0, 1),  utils::Point(1, 1),  utils::Point(2, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 4, 3}, {1, 2, 5, 4}, {3, 4, 7, 6}, {4, 5, 8, 7}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, -1), utils::Point(2, -1), utils::Point(2, 1),
                                              utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}};
  
  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 8);
}


TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapCutCorner)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // An edge of the first element on mesh2 cuts off the top left corner of el1.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(1, 1),
                                              utils::Point(0, 1),    utils::Point(0.75, 2),  utils::Point(1, 3),
                                              utils::Point(-0.5, 3), utils::Point(-0.5, -1), utils::Point(2, -1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},  {7, 8, 1, 0},  {1, 8, 2, -1}, {2, 8, 4, -1},
                                              {2, 4, 3, -1}, {3, 4, 5, -1}, {3, 5, 6, -1}, {0, 3, 6, 7}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-0.5, 0), utils::Point(0.75, 2),  utils::Point(1, 3),
                                              utils::Point(-0.5, 3), utils::Point(-0.5, -1), utils::Point(2, -1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {4, 5, 1, 0}};
  
  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 16);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapCutCornerTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // An edge of the first element on mesh2 cuts off the top left corner of el1.
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.75, 0.75),  utils::Point(-0.25, 0.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1},  {1, 3, 2, -1},  {0, 2, 4, -1}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0),    utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.75, 0.75),  utils::Point(-0.25, 0.75)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {0, 3, 4, -1}, {4, 3, 2, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 9);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeEdgeOverlapCutCorners)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // edges of an element on mesh2 cut off the top left and top right corners
  // of el1.  No other edges are intersected
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0),     utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(-0.5, -1), utils::Point(1.25, -0.5),
                                              utils::Point(2, 3), utils::Point(-0.5, 3)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7}, {0, 3, 7, 4}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-0.5, 0), utils::Point(0.75, 2),  utils::Point(1, 3),
                                              utils::Point(-0.5, 3), utils::Point(-0.5, -1), utils::Point(1.25, -.5),
                                              utils::Point(2, 3)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 2, 1}};
  
  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 22);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeEdgeOverlapCutCornersTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // edges of an element on mesh2 cut off the top left and top right corners
  // of el1.  No other edges are intersected
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),     utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.2, 0.9), utils::Point(0.9, 0.2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, -1}, {1, 4, 2, -1}, {4, 3, 2, -1}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0),     utils::Point(1, 0),     utils::Point(0, 1),
                                              utils::Point(0.2, 0.9), utils::Point(0.9, 0.2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, -1}, {0, 4, 3, -1}, {0, 3, 2, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 9);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapBisection2)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts     = {utils::Point(-1, 0), utils::Point(0, 0),  utils::Point(1, 0),
                                              utils::Point(3, 0),  utils::Point(-1, 1), utils::Point(0, 1),
                                              utils::Point(1, 1),  utils::Point(3, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, 0), utils::Point(0.5, 0), utils::Point(0.5, 1),
                                              utils::Point(-1, 1), utils::Point(3, 0),   utils::Point(3, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 8);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoEdgeOverlapBisection2Tri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0), utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, 0), utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1), utils::Point(0.5, 0), utils::Point(0.5, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 4, 3, -1}, {4, 5, 3, -1}, {4, 2, 5, -1}, {4, 1, 2, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 4);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoVertexOverlapDiagonal)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),  utils::Point(1, 0), utils::Point(1, 1),   utils::Point(0, 1),  utils::Point(-1, -1),
      utils::Point(2, -1), utils::Point(3, 0), utils::Point(3, 1.5), utils::Point(-1, 2), utils::Point(-2, 0.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, -1}, {1, 6, 7, 2},
                                              {2, 7, 8, 3}, {0, 3, 8, 9}, {4, 0, 9, -1}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1),  utils::Point(2, -1), utils::Point(-1, 2),
                                              utils::Point(-2, 0.5), utils::Point(3, 0),  utils::Point(3, 1.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 12);
}

TEST_F(MiddleGridTriangulatorParallelTester, TwoVertexOverlapDiagonalTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // One edge of mesh2 cuts el1 in half.  The top edges of mesh2 elements overlap
  // the top edge of el1, and similarly for the bottom edges
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),   utils::Point(1, 0),  utils::Point(1, 1), utils::Point(0, 1),
      utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7},
                                              {4, 0, 3, 7}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 3, -1}, {1, 2, 3, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 10);
}

TEST_F(MiddleGridTriangulatorParallelTester, FourVertexOverlapDiagonals)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // mesh2 edges cross both diagonals of element 1, dividing it into 4 triangles
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(0.5, -3), utils::Point(2, -1),
                                              utils::Point(3, 0.5), utils::Point(2, 2),    utils::Point(0.5, 3),
                                              utils::Point(-1, 2),  utils::Point(-3, 0.5), utils::Point(-1, -1)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0.5, 0.5), utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(2, 2),     utils::Point(-1, 2),  utils::Point(0.5, -3),
                                              utils::Point(3, 0.5),   utils::Point(0.5, 3), utils::Point(-3, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 16);
}

TEST_F(MiddleGridTriangulatorParallelTester, FourVertexOverlapDiagonalsTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // mesh2 edges cross both diagonals of element 1, dividing it into 4 triangles
  std::vector<utils::Point> mesh1Verts = {
      utils::Point(0, 0),   utils::Point(1, 0),  utils::Point(1, 1), utils::Point(0, 1),
      utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, 2}, {3, 2, 6, 7},
                                              {4, 0, 3, 7}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(2, -1), utils::Point(2, 2), utils::Point(-1, 2),
                                              utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 4, 3, -1}, {0, 1, 4, -1}, {1, 2, 4, -1}, {2, 3, 4, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 12);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeVertexOverlapDiagonals)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(0.5, -3), utils::Point(2, -1),
                                              utils::Point(3, 0.5), utils::Point(2, 2),    utils::Point(0.5, 3),
                                              utils::Point(-1, 2),  utils::Point(-3, 0.5), utils::Point(-1, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0.5, 0.5), utils::Point(-1, 0.25), utils::Point(2, -1),
                                              utils::Point(2, 2),     utils::Point(-1, 2),    utils::Point(0.5, -3),
                                              utils::Point(3, 0.5),   utils::Point(0.5, 3),   utils::Point(-3, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 18);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeVertexOverlapDiagonalsTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(-1, -1), utils::Point(2, -1),
                                              utils::Point(-1, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 2, -1},
                                              {3, 2, 6, -1}, {4, 0, 3, 6}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(2, -1), utils::Point(1, 1),
                                              utils::Point(-1, 2),  utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, -1}, {1, 2, 4, -1}, {4, 2, 3, -1}, {0, 4, 3, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 10);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeVertexNonCentroid)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // mesh2 edges cross 3 vertices and one edge of el1, but the interior point is not
  // on the centroid of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(0.5, -3), utils::Point(3, -2.0 / 3.0),
                                              utils::Point(3, 0.5), utils::Point(2, 2),    utils::Point(0.5, 3.5),
                                              utils::Point(-1, 4),  utils::Point(-3, 2),   utils::Point(-1, 0.125)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 0, 11, -1}, {4, 5, 1, 0},
                                              {5, 6, 2, 1}, {6, 7, 2, -1},  {2, 7, 8, -1},
                                              {2, 8, 9, 3}, {3, 9, 10, -1}, {0, 3, 10, 11}};

  std::vector<utils::Point> mesh2Verts = {
      utils::Point(0.25, 0.25), utils::Point(-1, 0.125), utils::Point(3, -2.0 / 3.0),
      utils::Point(2, 2),       utils::Point(-1, 4),     utils::Point(0.5, -3),
      utils::Point(3, 0.5),     utils::Point(0.5, 3.5),  utils::Point(-3, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{1, 5, 2, 0}, {2, 6, 3, 0}, {0, 3, 7, 4}, {1, 0, 4, 8}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 18);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeVertexNonCentroidTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // mesh2 edges cross 3 vertices and one edge of el1
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),   utils::Point(1, 0),    utils::Point(1, 1),
                                              utils::Point(0, 1),   utils::Point(-1, -1), utils::Point(1.75, -0.25),
                                              utils::Point(-0.25, 1.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 2, -1},
                                              {3, 2, 6, -1}, {4, 0, 3, 6}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(-1, -1), utils::Point(1.75, -0.25), utils::Point(1, 1),
                                              utils::Point(-0.25, 1.75),  utils::Point(0.25, 0.25)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, -1}, {1, 2, 4, -1}, {4, 2, 3, -1}, {0, 4, 3, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 10);
}

TEST_F(MiddleGridTriangulatorParallelTester, DoubleEdgeIntersection)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // Two edges of a mesh2 element intersect the same edge of the mesh1 element
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),     utils::Point(1, 0),       utils::Point(1, 1),
                                              utils::Point(0, 1),     utils::Point(-0.25, -1),  utils::Point(1.5, -1),
                                              utils::Point(1.5, 2.0), utils::Point(0.75, 2),    utils::Point(0.25, 2.5),
                                              utils::Point(-1, 0.75), utils::Point(-1.5, 0.25), utils::Point(-1, 0.25)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},   {4, 5, 1, 0},  {5, 6, 2, 1},   {2, 6, 7, -1},
                                              {2, 7, 8, 3},   {3, 8, 9, -1}, {3, 9, 11, -1}, {9, 10, 11, -1},
                                              {0, 3, 11, -1}, {0, 11, 4, -1}};

  std::vector<utils::Point> mesh2Verts = {
      utils::Point(0.5, 0.5),  utils::Point(-0.25, -1), utils::Point(0.5, -1),    utils::Point(-1, 0.25),
      utils::Point(1.5, -1),   utils::Point(1.5, 0.5),  utils::Point(1.5, 2),     utils::Point(0.75, 2),
      utils::Point(0.25, 2.5), utils::Point(-1, 0.75),  utils::Point(-1.5, 0.25),
  };
  std::vector<std::array<int, 4>> mesh2Els = {{1, 2, 0, 3}, {2, 4, 5, 0}, {0, 5, 6, 7}, {0, 7, 8, 9}, {3, 0, 9, 10}};
  
  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 30);
}

TEST_F(MiddleGridTriangulatorParallelTester, DoubleEdgeIntersectionTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // Two edges of a mesh2 element intersect the same edge of the mesh1 element
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),     utils::Point(1, 0),         utils::Point(1, 1),
                                              utils::Point(0, 1),     utils::Point(-0.25, 0.25),  utils::Point(-0.25, 0.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3},   {4, 0, 3, 5}};

  std::vector<utils::Point> mesh2Verts = {utils::Point(0, 0),     utils::Point(1, 0),         utils::Point(1, 1),
                                          utils::Point(0, 1),     utils::Point(-0.25, 0.25),  utils::Point(-0.25, 0.75),
                                          utils::Point(0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 6, 4, -1}, {0, 1, 6, -1}, {1, 2, 6, -1}, {6, 2, 3, -1}, {5, 6, 3, -1}, {4, 6, 5, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 10);
}

TEST_F(MiddleGridTriangulatorParallelTester, EdgeCrossThroughCornerVertex)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // Edge starts on midpoint of one edge of el1 and passes through a vertex that
  // is not adjacent to the edge
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),  utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1),  utils::Point(0, -1), utils::Point(2, -2),
                                              utils::Point(3, -1), utils::Point(3, 2),  utils::Point(0, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {5, 6, 1, -1}, {1, 6, 7, 2}, {2, 7, 8, 3}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, -1), utils::Point(2, -2), utils::Point(0.5, 1),
                                              utils::Point(0, 1),  utils::Point(3, -1), utils::Point(3, 2),
                                              utils::Point(0, 2)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 4, 5, 2}, {3, 2, 5, 6}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 11);
}

TEST_F(MiddleGridTriangulatorParallelTester, EdgeCrossThroughCornerVertexTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }  
  // Edge starts on midpoint of one edge of el1 and passes through a vertex that
  // is not adjacent to the edge
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),  utils::Point(1, 0),  utils::Point(1, 1),
                                              utils::Point(0, 1),  utils::Point(0, -1), utils::Point(1.5, -1),
                                              utils::Point(1.5, 2)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0}, {1, 5, 6, -1}, {1, 6, 2, -1}, {3, 2, 6, -1}};

  std::vector<utils::Point> mesh2Verts     = {utils::Point(0, -1), utils::Point(1.5, -1), utils::Point(1.5, 2),
                                              utils::Point(0, 1),  utils::Point(0.5, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 4, 3}, {1, 2, 4, -1}, {3, 4, 2, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 9);
}

TEST_F(MiddleGridTriangulatorParallelTester, ThreeElementsWithCutCorner)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // two elements cover most of el1, but a third element cuts off a corner
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),       utils::Point(1, 0),      utils::Point(1, 1),
                                              utils::Point(0, 1),       utils::Point(-1, 0),     utils::Point(0.5, -0.25),
                                              utils::Point(1.5, -0.25), utils::Point(1.5, 1.25), utils::Point(-0.25, 1.25),
                                              utils::Point(-0.25, 0.5), utils::Point(-0.5, 0.5)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {5, 6, 1, 0},  {1, 6, 7, 2}, {2, 7, 8, 3},
                                              {0, 3, 8, 9}, {0, 9, 10, 4}, {0, 4, 5, -1}};

  std::vector<utils::Point> mesh2Verts = {utils::Point(0.5, -0.25), utils::Point(-0.25, 0.5), utils::Point(-0.5, 0.5),
                                          utils::Point(-1, 0),      utils::Point(1.5, -0.25), utils::Point(1.5, 1.25),
                                          utils::Point(0.5, 1.25),  utils::Point(-0.25, 1.25)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {0, 4, 5, 6}, {0, 6, 7, 1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 22);
}


TEST_F(MiddleGridTriangulatorParallelTester, ThreeElementsWithCutCornerTri)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  // two elements cover most of el1, but a third element cuts off a corner
  std::vector<utils::Point> mesh1Verts     = {utils::Point(0, 0),       utils::Point(1, 0),      utils::Point(1, 1),
                                              utils::Point(0, 1),       utils::Point(-1, -1),    utils::Point(0.75, -0.5),
                                              utils::Point(-0.5, 0.75)};
  std::vector<std::array<int, 4>> mesh1Els = {{0, 1, 2, 3}, {4, 5, 1, 0},  {4, 0, 3, 6}};

  std::vector<utils::Point> mesh2Verts = {utils::Point(0.75, -0.5), utils::Point(1, 0),   utils::Point(0, 1),
                                          utils::Point(-0.5, 0.75), utils::Point(-1, -1), utils::Point(1, 1)};
  std::vector<std::array<int, 4>> mesh2Els = {{0, 1, 2, 3}, {1, 5, 2, -1}, {4, 0, 3, -1}};

  setup_mpmd(mesh1Verts, mesh1Els, mesh2Verts, mesh2Els);

  run_tests();

  if (m_meshInOnMesh1Procs)
    stk::middle_mesh::test_util::test_number_of_elements(m_meshInOnMesh1Procs, 10);
}


TEST_F(MiddleGridTriangulatorParallelTester, 2x2IdenticalSPMD)
{
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (!(commSize == 1 || commSize == 2 || commSize == 4))
    GTEST_SKIP();


  mesh::impl::MeshSpec spec1, spec2;
  spec1.xmin = 0; spec2.xmin = 0;
  spec1.xmax = 1; spec2.xmax = 1;
  spec1.ymin = 0; spec2.ymin = 0;
  spec1.ymax = 1; spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 2;
  spec1.numelY = 2; spec2.numelY = 2;

  setup_spmd(spec1, spec2, false);
  run_tests();
}

TEST_F(MiddleGridTriangulatorParallelTester, 2x2IdenticalSPMDTri)
{
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (!(commSize == 1 || commSize == 2 || commSize == 4))
    GTEST_SKIP();


  mesh::impl::MeshSpec spec1, spec2;
  spec1.xmin = 0; spec2.xmin = 0;
  spec1.xmax = 1; spec2.xmax = 1;
  spec1.ymin = 0; spec2.ymin = 0;
  spec1.ymax = 1; spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 2;
  spec1.numelY = 2; spec2.numelY = 2;

  setup_spmd(spec1, spec2, true);
  run_tests();
}

TEST_F(MiddleGridTriangulatorParallelTester, 2x3MPMD)
{
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (!(commSize == 4 || commSize == 8))
    GTEST_SKIP();


  mesh::impl::MeshSpec spec1, spec2;
  spec1.xmin = 0; spec2.xmin = 0;
  spec1.xmax = 1; spec2.xmax = 1;
  spec1.ymin = 0; spec2.ymin = 0;
  spec1.ymax = 1; spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 3;
  spec1.numelY = 2; spec2.numelY = 3;

  SplitCommTestUtil splitter(commSize/2, commSize/2);

  setup_mpmd(spec1, spec2, splitter);
  run_tests();
}

TEST_F(MiddleGridTriangulatorParallelTester, 2x3MPMDTri)
{
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (!(commSize == 4 || commSize == 8))
    GTEST_SKIP();


  mesh::impl::MeshSpec spec1, spec2;
  spec1.xmin = 0; spec2.xmin = 0;
  spec1.xmax = 1; spec2.xmax = 1;
  spec1.ymin = 0; spec2.ymin = 0;
  spec1.ymax = 1; spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 3;
  spec1.numelY = 2; spec2.numelY = 3;

  SplitCommTestUtil splitter(commSize/2, commSize/2);

  setup_mpmd(spec1, spec2, splitter, true);
  run_tests();
}


TEST_F(MiddleGridTriangulatorParallelTester, 2x4MPMD)
{
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (!(commSize == 4 || commSize == 8))
    GTEST_SKIP();


  mesh::impl::MeshSpec spec1, spec2;
  spec1.xmin = 0; spec2.xmin = 0;
  spec1.xmax = 1; spec2.xmax = 1;
  spec1.ymin = 0; spec2.ymin = 0;
  spec1.ymax = 1; spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 4;
  spec1.numelY = 2; spec2.numelY = 4;

  SplitCommTestUtil splitter(commSize/2, commSize/2);

  setup_mpmd(spec1, spec2, splitter);
  run_tests();
}

TEST_F(MiddleGridTriangulatorParallelTester, 2x4MPMDTri)
{
  int commSize = utils::impl::comm_size(MPI_COMM_WORLD);
  if (!(commSize == 4 || commSize == 8))
    GTEST_SKIP();


  mesh::impl::MeshSpec spec1, spec2;
  spec1.xmin = 0; spec2.xmin = 0;
  spec1.xmax = 1; spec2.xmax = 1;
  spec1.ymin = 0; spec2.ymin = 0;
  spec1.ymax = 1; spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 4;
  spec1.numelY = 2; spec2.numelY = 4;

  SplitCommTestUtil splitter(commSize/2, commSize/2);

  setup_mpmd(spec1, spec2, splitter, true);
  run_tests();
}