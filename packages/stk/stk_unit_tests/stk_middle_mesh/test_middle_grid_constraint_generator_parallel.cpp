#include "stk_middle_mesh/mesh_relational_data.hpp"
#include "stk_middle_mesh/edge_tracer_opts.hpp"
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
#include "util/parallel_search_test_util.hpp"
#include "util/nonconformal_interface_helpers.hpp"

namespace {

using namespace stk::middle_mesh;

class MiddleGridConstraintGeneratorParallelTester : public ::testing::Test
{
  protected:

    using BoundingBoxSearch = search::ElementToElementBoundingBoxSearch;

    void setup_spmd(const mesh::impl::MeshSpec& spec1, const mesh::impl::MeshSpec& /*spec2*/)
    {
      auto f = [](const utils::Point& pt) { return pt; };

      auto mesh1 = mesh::impl::create_mesh(spec1, f);
      auto mesh2 = mesh::impl::create_mesh(spec1, f);

      setup(mesh1, mesh2);
    }

    void setup_mpmd(const mesh::impl::MeshSpec& spec1, const mesh::impl::MeshSpec& spec2,
                    const SplitCommTestUtil& splitter)
    {
      auto f = [](const utils::Point& pt) { return pt; };

      std::shared_ptr<mesh::Mesh> mesh1, mesh2;

      if (splitter.get_color() == SplitCommColor::SEND)
      {
        mesh1 = mesh::impl::create_mesh(spec1, f, splitter.get_comm());
      }

      if (splitter.get_color() == SplitCommColor::RECV)
      {
        mesh2 = mesh::impl::create_mesh(spec2, f, splitter.get_comm());
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
    }


    std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec(
        std::shared_ptr<mesh::Mesh> sendMesh, std::shared_ptr<mesh::Mesh> recvMesh)
    {
      //TODO: add some tests where union communication is not MPI_COMM_WORLD
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
};

mesh::MeshEntityPtr find_entity(std::shared_ptr<mesh::Mesh> mesh, int dim, const utils::Point& pt, double tol=1e-13)
{
  mesh::MeshEntityPtr minEntity = nullptr;
  double minDist = std::numeric_limits<double>::max();

  for (auto entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      utils::Point centroid = mesh::compute_centroid(entity);
      double dist = std::sqrt(dot(centroid - pt, centroid - pt));
      if (dist < minDist && dist < tol)
      {
        minDist = dist;
        minEntity = entity;
      }
    }

  if (minEntity == nullptr)
    throw std::runtime_error("unable to find entity");

  return minEntity;
}
}


TEST_F(MiddleGridConstraintGeneratorParallelTester, 2x2IdenticalSPMD)
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

  setup_spmd(spec1, spec2);

  auto meshIn = m_meshInOnMesh1Procs;
  if (meshIn)
  {
    EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_vertices()),
              mesh::count_valid(m_mesh1->get_vertices()));

    EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_edges()),
              mesh::count_valid(m_mesh1->get_edges()));

    if (utils::impl::comm_size(MPI_COMM_WORLD) == 2)
    {
      if (utils::impl::comm_rank(meshIn->get_comm()) == 0)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,   0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5, 0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0,   0.5, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5, 0.5, 0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0,   1.0, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5, 1.0, 0})->count_remote_shared_entities(), 1);

        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 0,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 0.5, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 1.0, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 0.75, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 0.75, 0}));        
      } else
      {
        EXPECT_EQ(find_entity(meshIn, 0, {1.0, 0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5, 0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0, 0.5, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5, 0.5, 0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0, 1.0, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5, 1.0, 0})->count_remote_shared_entities(), 1);         

        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 0,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 0.5, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 1.0, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 0.75, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 0.75, 0}));          
      }
    }

    mesh::check_topology(meshIn, 1);
  }
}

TEST_F(MiddleGridConstraintGeneratorParallelTester, 2x3MPMD)
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

  auto meshIn = m_meshInOnMesh1Procs;
  if (meshIn)
  {
    if (utils::impl::comm_size(MPI_COMM_WORLD) == 4)
    {
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_vertices()), 15);
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_edges()), 22);

      if (utils::impl::comm_rank(meshIn->get_comm()) == 0)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,     0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     1.0/3.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 1.0/3.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0/3.0,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.5,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 0.5,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     2.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 2.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   2.0/3,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0,   0})->count_remote_shared_entities(), 1);

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 1.0,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 1.0,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 1.0/6,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 5.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 1.0/6,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 5.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/6,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 10.0/12, 0}));

      } else
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0/3,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 1.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   1.0/3,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 0.5,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0.5,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   2.0/3,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 2.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   2.0/3,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   1.0,   0})->count_remote_shared_entities(), 0);   

        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 1.0,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 1.0,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/6,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 1.0/6,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 5.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 1.0/6,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 5.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 10.0/12, 0}));            
      }
    } else if (utils::impl::comm_size(MPI_COMM_WORLD) == 8)
    {
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_vertices()), 9);
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_edges()), 12);

      int myrank = utils::impl::comm_rank(meshIn->get_comm());

      if (myrank == 0)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,     0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     1.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 1.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0/3,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 3);    
    
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 0.5,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 0.5,   0})); 

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 1.0/6,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 5.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 1.0/6,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 5.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/6,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/12, 0}));        

      }  else if (myrank == 1)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0/3,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 1.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   1.0/3,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 3);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0.5,   0})->count_remote_shared_entities(), 1); 

        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 0.5,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 0,     0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 1.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 0.5,   0})); 

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/6,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 1.0/6,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 5.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 1.0/6,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 5.0/12, 0}));                
      } else if (myrank == 2)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 3);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     2.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 2.0/3,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   2.0/3,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0/3, 1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0,   0})->count_remote_shared_entities(), 1); 

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/6, 1.0,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/12, 1.0,   0})); 

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/3, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 10.0/12, 0}));               
      } else
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 3);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0.5,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   2.0/3, 0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 2.0/3, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   2.0/3, 0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {2.0/3, 1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   1.0,   0})->count_remote_shared_entities(), 0);

        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/12, 1.0,   0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 0.5,   0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 2.0/3, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {10.0/12, 1.0,   0})); 

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {2.0/3, 10.0/12, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 7.0/12,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 10.0/12, 0}));                   
      }
    }

    mesh::check_topology(meshIn, 1);
  }
}


TEST_F(MiddleGridConstraintGeneratorParallelTester, 2x4MPMD)
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

  auto meshIn = m_meshInOnMesh1Procs;
  if (meshIn)
  {
    if (utils::impl::comm_size(MPI_COMM_WORLD) == 4)
    {
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_vertices()), 15);
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_edges()), 22);

      if (utils::impl::comm_rank(meshIn->get_comm()) == 0)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,     0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.25,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0.25,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.25,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.5,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0.5,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,    0.75,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25, 0.75,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,  0.75,   0})->count_remote_shared_entities(), 1);

        EXPECT_EQ(find_entity(meshIn, 0, {0,    1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25, 1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,  1.0,   0})->count_remote_shared_entities(), 1);

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 3.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 7.0/8, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 3.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 7.0/8, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 3.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/8, 0}));

      } else
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.25, 0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.25, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0.25, 0})->count_remote_shared_entities(), 0);    

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,  0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.5,  0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0.5,  0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.75,  0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.75,  0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   0.75,  0})->count_remote_shared_entities(), 0);   

        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.0,   1.0,   0})->count_remote_shared_entities(), 0);

        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 3.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/8, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 3.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 7.0/8, 0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 3.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 7.0/8, 0}));                                 

      }
    } else if (utils::impl::comm_size(MPI_COMM_WORLD) == 8)
    {
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_vertices()), 9);
      EXPECT_EQ(mesh::count_valid(m_meshInOnMesh1Procs->get_edges()), 12);

      int myrank = utils::impl::comm_rank(meshIn->get_comm());

      if (myrank == 0)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,     0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0,   0})->count_remote_shared_entities(), 1);  

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.25, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0.25, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.25, 0})->count_remote_shared_entities(), 1);  

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.5,  0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0.5,  0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,  0})->count_remote_shared_entities(), 3); 

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.5,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.5,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 1.0/8,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 3.0/8,    0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 3.0/8, 0}));     

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/8,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 3.0/8,  0}));

      }  else if (myrank == 1)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0.50,  0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.00,  0,   0})->count_remote_shared_entities(), 0);      
        
        EXPECT_EQ(find_entity(meshIn, 0, {0.50,  0.25, 0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.25, 0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.00,  0.25, 0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.50,  0.5,   0})->count_remote_shared_entities(), 3);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.00,  0.5,   0})->count_remote_shared_entities(), 1);  

        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.5,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.25, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.5,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 1.0/8,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 3.0/8,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 1.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 3.0/8, 0}));     

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 1.0/8,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 3.0/8,  0}));                      
      } else if (myrank == 2)
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.5,   0})->count_remote_shared_entities(), 3);  

        EXPECT_EQ(find_entity(meshIn, 0, {0,     0.75,  0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  0.75,  0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   0.75,  0})->count_remote_shared_entities(), 1);  

        EXPECT_EQ(find_entity(meshIn, 0, {0,     1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.25,  1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {0.5,   1.0,   0})->count_remote_shared_entities(), 1);  

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {3.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0, 5.0/8,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0, 7.0/8,    0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 5.0/8, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.25, 7.0/8, 0}));     

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/8,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/8,  0}));              
      } else
      {
        EXPECT_EQ(find_entity(meshIn, 0, {0.50,  0.5,   0})->count_remote_shared_entities(), 3);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.5,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {1.00,  0.5,   0})->count_remote_shared_entities(), 1);      
        
        EXPECT_EQ(find_entity(meshIn, 0, {0.50,  0.75,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  0.75,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.00,  0.75,   0})->count_remote_shared_entities(), 0);

        EXPECT_EQ(find_entity(meshIn, 0, {0.50,  1.0,   0})->count_remote_shared_entities(), 1);
        EXPECT_EQ(find_entity(meshIn, 0, {0.75,  1.0,   0})->count_remote_shared_entities(), 0);
        EXPECT_EQ(find_entity(meshIn, 0, {1.00,  1.0,   0})->count_remote_shared_entities(), 0);

        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {5.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.5,  0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 0.75, 0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {7.0/8, 1.0,  0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 5.0/8,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.5, 7.0/8,    0}));

        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 5.0/8,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {0.75, 7.0/8,    0}));     

        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 5.0/8,    0}));
        EXPECT_TRUE(find_entity(meshIn, 1, {1.0, 7.0/8,    0}));                  
      }
    }

    mesh::check_topology(meshIn, 1);
  }
}