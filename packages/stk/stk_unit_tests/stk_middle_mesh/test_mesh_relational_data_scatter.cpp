#include "gtest/gtest.h"
#include "stk_middle_mesh/predicates/intersection_common.hpp"
#include "stk_middle_mesh/mesh_relational_data_scatter.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_relational_data.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_scatter_spec.hpp"
#include "stk_middle_mesh/mesh_scatter.hpp"
#include "stk_middle_mesh/predicates/point_classifier_normal_wrapper.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"

namespace {

using namespace stk::middle_mesh;

class MeshRelationalDataScatterTesterBase : public ::testing::Test
{
  protected:
    MeshRelationalDataScatterTesterBase() :
      m_comm1(MPI_COMM_NULL),
      m_comm2(MPI_COMM_NULL),
      m_unionComm(MPI_COMM_WORLD)
    {}

    virtual ~MeshRelationalDataScatterTesterBase()
    {
      if (m_comm1 != MPI_COMM_NULL)
        MPI_Comm_free(&m_comm1);

      if (m_comm2 != MPI_COMM_NULL)
        MPI_Comm_free(&m_comm2);
    }

    void SetUp() override
    {
      int commsize = utils::impl::comm_size(m_unionComm);
      if (!(commsize == 1 || commsize == 2 || commsize == 4))
      {
        GTEST_SKIP();
      }

      split_comm();
      create_meshes();
      auto scatterSpec1 = create_scatter_spec1();
      auto scatterSpec2 = create_scatter_spec2();

      mesh::impl::MeshScatter scatter1(scatterSpec1, m_mesh1, m_comm2, true);
      m_mesh1ScatteredToMesh2      = scatter1.scatter();
      auto mesh1EntityDestinations = scatter1.get_entity_destinations();
      auto mesh1EntityOrigins      = scatter1.get_entity_origins();

      mesh::impl::MeshScatter scatter2(scatterSpec2, m_mesh2, m_comm1, true);
      m_mesh2ScatteredToMesh1      = scatter2.scatter();
      auto mesh2EntityDestinations = scatter2.get_entity_destinations();
      auto mesh2EntityOrigins      = scatter2.get_entity_origins();

      create_mesh_relational_data_input();

      nonconformal4::impl::MeshRelationalDataScatterInput scatterInput;
      scatterInput.mesh1                   = m_mesh1;
      scatterInput.mesh2                   = m_mesh2;
      scatterInput.mesh1ScatteredToMesh2   = m_mesh1ScatteredToMesh2;
      scatterInput.mesh2ScatteredToMesh1   = m_mesh2ScatteredToMesh1;
      scatterInput.mesh1EntityOrigins      = mesh1EntityOrigins;
      scatterInput.mesh1EntityDestinations = mesh1EntityDestinations;
      scatterInput.mesh2EntityOrigins      = mesh2EntityOrigins;
      scatterInput.mesh2EntityDestinations = mesh2EntityDestinations;

      m_meshRelationalDataScatter = std::make_shared<nonconformal4::impl::MeshRelationalDataScatter>(
                                      scatterInput,
                                      m_meshRelationalDataInput,
                                      m_pointClassifierOrigin, m_pointClassifierDest,
                                      m_unionComm);
    }


    void split_comm()
    {
      int myrank = utils::impl::comm_rank(m_unionComm);
      int commsize = utils::impl::comm_size(m_unionComm);

      if (commsize == 1)
      {
        MPI_Comm_dup(m_unionComm, &m_comm1);
        MPI_Comm_dup(m_unionComm, &m_comm2);
      } else if (commsize == 2 || commsize == 4)
      {
        int color = myrank < commsize/2 ? 0 : 1;
        MPI_Comm* newcomm = myrank < commsize/2 ? &m_comm1 : &m_comm2;
        MPI_Comm_split(m_unionComm, color, 0, newcomm);
      } else
      {
        throw std::runtime_error("MeshRelationalDataScatter tests only run with 1, 2, or 4 procs");
      }
    }

    virtual void create_meshes() = 0;

    virtual std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec1() = 0;

    virtual std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec2() = 0;
  
    void create_mesh_relational_data_input()
    {
      if (m_mesh2)
      {

        m_meshRelationalDataInput = std::make_shared<nonconformal4::impl::MeshRelationalData>(m_mesh1ScatteredToMesh2, m_mesh2, nullptr);
        m_pointClassifierOrigin = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(m_mesh2);

        // put enough data into the fields so the MeshRelationalScatter will complete
        // without segfaulting due to nullptrs
        // Each test can overwrite this data as needed
        auto& verts1ToFakeVerts = *(m_meshRelationalDataInput->verts1ToFakeVerts);
        for (auto& vert1 : m_mesh1ScatteredToMesh2->get_vertices())
          if (vert1)
          {
            verts1ToFakeVerts(vert1, 0, 0) = m_vertGenerator.get_vert(vert1->get_point_orig(0));
          }

        auto& verts2ToFakeVerts = *(m_meshRelationalDataInput->verts2ToFakeVerts);
        auto& verts2ClassOnMesh1 = *(m_meshRelationalDataInput->verts2ClassOnMesh1);
        mesh::MeshEntityPtr el1 = m_mesh1ScatteredToMesh2->get_elements()[0];
        for (auto& vert2 : m_mesh2->get_vertices())
          if (vert2)
          {
            verts2ToFakeVerts(vert2, 0, 0) = m_vertGenerator.get_vert(vert2->get_point_orig(0));
            verts2ClassOnMesh1(vert2, 0, 0) = predicates::impl::PointRecord(predicates::impl::PointClassification::Vert, 0, el1);
          }
      }
      
      if (m_mesh1)
      {
        m_pointClassifierDest = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(m_mesh2ScatteredToMesh1);
      }
    }

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh1ScatteredToMesh2;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_mesh2ScatteredToMesh1;
    nonconformal4::impl::FakeVertGenerator m_vertGenerator;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> m_meshRelationalDataInput;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierOrigin;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierDest;
    std::shared_ptr<nonconformal4::impl::MeshRelationalDataScatter> m_meshRelationalDataScatter;
    MPI_Comm m_comm1;
    MPI_Comm m_comm2;
    MPI_Comm m_unionComm;
};

class MeshRelationalDataScatterTester : public MeshRelationalDataScatterTesterBase
{
  protected:

    void create_meshes() override
    {
      mesh::impl::MeshSpec spec1, spec2;
      spec1.xmin = 0;   spec1.xmax = 1;
      spec1.ymin = 0;   spec1.ymax = 1;
      spec1.numelX = 2; spec1.numelY = 2;
      
      spec2.xmin = 0;   spec2.xmax = 1;
      spec2.ymin = 0;   spec2.ymax = 1;
      spec2.numelX = 3; spec2.numelY = 3;   

      auto f = [](const utils::Point& pt) { return pt; };
      if (m_comm1 != MPI_COMM_NULL)
      {
        m_mesh1 = create_mesh(spec1, f, m_comm1);
      }   

      if (m_comm2 != MPI_COMM_NULL)
      {
        m_mesh2 = create_mesh(spec2, f, m_comm2);
      }            
    }

    std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec1() override
    {
      int myrank = utils::impl::comm_rank(m_unionComm);
      int commsize = utils::impl::comm_size(m_unionComm);
      auto scatterSpec1 = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_mesh1);
      if (m_mesh1)
      {
        if (commsize == 1)
        {
          for (auto el : m_mesh1->get_elements())
            if (el)
              scatterSpec1->add_destination(el, 0);
        } else if (commsize == 2)
        {
          for (auto el : m_mesh1->get_elements())
            if (el)
              scatterSpec1->add_destination(el, 1);
        } else
        {
          for (auto el : m_mesh1->get_elements())
            if (el)
            {
              scatterSpec1->add_destination(el, 2);
              if (myrank == 1)
                scatterSpec1->add_destination(el, 3);
            }
        }
      }

      return scatterSpec1;
    }

    std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec2() override
    {
      int myrank = utils::impl::comm_rank(m_unionComm);
      int commsize = utils::impl::comm_size(m_unionComm);
      auto scatterSpec2 = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_mesh2);
      if (m_mesh2)
      {
        if (commsize == 1)
        {
          for (auto el : m_mesh2->get_elements())
            if (el)
              scatterSpec2->add_destination(el, 0);
        } else if (commsize == 2)
        {
          for (auto el : m_mesh2->get_elements())
            if (el)
              scatterSpec2->add_destination(el, 0);          
        } else
        {
          for (auto el : m_mesh2->get_elements())
            if (el)
            {
              scatterSpec2->add_destination(el, myrank - 2);
              if (std::abs(mesh::compute_centroid(el).x - 0.5) < 1e-12)
                scatterSpec2->add_destination(el, 1);
            }          
        }
      }

      return scatterSpec2;
    }      
};

class MeshRelationalDataScatterTesterTri : public MeshRelationalDataScatterTesterBase
{
  protected:

    void create_meshes() override
    {
      mesh::impl::MeshSpec spec1, spec2;
      spec1.xmin = 0;   spec1.xmax = 1;
      spec1.ymin = 0;   spec1.ymax = 1;
      spec1.numelX = 2; spec1.numelY = 2;
      
      spec2.xmin = 0;   spec2.xmax = 1;
      spec2.ymin = 0;   spec2.ymax = 1;
      spec2.numelX = 3; spec2.numelY = 3;   

      auto f = [](const utils::Point& pt) { return pt; };
      if (m_comm1 != MPI_COMM_NULL)
      {
        m_mesh1 = create_mesh(spec1, f, m_comm1, true);
      }   

      if (m_comm2 != MPI_COMM_NULL)
      {
        m_mesh2 = create_mesh(spec2, f, m_comm2, true);
      }            
    }

    std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec1() override
    {
      int myrank = utils::impl::comm_rank(m_unionComm);
      int commsize = utils::impl::comm_size(m_unionComm);
      auto scatterSpec1 = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_mesh1);
      if (m_mesh1)
      {
        if (commsize == 1)
        {
          for (auto el : m_mesh1->get_elements())
            if (el)
              scatterSpec1->add_destination(el, 0);
        } else if (commsize == 2)
        {
          for (auto el : m_mesh1->get_elements())
            if (el)
              scatterSpec1->add_destination(el, 1);
        } else
        {
          for (auto el : m_mesh1->get_elements())
            if (el)
            {
              scatterSpec1->add_destination(el, 2);
              if (myrank == 1)
                scatterSpec1->add_destination(el, 3);
            }
        }
      }

      return scatterSpec1;
    }

    std::shared_ptr<mesh::impl::MeshScatterSpec> create_scatter_spec2() override
    {
      int myrank = utils::impl::comm_rank(m_unionComm);
      int commsize = utils::impl::comm_size(m_unionComm);
      auto scatterSpec2 = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_mesh2);
      if (m_mesh2)
      {
        if (commsize == 1)
        {
          for (auto el : m_mesh2->get_elements())
            if (el)
              scatterSpec2->add_destination(el, 0);
        } else if (commsize == 2)
        {
          for (auto el : m_mesh2->get_elements())
            if (el)
              scatterSpec2->add_destination(el, 0);          
        } else
        {
          for (auto el : m_mesh2->get_elements())
            if (el)
            {
              scatterSpec2->add_destination(el, myrank - 2);
              if (mesh::compute_centroid(el).x > 1.0/3 && mesh::compute_centroid(el).x < 2.0/3)
                scatterSpec2->add_destination(el, 1);
            }          
        }
      }

      return scatterSpec2;
    }      
};

mesh::MeshEntityPtr find_entity(std::shared_ptr<mesh::Mesh> mesh, int dim, const utils::Point& centroid, double tol)
{
  double minDist = std::numeric_limits<double>::max();
  mesh::MeshEntityPtr minEntity = nullptr;

  for (auto& entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      auto disp = mesh::compute_centroid(entity) - centroid;
      double dist = std::sqrt(dot(disp, disp));
      if (dist < tol && dist < minDist)
      {
        minDist = dist;
        minEntity = entity;
      }
    }

  return minEntity;
}

}  // namespace


TEST_F(MeshRelationalDataScatterTester, Mesh1Verts)
{
  if (m_mesh1ScatteredToMesh2)
  {
    auto& verts1ToFakeVerts = *(m_meshRelationalDataInput->verts1ToFakeVerts);

    for (auto& vert : m_mesh1ScatteredToMesh2->get_vertices())
      if (vert)
        verts1ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));
  }

  std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh1)
  {
    std::set<int> fakeVertIds;
    auto& verts1ToFakeVerts = *(meshRelationalDataOutput->verts1ToFakeVerts);
    int minFakeVertId = std::numeric_limits<int>::max();
    int maxFakeVertId = std::numeric_limits<int>::min();
    int nverts = 0;
    for (auto& vert : m_mesh1->get_vertices())
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);

        EXPECT_EQ(fakeVertIds.count(fv.id), 0u);
        fakeVertIds.insert(fv.id);

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);

        minFakeVertId = std::min(minFakeVertId, fv.id);
        maxFakeVertId = std::max(maxFakeVertId, fv.id);
        nverts++;
      }

    EXPECT_EQ(fakeVertIds.size(), size_t(mesh::count_valid(m_mesh1->get_vertices())));
    EXPECT_EQ(maxFakeVertId - minFakeVertId + 1, nverts);
  }
}


TEST_F(MeshRelationalDataScatterTester, Mesh2Verts)
{
  if (m_mesh2)
  {
    auto& verts2ToFakeVerts  = *(m_meshRelationalDataInput->verts2ToFakeVerts);
    auto& verts2ClassOnMesh1 = *(m_meshRelationalDataInput->verts2ClassOnMesh1);

    for (auto& vert : m_mesh2->get_vertices())
      if (vert)
      {
        verts2ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));

        for (auto& el : m_mesh1ScatteredToMesh2->get_elements())
          if (el)
          {
            mesh::MeshEntityPtr el2 = vert->get_up(0)->get_up(0);
            predicates::impl::PointRecord record = m_pointClassifierOrigin->classify(el, el2, vert->get_point_orig(0));
            if (record.type != predicates::impl::PointClassification::Exterior)
            {
              verts2ClassOnMesh1(vert, 0, 0) = record;   
            }           
          }
      }
  }

  std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    std::set<int> fakeVertIds;
    auto& verts2ToFakeVerts = *(meshRelationalDataOutput->verts2ToFakeVerts);
    auto& verts2ClassOnMesh1 = *(meshRelationalDataOutput->verts2ClassOnMesh1);

    int minFakeVertId = std::numeric_limits<int>::max();
    int maxFakeVertId = std::numeric_limits<int>::min();
    int nverts = 0;    
    for (auto& vert : m_mesh2ScatteredToMesh1->get_vertices())
    {
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);

        EXPECT_EQ(fakeVertIds.count(fv.id), 0u);
        fakeVertIds.insert(fv.id);

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);

        minFakeVertId = std::min(minFakeVertId, fv.id);
        maxFakeVertId = std::max(maxFakeVertId, fv.id);
        nverts++;

        predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert, 0, 0);
        if (record.el)
        {       
          mesh::MeshEntityPtr el2 = vert->get_up(0)->get_up(0);
          predicates::impl::PointRecord recordLocal = m_pointClassifierDest->classify(record.el, el2, vert->get_point_orig(0));

          EXPECT_EQ(record.type, recordLocal.type);
          EXPECT_EQ(record.el, recordLocal.el);
          EXPECT_EQ(record.id, recordLocal.id);
          EXPECT_EQ(m_pointClassifierDest->compute_xyz_coords(record), m_pointClassifierDest->compute_xyz_coords(recordLocal));
        }

      }
    }

    EXPECT_EQ(fakeVertIds.size(), size_t(mesh::count_valid(m_mesh2ScatteredToMesh1->get_vertices())));
    EXPECT_EQ(maxFakeVertId - minFakeVertId + 1, nverts);

    size_t numFakeVertsExpected = mesh::count_valid(m_mesh2ScatteredToMesh1->get_vertices()) +
                                  mesh::count_valid(m_mesh1->get_vertices());
    EXPECT_EQ(meshRelationalDataOutput->fakeVertsToVertsIn.size(), numFakeVertsExpected);
  }
}

TEST_F(MeshRelationalDataScatterTesterTri, Mesh2Verts)
{
  if (m_mesh2)
  {
    auto& verts2ToFakeVerts  = *(m_meshRelationalDataInput->verts2ToFakeVerts);
    auto& verts2ClassOnMesh1 = *(m_meshRelationalDataInput->verts2ClassOnMesh1);

    for (auto& vert : m_mesh2->get_vertices())
      if (vert)
      {
        verts2ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));

        for (auto& el : m_mesh1ScatteredToMesh2->get_elements())
          if (el)
          {
            mesh::MeshEntityPtr el2 = vert->get_up(0)->get_up(0);
            predicates::impl::PointRecord record = m_pointClassifierOrigin->classify(el, el2, vert->get_point_orig(0));
            if (record.type != predicates::impl::PointClassification::Exterior)
            {
              verts2ClassOnMesh1(vert, 0, 0) = record;   
            }           
          }
      }
  }

  std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    std::set<int> fakeVertIds;
    auto& verts2ToFakeVerts = *(meshRelationalDataOutput->verts2ToFakeVerts);
    auto& verts2ClassOnMesh1 = *(meshRelationalDataOutput->verts2ClassOnMesh1);

    int minFakeVertId = std::numeric_limits<int>::max();
    int maxFakeVertId = std::numeric_limits<int>::min();
    int nverts = 0;    
    for (auto& vert : m_mesh2ScatteredToMesh1->get_vertices())
    {
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);

        EXPECT_EQ(fakeVertIds.count(fv.id), 0u);
        fakeVertIds.insert(fv.id);

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);

        minFakeVertId = std::min(minFakeVertId, fv.id);
        maxFakeVertId = std::max(maxFakeVertId, fv.id);
        nverts++;

        predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert, 0, 0);
        if (record.el)
        {       
          mesh::MeshEntityPtr el2 = vert->get_up(0)->get_up(0);
          predicates::impl::PointRecord recordLocal = m_pointClassifierDest->classify(record.el, el2, vert->get_point_orig(0));

          EXPECT_EQ(record.type, recordLocal.type);
          EXPECT_EQ(record.el, recordLocal.el);
          EXPECT_EQ(record.id, recordLocal.id);
          EXPECT_EQ(m_pointClassifierDest->compute_xyz_coords(record), m_pointClassifierDest->compute_xyz_coords(recordLocal));
        }

      }
    }

    EXPECT_EQ(fakeVertIds.size(), size_t(mesh::count_valid(m_mesh2ScatteredToMesh1->get_vertices())));
    EXPECT_EQ(maxFakeVertId - minFakeVertId + 1, nverts);

    size_t numFakeVertsExpected = mesh::count_valid(m_mesh2ScatteredToMesh1->get_vertices()) +
                                  mesh::count_valid(m_mesh1->get_vertices());
    EXPECT_EQ(meshRelationalDataOutput->fakeVertsToVertsIn.size(), numFakeVertsExpected);
  }
}


TEST_F(MeshRelationalDataScatterTester, Mesh2ClassificationNearBoundary)
{
  if (m_mesh2)
  {
    auto& verts2ClassOnMesh1 = *(m_meshRelationalDataInput->verts2ClassOnMesh1);

    auto vert2 = find_entity(m_mesh2, 0, {2.0/3.0, 2.0/3.0, 0}, 1e-13);
    if (vert2)
    {
      auto el1   = find_entity(m_mesh1ScatteredToMesh2, 2, {0.75, 0.25}, 1e-13);
      verts2ClassOnMesh1(vert2, 0, 0) = predicates::impl::PointRecord(predicates::impl::PointClassification::Vert, 0, el1);
    }
  }

  std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    auto& verts2ClassOnMesh1 = *(meshRelationalDataOutput->verts2ClassOnMesh1);
    auto vert2 = find_entity(m_mesh2ScatteredToMesh1, 0, {2.0/3.0, 2.0/3.0, 0}, 1e-13);
    predicates::impl::PointRecord record = verts2ClassOnMesh1(vert2, 0, 0);
    EXPECT_NE(record.el, nullptr);
  }
}

TEST_F(MeshRelationalDataScatterTester, Mesh1VertOnMesh2Vert)
{
  if (m_mesh2)
  {
    auto& verts1ToFakeVerts = *(m_meshRelationalDataInput->verts1ToFakeVerts);
    auto& verts2ToFakeVerts = *(m_meshRelationalDataInput->verts2ToFakeVerts);

    for (auto& vert : m_mesh1ScatteredToMesh2->get_vertices())
      if (vert)
      {
        verts1ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));
      }

    for (auto& vert : m_mesh2->get_vertices())
      if (vert)
      {
        if (utils::impl::comm_rank(m_mesh2->get_comm()) == 0 && vert->get_id() == 0)
        {
          verts2ToFakeVerts(vert, 0, 0) = verts1ToFakeVerts(m_mesh1ScatteredToMesh2->get_vertices()[0], 0, 0);
        } else
        {
          verts2ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));
        }
      }
  }

  std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    std::set<int> fakeVertIds;
    auto& verts1ToFakeVerts = *(meshRelationalDataOutput->verts1ToFakeVerts);
    auto& verts2ToFakeVerts = *(meshRelationalDataOutput->verts2ToFakeVerts);
    for (auto& vert : m_mesh1->get_vertices())
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);
        fakeVertIds.insert(fv.id);  

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);          
      }

    for (auto& vert : m_mesh2ScatteredToMesh1->get_vertices())
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);

        fakeVertIds.insert(fv.id);

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);
      }

    int numUniqueIds = fakeVertIds.size();
    int numUniqueIdsExpected = mesh::count_valid(m_mesh1->get_vertices()) +
                               mesh::count_valid(m_mesh2ScatteredToMesh1->get_vertices());

    if (utils::impl::comm_rank(m_mesh1->get_comm()) == 0)
      --numUniqueIdsExpected;

    numUniqueIds = stk::get_global_sum(m_mesh1->get_comm(), numUniqueIds);
    numUniqueIdsExpected = stk::get_global_sum(m_mesh1->get_comm(), numUniqueIdsExpected);
    EXPECT_EQ(numUniqueIds, numUniqueIdsExpected);
  }
}


TEST_F(MeshRelationalDataScatterTester, Mesh2VertOnMesh1Vert)
{
  if (m_mesh2)
  {
    auto& verts1ToFakeVerts = *(m_meshRelationalDataInput->verts1ToFakeVerts);
    auto& verts2ToFakeVerts = *(m_meshRelationalDataInput->verts2ToFakeVerts);

    for (auto& vert : m_mesh2->get_vertices())
      if (vert)
      {
        verts2ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));
      }

    for (auto& vert : m_mesh1ScatteredToMesh2->get_vertices())
      if (vert)
      {
        if (utils::impl::comm_rank(m_mesh1ScatteredToMesh2->get_comm()) == 0 && vert->get_id() == 0)
        {
          verts1ToFakeVerts(vert, 0, 0) = verts2ToFakeVerts(m_mesh2->get_vertices()[0], 0, 0);
        } else
        {
          verts1ToFakeVerts(vert, 0, 0) = m_vertGenerator.get_vert(vert->get_point_orig(0));
        }
      }
  }

  std::shared_ptr<nonconformal4::impl::MeshRelationalData> meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    std::set<int> fakeVertIds;
    auto& verts1ToFakeVerts = *(meshRelationalDataOutput->verts1ToFakeVerts);
    auto& verts2ToFakeVerts = *(meshRelationalDataOutput->verts2ToFakeVerts);
    for (auto& vert : m_mesh1->get_vertices())
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts1ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);
        fakeVertIds.insert(fv.id);  

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);          
      }
      
    for (auto& vert : m_mesh2ScatteredToMesh1->get_vertices())
      if (vert)
      {
        nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(vert, 0, 0);
        utils::Point pt = vert->get_point_orig(0);

        fakeVertIds.insert(fv.id);

        double dist = std::sqrt(dot(fv.pt - pt, fv.pt - pt));
        EXPECT_NEAR(dist, 0, 1e-13);
      }

    int numUniqueIds = fakeVertIds.size();
    int numUniqueIdsExpected = mesh::count_valid(m_mesh1->get_vertices()) +
                               mesh::count_valid(m_mesh2ScatteredToMesh1->get_vertices());

    if (utils::impl::comm_rank(m_mesh1->get_comm()) == 0)
      --numUniqueIdsExpected;

    numUniqueIds = stk::get_global_sum(m_mesh1->get_comm(), numUniqueIds);
    numUniqueIdsExpected = stk::get_global_sum(m_mesh1->get_comm(), numUniqueIdsExpected);
    EXPECT_EQ(numUniqueIds, numUniqueIdsExpected);
  }
}


TEST_F(MeshRelationalDataScatterTester, mesh1EdgeFields)
{
  utils::Point pt(0, 1.0/3, 0);
  double xi = 2.0/3.0;
  if (m_mesh1ScatteredToMesh2)
  {
    auto& mesh1EdgesToSplit = *(m_meshRelationalDataInput->mesh1EdgesToSplit);
    mesh::MeshEntityPtr edge1 = find_entity(m_mesh1ScatteredToMesh2, 1, {0.0, 0.25, 0}, 1e-13);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2, 1, {0, 1.0/6, 0}, 1e-12);

    if (edge1)
    {
      nonconformal4::impl::FakeVert fv = m_vertGenerator.get_vert(pt);
      nonconformal4::impl::EdgeSplitRecord edgeSplit{fv, xi, edge2};
      mesh1EdgesToSplit.insert(edge1, 0, edgeSplit);
    }
  }

  auto meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh1)
  {
    auto& mesh1EdgesToSplit = *(meshRelationalDataOutput->mesh1EdgesToSplit);
    mesh::MeshEntityPtr edge1 = find_entity(m_mesh1, 1, {0, 0.25, 0}, 1e-13);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2ScatteredToMesh1, 1, {0, 1.0/6, 0}, 1e-12);
    
    if (edge1)
    {
      EXPECT_EQ(mesh1EdgesToSplit.get_num_comp(edge1, 0), 1);
      nonconformal4::impl::EdgeSplitRecord edgeSplit = mesh1EdgesToSplit(edge1, 0, 0);

      double dist = std::sqrt(dot(edgeSplit.vert.pt - pt, edgeSplit.vert.pt - pt));
      EXPECT_NEAR(dist, 0, 1e-12);
      EXPECT_NEAR(edgeSplit.xi, xi, 1e-13);
      EXPECT_EQ(edgeSplit.otherMeshEntity, edge2);
    }

  }
}

TEST_F(MeshRelationalDataScatterTester, mesh1EdgeFieldsSenderEntityShared)
{
  utils::Point pt(0, 1.0/3, 0);
  double xi = 0.5;
  if (m_mesh1ScatteredToMesh2)
  {
    auto& mesh1EdgesToSplit = *(m_meshRelationalDataInput->mesh1EdgesToSplit);
    mesh::MeshEntityPtr edge1 = find_entity(m_mesh1ScatteredToMesh2, 1, {0.75, 0.5, 0}, 1e-13);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2, 1, {2.0/3.0, 0.5, 0}, 1e-13);

    if (edge1)
    {
      nonconformal4::impl::FakeVert fv = m_vertGenerator.get_vert(pt);
      nonconformal4::impl::EdgeSplitRecord edgeSplit{fv, xi, edge2};
      mesh1EdgesToSplit.insert(edge1, 0, edgeSplit);
    }
  }

  auto meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh1)
  {
    auto& mesh1EdgesToSplit = *(meshRelationalDataOutput->mesh1EdgesToSplit);
    mesh::MeshEntityPtr edge1 = find_entity(m_mesh1, 1, {0.75, 0.5, 0}, 1e-13);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2ScatteredToMesh1, 1, {2.0/3.0, 0.5, 0}, 1e-13);
    
    if (edge1)
    {
      EXPECT_EQ(mesh1EdgesToSplit.get_num_comp(edge1, 0), 1);
      nonconformal4::impl::EdgeSplitRecord edgeSplit = mesh1EdgesToSplit(edge1, 0, 0);

      double dist = std::sqrt(dot(edgeSplit.vert.pt - pt, edgeSplit.vert.pt - pt));
      EXPECT_NEAR(dist, 0, 1e-12);
      EXPECT_NEAR(edgeSplit.xi, xi, 1e-13);
      EXPECT_EQ(edgeSplit.otherMeshEntity, edge2);
    }

  }
}


TEST_F(MeshRelationalDataScatterTester, mesh1EdgeFieldsWithRepeatedVert)
{
  utils::Point pt(1.0/3, 0, 0);
  double xi = 2.0/3.0;
  if (m_mesh1ScatteredToMesh2)
  {
    auto& mesh1EdgesToSplit = *(m_meshRelationalDataInput->mesh1EdgesToSplit);
    mesh::MeshEntityPtr edge1 = find_entity(m_mesh1ScatteredToMesh2, 1, {0.25, 0, 0}, 1e-13);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2, 1, {1.0/6, 0, 0}, 1e-12);
    auto& verts2ToFakeVerts = *(m_meshRelationalDataInput->verts2ToFakeVerts);

    if (edge1)
    {
      nonconformal4::impl::FakeVert fv = verts2ToFakeVerts(edge2->get_down(1), 0, 0);
      nonconformal4::impl::EdgeSplitRecord edgeSplit{fv, xi, edge2};
      mesh1EdgesToSplit.insert(edge1, 0, edgeSplit);
    }
  }

  auto meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh1)
  {
    auto& mesh1EdgesToSplit = *(meshRelationalDataOutput->mesh1EdgesToSplit);
    mesh::MeshEntityPtr edge1 = find_entity(m_mesh1, 1, {0.25, 0, 0}, 1e-13);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2ScatteredToMesh1, 1, {1.0/6, 0, 0}, 1e-12);
    auto& verts2ToFakeVerts   = *(meshRelationalDataOutput->verts2ToFakeVerts);
    
    if (edge1)
    {
      EXPECT_EQ(mesh1EdgesToSplit.get_num_comp(edge1, 0), 1);
      nonconformal4::impl::EdgeSplitRecord edgeSplit = mesh1EdgesToSplit(edge1, 0, 0);

      double dist = std::sqrt(dot(edgeSplit.vert.pt - pt, edgeSplit.vert.pt - pt));
      EXPECT_NEAR(dist, 0, 1e-12);
      EXPECT_NEAR(edgeSplit.xi, xi, 1e-13);
      EXPECT_EQ(edgeSplit.otherMeshEntity, edge2);
      EXPECT_EQ(edgeSplit.vert.id, verts2ToFakeVerts(edge2->get_down(1), 0, 0).id);
    }

  }
}


TEST_F(MeshRelationalDataScatterTester, mesh2EdgeFields)
{
  utils::Point pt1(1.0/3, 0, 0), pt2(1.0/3, 1.0/3, 0);
  if (m_mesh2)
  {
    auto& edges2ToFakeVertsIn = *(m_meshRelationalDataInput->edges2ToFakeVertsIn);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2, 1, {1.0/3, 1.0/6, 0}, 1e-12);

    if (edge2)
    {
      nonconformal4::impl::FakeVert fv1 = m_vertGenerator.get_vert(edge2->get_down(0)->get_point_orig(0));
      nonconformal4::impl::FakeVert fv2 = m_vertGenerator.get_vert(edge2->get_down(1)->get_point_orig(0));
      edges2ToFakeVertsIn.insert(edge2, 0, {fv1, 0});
      edges2ToFakeVertsIn.insert(edge2, 0, {fv2, 1});
    }
  }

  auto meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    auto& edges2ToFakeVertsIn = *(meshRelationalDataOutput->edges2ToFakeVertsIn);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2ScatteredToMesh1, 1, {1.0/3, 1.0/6, 0}, 1e-12);
    
    if (edge2)
    {
      EXPECT_EQ(edges2ToFakeVertsIn.get_num_comp(edge2, 0), 2);
      nonconformal4::impl::VertOnEdge vert1 = edges2ToFakeVertsIn(edge2, 0, 0);
      nonconformal4::impl::VertOnEdge vert2 = edges2ToFakeVertsIn(edge2, 0, 1);

      double dist1 = std::sqrt(dot(vert1.vert.pt - pt1, vert1.vert.pt - pt1));
      double dist2 = std::sqrt(dot(vert2.vert.pt - pt2, vert2.vert.pt - pt2));

      EXPECT_NEAR(dist1, 0, 1e-12);
      EXPECT_NEAR(dist2, 0, 1e-12);
      EXPECT_NEAR(vert1.xi, 0, 1e-13);
      EXPECT_NEAR(vert2.xi, 1, 1e-13);
    }
  }
}

TEST_F(MeshRelationalDataScatterTester, mesh2EdgeFieldsSenderEntityShared)
{
  utils::Point pt1(2.0/3, 1.0/3, 0), pt2(2.0/3, 2.0/3, 0);
  if (m_mesh2)
  {
    auto& edges2ToFakeVertsIn = *(m_meshRelationalDataInput->edges2ToFakeVertsIn);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2, 1, {2.0/3, 0.5, 0}, 1e-12);

    if (edge2)
    {
      nonconformal4::impl::FakeVert fv1 = m_vertGenerator.get_vert(edge2->get_down(0)->get_point_orig(0));
      nonconformal4::impl::FakeVert fv2 = m_vertGenerator.get_vert(edge2->get_down(1)->get_point_orig(0));
      edges2ToFakeVertsIn.insert(edge2, 0, {fv1, 0});
      edges2ToFakeVertsIn.insert(edge2, 0, {fv2, 1});
    }
  }

  auto meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    auto& edges2ToFakeVertsIn = *(meshRelationalDataOutput->edges2ToFakeVertsIn);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2ScatteredToMesh1, 1, {2.0/3, 0.5, 0}, 1e-12);
    
    if (edge2)
    {
      EXPECT_EQ(edges2ToFakeVertsIn.get_num_comp(edge2, 0), 2);
      nonconformal4::impl::VertOnEdge vert1 = edges2ToFakeVertsIn(edge2, 0, 0);
      nonconformal4::impl::VertOnEdge vert2 = edges2ToFakeVertsIn(edge2, 0, 1);

      double dist1 = std::sqrt(dot(vert1.vert.pt - pt1, vert1.vert.pt - pt1));
      double dist2 = std::sqrt(dot(vert2.vert.pt - pt2, vert2.vert.pt - pt2));

      EXPECT_NEAR(dist1, 0, 1e-12);
      EXPECT_NEAR(dist2, 0, 1e-12);
      EXPECT_NEAR(vert1.xi, 0, 1e-13);
      EXPECT_NEAR(vert2.xi, 1, 1e-13);
    }
  }
}

TEST_F(MeshRelationalDataScatterTester, mesh2EdgeFieldsWithRepeatedVert)
{
  utils::Point pt1(0, 0, 0), pt2(1.0/3, 0, 0);
  if (m_mesh2)
  {
    auto& verts2ToFakeVerts   = *(m_meshRelationalDataInput->verts2ToFakeVerts);
    auto& edges2ToFakeVertsIn = *(m_meshRelationalDataInput->edges2ToFakeVertsIn);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2, 1, {1.0/6, 0, 0}, 1e-12);

    if (edge2)
    {
      nonconformal4::impl::FakeVert fv1 = verts2ToFakeVerts(edge2->get_down(0), 0, 0);
      nonconformal4::impl::FakeVert fv2 = verts2ToFakeVerts(edge2->get_down(1), 0, 0);
      edges2ToFakeVertsIn.insert(edge2, 0, {fv1, 0});
      edges2ToFakeVertsIn.insert(edge2, 0, {fv2, 1});
    }
  }

  auto meshRelationalDataOutput = m_meshRelationalDataScatter->scatter();

  if (m_mesh2ScatteredToMesh1)
  {
    auto& verts2ToFakeVerts   = *(meshRelationalDataOutput->verts2ToFakeVerts);
    auto& edges2ToFakeVertsIn = *(meshRelationalDataOutput->edges2ToFakeVertsIn);
    mesh::MeshEntityPtr edge2 = find_entity(m_mesh2ScatteredToMesh1, 1, {1.0/6, 0, 0}, 1e-12);
    
    if (edge2)
    {
      EXPECT_EQ(edges2ToFakeVertsIn.get_num_comp(edge2, 0), 2);
      nonconformal4::impl::VertOnEdge vert1 = edges2ToFakeVertsIn(edge2, 0, 0);
      nonconformal4::impl::VertOnEdge vert2 = edges2ToFakeVertsIn(edge2, 0, 1);

      double dist1 = std::sqrt(dot(vert1.vert.pt - pt1, vert1.vert.pt - pt1));
      double dist2 = std::sqrt(dot(vert2.vert.pt - pt2, vert2.vert.pt - pt2));

      EXPECT_NEAR(dist1, 0, 1e-12);
      EXPECT_NEAR(dist2, 0, 1e-12);
      EXPECT_NEAR(vert1.xi, 0, 1e-13);
      EXPECT_NEAR(vert2.xi, 1, 1e-13);
      EXPECT_EQ(vert1.vert.id, verts2ToFakeVerts(edge2->get_down(0), 0, 0).id);
      EXPECT_EQ(vert2.vert.id, verts2ToFakeVerts(edge2->get_down(1), 0, 0).id);

    }
  }
}