#include "gtest/gtest.h"

#include "stk_middle_mesh/field_scatter_from_root.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh/mesh_scatter_from_root.hpp"

namespace {

using namespace stk::middle_mesh;

class FieldScatterFromRootTester
{
  public:
    FieldScatterFromRootTester(MPI_Comm unionComm, int rootRankOnUnionComm, MPI_Comm meshComm, std::vector<int> meshCommRanksOnUnionComm) :
      m_unionComm(unionComm),
      m_rootRankOnUnionComm(rootRankOnUnionComm),
      m_meshComm(meshComm),
      m_meshCommRanksOnUnionComm(meshCommRanksOnUnionComm)
    {
      if (utils::impl::comm_rank(m_unionComm) == m_rootRankOnUnionComm)
      {
        create_serial_mesh();
      }      
    }

    void runtest()
    {
      scatter_mesh_and_field();
    }

    void create_serial_mesh()
    {
      mesh::impl::MeshSpec spec;
      spec.xmin = 0; spec.ymin = 0;
      spec.xmax = 1; spec.ymax = 1;
      spec.numelX = 4; spec.numelY = 4;

      auto f = [](const utils::Point& pt) { return pt; };

      m_meshSerial = mesh::impl::create_mesh(spec, f, MPI_COMM_SELF);
    }

    mesh::FieldPtr<utils::Point> create_serial_field()
    {
      auto field = mesh::create_field<utils::Point>(m_meshSerial, mesh::impl::FieldShape(1, 1, 1), 1);

      for (int dim=0; dim < 3; ++dim)
        for (auto& entity : m_meshSerial->get_mesh_entities(dim))
          if (entity)
            (*field)(entity, 0, 0) = mesh::compute_centroid(entity);

      return field;
    }

    mesh::FieldPtr<int> create_element_destination_field()
    {
      auto field = mesh::create_field<int>(m_meshSerial, mesh::impl::FieldShape(0, 0, 1), 1);

      size_t idx = 0;
      for (auto& el : m_meshSerial->get_elements())
        if (el)
        {
          (*field)(el, 0, 0) = m_meshCommRanksOnUnionComm[idx];
          idx = (idx + 1) % m_meshCommRanksOnUnionComm.size();
        }

      return field;
    }

    void scatter_mesh_and_field()
    {
      bool amIRoot     = utils::impl::comm_rank(m_unionComm) == m_rootRankOnUnionComm;
      bool amIReceiver = m_meshComm != MPI_COMM_NULL;

      mesh::FieldPtr<int> elementDestinationRanks;
      if (amIRoot)
        elementDestinationRanks = create_element_destination_field();

      mesh::impl::MeshScatterFromRoot scatterer(m_unionComm, m_meshSerial, m_meshComm,
                                                elementDestinationRanks); 
      std::shared_ptr<mesh::Mesh> meshParallel = scatterer.scatter();
      mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityDestinations = scatterer.get_entity_destinations();


      mesh::FieldPtr<utils::Point> fieldSrc, fieldDest;
      if (amIRoot)
        fieldSrc = create_serial_field();

      if (amIReceiver)
        fieldDest = mesh::create_field<utils::Point>(meshParallel, mesh::impl::FieldShape(1, 1, 1), 1);

      mesh::impl::FieldScatterFromRoot<utils::Point> fieldScatterer(m_unionComm, m_rootRankOnUnionComm,
                                                            entityDestinations,
                                                            fieldSrc, fieldDest);
      fieldScatterer.scatter();

      if (amIReceiver)
        check_field(fieldDest);
    }

    void check_field(mesh::FieldPtr<utils::Point> fieldDest)
    {
      for (int dim=0; dim < 3; ++dim)
        for (auto& entity : fieldDest->get_mesh()->get_mesh_entities(dim))
          if (entity)
          {
            utils::Point centroidLocal = mesh::compute_centroid(entity);
            utils::Point centroidRecv  = (*fieldDest)(entity, 0, 0);

            for (int dim2=0; dim2 < 3; ++dim2)
              EXPECT_NEAR(centroidLocal[dim2], centroidRecv[dim2], 1e-13);
          }
    }

  private:
    MPI_Comm m_unionComm;
    int m_rootRankOnUnionComm;
    MPI_Comm m_meshComm;
    std::vector<int> m_meshCommRanksOnUnionComm;
    std::shared_ptr<mesh::Mesh> m_meshSerial;
};

}

TEST(FieldScatterFromRoot, RootOnMeshComm)
{
  std::vector<int> meshCommRanks;
  for (int i=0; i < utils::impl::comm_size(MPI_COMM_WORLD); ++i)
    meshCommRanks.push_back(i);

  FieldScatterFromRootTester tester(MPI_COMM_WORLD, 0, MPI_COMM_WORLD, meshCommRanks);
  tester.runtest();
}

TEST(FieldScatterFromRoot, RootNotOnMeshComm)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  //TODO: FIX THIS
  std::vector<int> meshCommRanks;
  for (int i=1; i < utils::impl::comm_size(MPI_COMM_WORLD); ++i)
    meshCommRanks.push_back(i);

  MPI_Comm meshComm;
  int color = utils::impl::comm_rank(MPI_COMM_WORLD) == 0 ? MPI_UNDEFINED : 0;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

  
  FieldScatterFromRootTester tester(MPI_COMM_WORLD, 0, meshComm, meshCommRanks);
  tester.runtest();

  if (meshComm != MPI_COMM_NULL)
    MPI_Comm_free(&meshComm);
}

TEST(FieldScatterFromRoot, ExtraneousProcesses)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) < 2)
    GTEST_SKIP();

  int color = utils::impl::comm_rank(MPI_COMM_WORLD) < (utils::impl::comm_size(MPI_COMM_WORLD) - 1) ? 0 : MPI_UNDEFINED;
  MPI_Comm meshComm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

  std::vector<int> meshCommRanks;
  for (int i=0; i < utils::impl::comm_size(MPI_COMM_WORLD) - 1; ++i)
    meshCommRanks.push_back(i);

  
  FieldScatterFromRootTester tester(MPI_COMM_WORLD, 0, meshComm, meshCommRanks);
  tester.runtest();

  if (meshComm != MPI_COMM_NULL)
    MPI_Comm_free(&meshComm);
}