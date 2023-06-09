#include <functional>
#include "gtest/gtest.h"
#include "stk_transfer/ConservativeTransfer.hpp"
#include "stk_unit_test_utils/ConservativeTransferUserExample.hpp"

namespace {

using namespace stk::middle_mesh;

class ConservativeTransferTester : public ::testing::Test
{
  public:
    ~ConservativeTransferTester() { free_comms(); }

    //void setup(int nprocs1, int nprocs2, int nprocs2Offset)
    void setup(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range)
    {
      free_comms();

      inputMesh1 = nullptr;
      inputMesh2 = nullptr;

      check_proc_range(proc1Range);
      check_proc_range(proc2Range);
      create_mesh_comms(proc1Range, proc2Range);
      create_input_meshes();

      if (inputMesh1)
        transferCallback1 = std::make_shared<ConservativeTransferUserForTest>(inputMesh1);

      if (inputMesh2)
        transferCallback2 = std::make_shared<ConservativeTransferUserForTest>(inputMesh2);

      conservativeTransfer = std::make_shared<stk::transfer::ConservativeTransfer>(MPI_COMM_WORLD, inputMesh1, inputMesh2,
                                                                                   transferCallback1, transferCallback2);
    }

    void test_exactness(std::function<double(const utils::Point&)> func)
    {
      mesh::FieldPtr<double> sendFieldPtr, recvFieldPtr;
      if (inputMesh1)
      {
        sendFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
        set_field(sendFieldPtr, func);
      }

      if (inputMesh2)
        recvFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);

      conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr);
      conservativeTransfer->finish_transfer();

      if (inputMesh2)
      {
        auto& recvField = *recvFieldPtr;
        for (auto vert : inputMesh2->get_vertices())
        {
          double valExpected = func(vert->get_point_orig(0));
          EXPECT_NEAR(recvField(vert, 0, 0), valExpected, 1e-12);
        }
      }

    }

    void test_exactness_bidirectional(std::function<double(const utils::Point&)> func1, 
                                      std::function<double(const utils::Point&)> func2)
    {
      assert(inputMesh1 || inputMesh2);
      mesh::FieldPtr<double> sendFieldPtr, recvFieldPtr;
      
      if (inputMesh1)
      {
        sendFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
        recvFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
        set_field(sendFieldPtr, func1);
      }

      if (inputMesh2)
      {
        sendFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
        recvFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
        set_field(sendFieldPtr, func2);
      }

      conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr);
      conservativeTransfer->finish_transfer();

      auto func = inputMesh1 ? func2 : func1;
      auto mesh = inputMesh1 ? inputMesh1 : inputMesh2;
      auto& recvField = *recvFieldPtr;
      for (auto vert : mesh->get_vertices())
      {
        double valExpected = func(vert->get_point_orig(0));
        EXPECT_NEAR(recvField(vert, 0, 0), valExpected, 1e-12);
      }
    }    

    void test_conservation(std::function<double(const utils::Point&)> func)
    {
      mesh::FieldPtr<double> sendFieldPtr, recvFieldPtr;
      double sendIntegral, recvIntegral;
      if (inputMesh1)
      {
        sendFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
        set_field(sendFieldPtr, func);
        sendIntegral = transferCallback1->integrate_function(sendFieldPtr);
      }

      if (inputMesh2)
        recvFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);

      conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr);
      conservativeTransfer->finish_transfer();

      if (inputMesh2)
        recvIntegral = transferCallback2->integrate_function(recvFieldPtr);

      int root = 0;

      MPI_Request sendFieldSendReq, sendFieldRecvReq, recvFieldSendReq, recvFieldRecvReq;
      auto sendFieldTag = stk::get_mpi_tag_manager().get_tag(MPI_COMM_WORLD);
      auto recvFieldTag = stk::get_mpi_tag_manager().get_tag(MPI_COMM_WORLD);


      if (inputMesh1)
        MPI_Isend(&sendIntegral, 1, MPI_DOUBLE, root, sendFieldTag, MPI_COMM_WORLD, &sendFieldSendReq);

      if (inputMesh2)
        MPI_Isend(&recvIntegral, 1, MPI_DOUBLE, root, recvFieldTag, MPI_COMM_WORLD, &recvFieldSendReq);

      if (utils::impl::comm_rank(MPI_COMM_WORLD) == root)
      {
        double sendIntegralRoot, recvIntegralRoot;

        MPI_Irecv(&sendIntegralRoot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, sendFieldTag, MPI_COMM_WORLD, &sendFieldRecvReq);
        MPI_Irecv(&recvIntegralRoot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, recvFieldTag, MPI_COMM_WORLD, &recvFieldRecvReq);

        MPI_Wait(&sendFieldRecvReq, MPI_STATUS_IGNORE);
        MPI_Wait(&recvFieldRecvReq, MPI_STATUS_IGNORE);

        EXPECT_NEAR(sendIntegralRoot, recvIntegralRoot, 1e-12);
      }

      if (inputMesh1)
        MPI_Wait(&sendFieldSendReq, MPI_STATUS_IGNORE);

      if (inputMesh2)
        MPI_Wait(&recvFieldSendReq, MPI_STATUS_IGNORE);
    }

    void test_conservation_bidirectional(std::function<double(const utils::Point&)> func1,
                                         std::function<double(const utils::Point&)> func2)
    {
      assert(inputMesh1 || inputMesh2);

      mesh::FieldPtr<double> sendFieldPtr, recvFieldPtr;
      double sendIntegral, recvIntegral;
      if (inputMesh1)
      {
        sendFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
        recvFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
        set_field(sendFieldPtr, func1);
        sendIntegral = transferCallback1->integrate_function(sendFieldPtr);
      }

      if (inputMesh2)
      {
        sendFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
        recvFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
        set_field(sendFieldPtr, func2);
        sendIntegral = transferCallback2->integrate_function(sendFieldPtr);
      }

      conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr);
      conservativeTransfer->finish_transfer();

      if (inputMesh1)
        recvIntegral = transferCallback1->integrate_function(recvFieldPtr);

      if (inputMesh2)
        recvIntegral = transferCallback2->integrate_function(recvFieldPtr);

      int root = 0;

      MPI_Request app1SendReq, app1RecvReq, app2SendReq, app2RecvReq;
      auto app1Tag = stk::get_mpi_tag_manager().get_tag(MPI_COMM_WORLD);
      auto app2Tag = stk::get_mpi_tag_manager().get_tag(MPI_COMM_WORLD);


      std::array<double, 2> integralVals = {sendIntegral, recvIntegral};
      if (inputMesh1)
        MPI_Isend(integralVals.data(), 2, MPI_DOUBLE, root, app1Tag, MPI_COMM_WORLD, &app1SendReq);

      if (inputMesh2)
        MPI_Isend(integralVals.data(), 2, MPI_DOUBLE, root, app2Tag, MPI_COMM_WORLD, &app2SendReq);

      if (utils::impl::comm_rank(MPI_COMM_WORLD) == root)
      {
        std::array<double, 2> app1IntegralVals, app2IntegralVals;

        MPI_Irecv(app1IntegralVals.data(), 2, MPI_DOUBLE, MPI_ANY_SOURCE, app1Tag, MPI_COMM_WORLD, &app1RecvReq);
        MPI_Irecv(app2IntegralVals.data(), 2, MPI_DOUBLE, MPI_ANY_SOURCE, app2Tag, MPI_COMM_WORLD, &app2RecvReq);

        MPI_Wait(&app1RecvReq, MPI_STATUS_IGNORE);
        MPI_Wait(&app2RecvReq, MPI_STATUS_IGNORE);

        EXPECT_NEAR(app1IntegralVals[0], app2IntegralVals[1], 1e-12);
        EXPECT_NEAR(app1IntegralVals[1], app2IntegralVals[0], 1e-12);
      }

      if (inputMesh1)
        MPI_Wait(&app1SendReq, MPI_STATUS_IGNORE);

      if (inputMesh2)
        MPI_Wait(&app2SendReq, MPI_STATUS_IGNORE);
    }    

    std::shared_ptr<mesh::Mesh> inputMesh1;
    std::shared_ptr<mesh::Mesh> inputMesh2;
    std::shared_ptr<ConservativeTransferUserForTest> transferCallback1;
    std::shared_ptr<ConservativeTransferUserForTest> transferCallback2;
    std::shared_ptr<stk::transfer::ConservativeTransfer> conservativeTransfer;


  private:

    void check_proc_range(std::pair<int, int> procRange)
    {
      std::cout << "proc range = " << procRange.first << ", " << procRange.second << std::endl;
      if ((procRange.second - procRange.first + 1) <= 0)
        throw std::runtime_error("nprocs cannot be zero");

      if (!(procRange.first >= 0 && procRange.second < utils::impl::comm_size(MPI_COMM_WORLD)))
        throw std::runtime_error("process range is invalid");
    }

    bool create_mesh_comms(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range)
    {
      int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
      int color1  = myrank >= proc1Range.first  && myrank <= proc1Range.second? 0 : MPI_UNDEFINED;
      MPI_Comm_split(MPI_COMM_WORLD, color1, 0, &m_meshComm1);

      int color2  = myrank >= proc2Range.first  && myrank <= proc2Range.second? 0 : MPI_UNDEFINED;
      MPI_Comm_split(MPI_COMM_WORLD, color2, 0, &m_meshComm2);

      return color1 == 0 || color2 == 0;
    }

    void free_comms()
    {
      if (m_meshComm1 != MPI_COMM_NULL)
        MPI_Comm_free(&m_meshComm1);

      if (m_meshComm2 != MPI_COMM_NULL)
        MPI_Comm_free(&m_meshComm2);
    }

    void create_input_meshes()
    {
      mesh::impl::MeshSpec spec1, spec2;
      spec1.xmin   = 0;
      spec2.xmin   = 0;
      spec1.xmax   = 1;
      spec2.xmax   = 1;
      spec1.ymin   = 0;
      spec2.ymin   = 0;
      spec1.ymax   = 1;
      spec2.ymax   = 1;
      spec1.numelX = 4;
      spec2.numelX = 5;
      spec1.numelY = 4;
      spec2.numelY = 5;

      auto f     = [](const utils::Point& pt) { return pt; };

      if (m_meshComm1 != MPI_COMM_NULL)
        inputMesh1 = create_mesh(spec1, f, m_meshComm1);

      if (m_meshComm2 != MPI_COMM_NULL)
        inputMesh2 = create_mesh(spec2, f, m_meshComm2);
    }

    void set_field(mesh::FieldPtr<double> fieldPtr, std::function<double(const utils::Point&)> func)
    {
      auto& field = *fieldPtr;
      for (auto vert : field.get_mesh()->get_vertices())
        if (vert)
        {
          field(vert, 0, 0) = func(vert->get_point_orig(0));
        }
    }

    MPI_Comm m_meshComm1 = MPI_COMM_NULL;
    MPI_Comm m_meshComm2 = MPI_COMM_NULL;
};

}


TEST_F(ConservativeTransferTester, SPMDLinear)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {0, 0});
  test_exactness(f);
  test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDExponential)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return std::exp(pt.x + pt.y); };

  setup({0, 0}, {0, 0});
  test_conservation(f);
}

TEST_F(ConservativeTransferTester, MPMDLinear)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto f1 = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  auto f2 = [](const utils::Point& pt) { return pt.x + 2*pt.y + 1; };


  setup({0, 0}, {1, 1});
  test_exactness(f1);
  test_conservation(f1);

  std::cout << "testing bidirectional " << std::endl;
  test_exactness_bidirectional(f1, f2);
  test_conservation_bidirectional(f1, f2);
}

TEST_F(ConservativeTransferTester, MPMDExponential)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto f1 = [](const utils::Point& pt) { return std::exp(pt.x + pt.y); };
  auto f2 = [](const utils::Point& pt) { return pt.x + 2*pt.y + 1; };


  setup({0, 0}, {1, 1});
  test_conservation(f1);

  test_conservation_bidirectional(f1, f2);
}

TEST_F(ConservativeTransferTester, SPMDExtraneousProcess)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {0, 0});
  test_exactness(f);
  test_conservation(f);
}

TEST_F(ConservativeTransferTester, MPMDExtraneousProcess)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {1, 1});
  test_exactness(f);
  test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDRootExtraneousProcess)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({1, 1}, {1, 1});
  test_exactness(f);
  test_conservation(f);
}

TEST_F(ConservativeTransferTester, StartBeforeFinishError)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  setup({0, 0}, {0, 0});
  mesh::FieldPtr<double> sendFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
  mesh::FieldPtr<double> recvFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
  conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr);
  EXPECT_ANY_THROW(conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr));
  conservativeTransfer->finish_transfer();
  
  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  test_exactness(f);
}

TEST_F(ConservativeTransferTester, FinishBeforeStartError)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  setup({0, 0}, {0, 0});
  mesh::FieldPtr<double> sendFieldPtr = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
  mesh::FieldPtr<double> recvFieldPtr = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
  EXPECT_ANY_THROW(conservativeTransfer->finish_transfer());

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  test_exactness(f);

  EXPECT_ANY_THROW(conservativeTransfer->finish_transfer());
}