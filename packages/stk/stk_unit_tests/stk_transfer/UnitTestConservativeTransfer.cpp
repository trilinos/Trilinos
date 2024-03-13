#include <functional>
#include "gtest/gtest.h"
#include "stk_transfer/ConservativeTransfer.hpp"
#include "stk_transfer/ConservativeTransferUser.hpp"
#include "stk_unit_test_utils/ConservativeTransferUserExample.hpp"

#include "stk_unit_test_utils/stk_transfer_fixtures/ConservativeTransferFixture.hpp"

namespace {

using namespace stk::middle_mesh;

class ConservativeTransferTester : public ::testing::Test
{
  public:
    void setup(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range, ApplicationInterfaceType type)
    {
      comms  = stk::transfer::make_comm_splitter(proc1Range, proc2Range);
      meshes = std::make_shared<stk::transfer::MeshSetup>(comms->get_comm1(), comms->get_comm2());

      std::shared_ptr<stk::transfer::ConservativeTransferUser> transferCallback1, transferCallback2;
      if (meshes->get_mesh1())
        transferCallback1 = std::make_shared<ConservativeTransferUserForTest>(meshes->get_mesh1());

      if (meshes->get_mesh2())
        transferCallback2 = std::make_shared<ConservativeTransferUserForTest>(meshes->get_mesh2());

      tests  = std::make_shared<stk::transfer::ConservativeTransferTests>(
                meshes->get_mesh1(), meshes->get_mesh2(), transferCallback1, transferCallback2, type);
    }

   std::shared_ptr<stk::transfer::CommSplitter> comms;
   std::shared_ptr<stk::transfer::MeshSetup> meshes;
   std::shared_ptr<stk::transfer::ConservativeTransferTests> tests;
};

}


TEST_F(ConservativeTransferTester, SPMDLinear)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDLinearParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::Parallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDLinearVector)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return utils::Point{pt.x + 2*pt.y, 3*pt.x + 1.2*pt.y, pt.y + 2.1*pt.z}; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness_vector(f);
  //tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDLinearVectorParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return utils::Point{pt.x + 2*pt.y, 3*pt.x + 1.2*pt.y, pt.y + 2.1*pt.z}; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::Parallel);
  tests->test_exactness_vector(f);
  //tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDLinearScalarThenVector)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto fscalar = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  auto fvector = [](const utils::Point& pt) { return utils::Point{pt.x + 2*pt.y, 3*pt.x + 1.2*pt.y, pt.y + 2.1*pt.z}; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness(fscalar);
  tests->test_exactness_vector(fvector);
}

TEST_F(ConservativeTransferTester, SPMDLinearScalarThenVectorParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto fscalar = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  auto fvector = [](const utils::Point& pt) { return utils::Point{pt.x + 2*pt.y, 3*pt.x + 1.2*pt.y, pt.y + 2.1*pt.z}; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::Parallel);
  tests->test_exactness(fscalar);
  tests->test_exactness_vector(fvector);
}

TEST_F(ConservativeTransferTester, SPMDExponential)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return std::exp(pt.x + pt.y); };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDExponentialParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return std::exp(pt.x + pt.y); };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::Parallel);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, MPMDLinear)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto f1 = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  auto f2 = [](const utils::Point& pt) { return pt.x + 2*pt.y + 1; };


  setup({0, 0}, {1, 1}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness(f1);
  tests->test_conservation(f1);

  tests->test_exactness_bidirectional(f1, f2);
  tests->test_conservation_bidirectional(f1, f2);
}

TEST_F(ConservativeTransferTester, MPMDLinearParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto f1 = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  auto f2 = [](const utils::Point& pt) { return pt.x + 2*pt.y + 1; };


  setup({0, 0}, {1, 1}, ApplicationInterfaceType::Parallel);
  tests->test_exactness(f1);
  tests->test_conservation(f1);

  tests->test_exactness_bidirectional(f1, f2);
  tests->test_conservation_bidirectional(f1, f2);
}

TEST_F(ConservativeTransferTester, MPMDExponential)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto f1 = [](const utils::Point& pt) { return std::exp(pt.x + pt.y); };
  auto f2 = [](const utils::Point& pt) { return pt.x + 2*pt.y + 1; };


  setup({0, 0}, {1, 1}, ApplicationInterfaceType::FakeParallel);
  tests->test_conservation(f1);
  tests->test_conservation_bidirectional(f1, f2);
}

TEST_F(ConservativeTransferTester, MPMDExponentialParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  auto f1 = [](const utils::Point& pt) { return std::exp(pt.x + pt.y); };
  auto f2 = [](const utils::Point& pt) { return pt.x + 2*pt.y + 1; };


  setup({0, 0}, {1, 1}, ApplicationInterfaceType::Parallel);
  tests->test_conservation(f1);
  tests->test_conservation_bidirectional(f1, f2);
}

TEST_F(ConservativeTransferTester, SPMDExtraneousProcess)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDExtraneousProcessParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::Parallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, MPMDExtraneousProcess)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {1, 1}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, MPMDExtraneousProcessParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({0, 0}, {1, 1}, ApplicationInterfaceType::Parallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDRootExtraneousProcess)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({1, 1}, {1, 1}, ApplicationInterfaceType::FakeParallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, SPMDRootExtraneousProcessParallel)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };

  setup({1, 1}, {1, 1}, ApplicationInterfaceType::Parallel);
  tests->test_exactness(f);
  tests->test_conservation(f);
}

TEST_F(ConservativeTransferTester, StartBeforeFinishError)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  mesh::FieldPtr<double> sendFieldPtr = mesh::create_field<double>(meshes->get_mesh1(), mesh::FieldShape(1, 0, 0), 1);
  mesh::FieldPtr<double> recvFieldPtr = mesh::create_field<double>(meshes->get_mesh2(), mesh::FieldShape(1, 0, 0), 1);
  tests->conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr);
  EXPECT_ANY_THROW(tests->conservativeTransfer->start_transfer(sendFieldPtr, recvFieldPtr));
  tests->conservativeTransfer->finish_transfer();
  
  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  tests->test_exactness(f);
}

TEST_F(ConservativeTransferTester, FinishBeforeStartError)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  setup({0, 0}, {0, 0}, ApplicationInterfaceType::FakeParallel);
  EXPECT_ANY_THROW(tests->conservativeTransfer->finish_transfer());

  auto f = [](const utils::Point& pt) { return pt.x + 2*pt.y; };
  tests->test_exactness(f);

  EXPECT_ANY_THROW(tests->conservativeTransfer->finish_transfer());
}