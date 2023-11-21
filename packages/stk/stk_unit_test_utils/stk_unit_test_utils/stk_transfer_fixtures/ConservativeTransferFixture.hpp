#ifndef STK_TRANSFER_CONSERVATIVE_TRANSFER_TEST_H
#define STK_TRANSFER_CONSERVATIVE_TRANSFER_TEST_H

#include "gtest/gtest.h"

#include "stk_transfer/ConservativeTransfer.hpp"
#include "stk_unit_test_utils/ConservativeTransferUserExample.hpp"
#include <stk_middle_mesh/mesh.hpp>

namespace stk {
namespace transfer {

class CommSplitter
{
  public:
    CommSplitter(const CommSplitter&) = delete;

    CommSplitter& operator=(const CommSplitter&) = delete;

    ~CommSplitter();

    MPI_Comm get_comm1() const { return m_meshComm1; }

    MPI_Comm get_comm2() const { return m_meshComm2; }

  private:

    void check_proc_range(std::pair<int, int> procRange);

    bool create_mesh_comms(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range);

    void free_comms();

    CommSplitter(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range)
    {
      check_proc_range(proc1Range);
      check_proc_range(proc2Range);
      create_mesh_comms(proc1Range, proc2Range);
    }

    friend std::shared_ptr<CommSplitter> make_comm_splitter(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range);

    MPI_Comm m_meshComm1 = MPI_COMM_NULL;
    MPI_Comm m_meshComm2 = MPI_COMM_NULL;    
};

inline std::shared_ptr<CommSplitter> make_comm_splitter(std::pair<int, int> proc1Range, std::pair<int, int> proc2Range)
{
  return std::shared_ptr<CommSplitter>(new CommSplitter(proc1Range, proc2Range));
}

class MeshSetupBase
{
  public:
    virtual ~MeshSetupBase() = default;

    virtual std::shared_ptr<stk::middle_mesh::mesh::Mesh> get_mesh1() const = 0;

    virtual std::shared_ptr<stk::middle_mesh::mesh::Mesh> get_mesh2() const = 0; 
};

class MeshSetup : public MeshSetupBase
{
  public:
    MeshSetup(MPI_Comm comm1, MPI_Comm comm2)
    {
      create_input_meshes(comm1, comm2);
    }

    std::shared_ptr<stk::middle_mesh::mesh::Mesh> get_mesh1() const override { return m_inputMesh1; }

    std::shared_ptr<stk::middle_mesh::mesh::Mesh> get_mesh2() const override { return m_inputMesh2; }

  private:
    void create_input_meshes(MPI_Comm comm1, MPI_Comm comm2);

    std::shared_ptr<stk::middle_mesh::mesh::Mesh> m_inputMesh1;
    std::shared_ptr<stk::middle_mesh::mesh::Mesh> m_inputMesh2;
};


class ConservativeTransferTests
{
  public:
   ConservativeTransferTests(std::shared_ptr<stk::middle_mesh::mesh::Mesh> mesh1,
       std::shared_ptr<stk::middle_mesh::mesh::Mesh> mesh2,
       std::shared_ptr<ConservativeTransferUser> transferCallback1,
       std::shared_ptr<ConservativeTransferUser> transferCallback2,
       stk::middle_mesh::ApplicationInterfaceType type = stk::middle_mesh::ApplicationInterfaceType::FakeParallel);

   void test_exactness(std::function<double(const utils::Point&)> func);

   void test_exactness_vector(std::function<utils::Point(const utils::Point&)> func);

   void test_exactness_bidirectional(
       std::function<double(const utils::Point&)> func1, std::function<double(const utils::Point&)> func2);

   void test_conservation(std::function<double(const utils::Point&)> func);

   void test_conservation_bidirectional(
       std::function<double(const utils::Point&)> func1, std::function<double(const utils::Point&)> func2);

   std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh1;
   std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh2;
   std::shared_ptr<ConservativeTransferUser> transferCallback1;
   std::shared_ptr<ConservativeTransferUser> transferCallback2;
   std::shared_ptr<stk::transfer::ConservativeTransfer> conservativeTransfer;


  private:

    void set_field(stk::middle_mesh::mesh::FieldPtr<double> fieldPtr, std::function<double(const utils::Point&)> func);

    void set_field(stk::middle_mesh::mesh::FieldPtr<double> fieldPtr, std::function<utils::Point(const utils::Point&)> func);

};

}
}

#endif