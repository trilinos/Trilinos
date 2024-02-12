#ifndef STK_TRANSFER_CONSERVATIVE_TRANSFER
#define STK_TRANSFER_CONSERVATIVE_TRANSFER

#include "stk_middle_mesh/application_interface.hpp"
#include "stk_middle_mesh/communication_api.hpp"
#include "stk_transfer/ConservativeTransferUser.hpp"

namespace stk {
namespace transfer {

class ConservativeTransfer
{
  public:
    ConservativeTransfer(MPI_Comm unionComm, std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh1,
                         std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh2,
                         std::shared_ptr<ConservativeTransferUser> transferCallback1,
                         std::shared_ptr<ConservativeTransferUser> transferCallback2,
                         stk::middle_mesh::ApplicationInterfaceType interfaceType = stk::middle_mesh::ApplicationInterfaceType::FakeParallel);

    ConservativeTransfer(MPI_Comm unionComm, std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh1,
                         std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh2,
                         std::shared_ptr<stk::middle_mesh::ApplicationInterface> interface,
                         std::shared_ptr<ConservativeTransferUser> transferCallback1,
                         std::shared_ptr<ConservativeTransferUser> transferCallback2);

    ConservativeTransfer(const ConservativeTransfer& rhs) = delete;

    ConservativeTransfer& operator=(const ConservativeTransfer& rhs) = delete;

    void start_transfer(stk::middle_mesh::mesh::FieldPtr<double> functionValsSend,
                        stk::middle_mesh::mesh::FieldPtr<double> functionValsRecv);

    void finish_transfer();

  private:

    void setup_for_transfer(std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh1,
                            std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh2,
                            std::shared_ptr<stk::middle_mesh::ApplicationInterface> interface,                           
                            std::shared_ptr<ConservativeTransferUser> transferCallback1,
                            std::shared_ptr<ConservativeTransferUser> transferCallback2);

    struct MeshData
    {
      std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh;
      std::shared_ptr<ConservativeTransferUser> transferCallback;
      std::shared_ptr<stk::middle_mesh::mesh::Mesh> middleMesh;
      stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::mesh::MeshEntityPtr> middleMeshToInputMesh;
      stk::middle_mesh::mesh::VariableSizeFieldPtr<stk::middle_mesh::mesh::MeshEntityPtr> inputMeshToMiddleMesh;
      stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::mesh::RemoteSharedEntity> remoteInfo;

      stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::utils::Point> middleMeshQuadPointsOnInputMesh;
      stk::middle_mesh::mesh::FieldPtr<double> quadVals;
    };


    MeshData get_mesh_data(stk::middle_mesh::mesh::FieldPtr<double> functionValsSend);

    void check_xi_coords_same_on_all_procs_debug_only(std::shared_ptr<ConservativeTransferUser> transferCallback1,
                                                      std::shared_ptr<ConservativeTransferUser> transferCallback2);

    int check_number_of_points(std::shared_ptr<stk::middle_mesh::XiCoordinates> xiCoords, middle_mesh::mesh::MeshEntityType type);

    void check_point_coords(std::shared_ptr<stk::middle_mesh::XiCoordinates> xiCoords, middle_mesh::mesh::MeshEntityType type,
                            int npts);

    MPI_Comm m_unionComm;
    std::shared_ptr<stk::middle_mesh::MiddleMeshFieldCommunication<double>> m_exchanger;
    MeshData m_mesh1Data;
    MeshData m_mesh2Data;
    
    stk::middle_mesh::mesh::FieldPtr<double> m_functionValsRecv;
    bool m_transferInProgress = false;
};

}
}


#endif