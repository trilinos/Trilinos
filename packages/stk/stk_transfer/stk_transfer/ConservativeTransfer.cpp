#include "ConservativeTransfer.hpp"

namespace stk {
namespace transfer {

using namespace stk::middle_mesh;

ConservativeTransfer::ConservativeTransfer(MPI_Comm unionComm,
                                           std::shared_ptr<mesh::Mesh> inputMesh1,
                                           std::shared_ptr<mesh::Mesh> inputMesh2,
                                           std::shared_ptr<ConservativeTransferUser> transferCallback1,
                                           std::shared_ptr<ConservativeTransferUser> transferCallback2,
                                           stk::middle_mesh::ApplicationInterfaceType interfaceType) :
  m_unionComm(unionComm)
{
  if (unionComm == MPI_COMM_NULL)
    throw std::runtime_error("unionComm cannot be null");

  if ((inputMesh1 && !transferCallback1) || (transferCallback1 && !inputMesh1))
    throw std::runtime_error("must provide both inputMesh1 and transferCallback1, or neither");

  if ((inputMesh2 && !transferCallback2) || (transferCallback2 && !inputMesh2))
    throw std::runtime_error("must provide both inputMesh2 and transferCallback2, or neither");    

  std::shared_ptr<XiCoordinates> xiCoords;
  if (transferCallback1)
    xiCoords = transferCallback1->create_xi_coords();
  else if (transferCallback2)
    xiCoords = transferCallback2->create_xi_coords();

  std::shared_ptr<ApplicationInterface> interface = application_interface_factory(
    interfaceType, inputMesh1, inputMesh2, m_unionComm, xiCoords);  

  setup_for_transfer(inputMesh1, inputMesh2, interface, transferCallback1, transferCallback2);
}

ConservativeTransfer::ConservativeTransfer(MPI_Comm unionComm, std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh1,
                      std::shared_ptr<stk::middle_mesh::mesh::Mesh> inputMesh2,
                      std::shared_ptr<stk::middle_mesh::ApplicationInterface> interface,
                      std::shared_ptr<ConservativeTransferUser> transferCallback1,
                      std::shared_ptr<ConservativeTransferUser> transferCallback2) :
  m_unionComm(unionComm)
{
  if (unionComm == MPI_COMM_NULL)
    throw std::runtime_error("unionComm cannot be null");

  if ((inputMesh1 && !transferCallback1) || (transferCallback1 && !inputMesh1))
    throw std::runtime_error("must provide both inputMesh1 and transferCallback1, or neither");

  if ((inputMesh2 && !transferCallback2) || (transferCallback2 && !inputMesh2))
    throw std::runtime_error("must provide both inputMesh2 and transferCallback2, or neither");

  setup_for_transfer(inputMesh1, inputMesh2, interface, transferCallback1, transferCallback2);  
}                      

void ConservativeTransfer::setup_for_transfer(std::shared_ptr<mesh::Mesh> inputMesh1,
                                              std::shared_ptr<mesh::Mesh> inputMesh2,
                                              std::shared_ptr<stk::middle_mesh::ApplicationInterface> interface,
                                              std::shared_ptr<ConservativeTransferUser> transferCallback1,
                                              std::shared_ptr<ConservativeTransferUser> transferCallback2)
{
  check_xi_coords_same_on_all_procs_debug_only(transferCallback1, transferCallback2);

  interface->create_middle_grid();
  if (inputMesh1)
  {
    assert(transferCallback1);
    m_mesh1Data.inputMesh                        = inputMesh1;
    m_mesh1Data.transferCallback                 = transferCallback1;
    m_mesh1Data.middleMesh                       = interface->get_middle_grid_for_mesh1();
    m_mesh1Data.middleMeshToInputMesh            = interface->get_mesh1_classification();
    m_mesh1Data.inputMeshToMiddleMesh            = interface->compute_mesh1_inverse_classification();
    m_mesh1Data.remoteInfo                       = interface->get_remote_info_mesh_one_to_two();
    m_mesh1Data.middleMeshQuadPointsOnInputMesh  = interface->get_xi_points_on_mesh1();

    m_mesh1Data.transferCallback->set_middle_mesh(m_mesh1Data.middleMesh, m_mesh1Data.inputMeshToMiddleMesh,
                                                  m_mesh1Data.middleMeshQuadPointsOnInputMesh);
  }

  if (inputMesh2)
  {
    assert(transferCallback2);
    m_mesh2Data.inputMesh                        = inputMesh2;
    m_mesh2Data.transferCallback                 = transferCallback2;
    m_mesh2Data.middleMesh                       = interface->get_middle_grid_for_mesh2();
    m_mesh2Data.middleMeshToInputMesh            = interface->get_mesh2_classification();
    m_mesh2Data.inputMeshToMiddleMesh            = interface->compute_mesh2_inverse_classification();
    m_mesh2Data.remoteInfo                       = interface->get_remote_info_mesh_two_to_one();
    m_mesh2Data.middleMeshQuadPointsOnInputMesh  = interface->get_xi_points_on_mesh2();

    m_mesh2Data.transferCallback->set_middle_mesh(m_mesh2Data.middleMesh, m_mesh2Data.inputMeshToMiddleMesh,
                                                  m_mesh2Data.middleMeshQuadPointsOnInputMesh);
  }

  m_exchanger = std::make_shared<MiddleMeshFieldCommunication<double>>(m_unionComm, m_mesh1Data.middleMesh, m_mesh2Data.middleMesh,
                                                                       m_mesh1Data.remoteInfo, m_mesh2Data.remoteInfo);   
}


void ConservativeTransfer::start_transfer(mesh::FieldPtr<double> functionValsSend, mesh::FieldPtr<double> functionValsRecv)
{
  if (m_transferInProgress)
    throw std::runtime_error("cannot start a new transfer while another transfer is still in progress.  Did you forget to call finish_transfer()?");

  m_transferInProgress = true;
  m_functionValsRecv   = functionValsRecv;
  MeshData sendData    = get_mesh_data(functionValsSend);
  MeshData recvData    = get_mesh_data(functionValsRecv);

  if (functionValsSend)      
  {
    sendData.transferCallback->interpolate_to_quad_pts(functionValsSend, sendData.quadVals);
  }

  m_exchanger->start_exchange(sendData.quadVals, recvData.quadVals);
}

void ConservativeTransfer::finish_transfer()
{
  if (!m_transferInProgress)
    throw std::runtime_error("cannot finish a transfer before it is started");

  mesh::FieldPtr<double> functionValsRecv = m_functionValsRecv;
  MeshData recvData = get_mesh_data(functionValsRecv);

  m_exchanger->finish_exchange(recvData.quadVals);

  if (functionValsRecv)
  {
    mesh::FieldPtr<double> linearSystemRhs = mesh::create_field<double>(recvData.inputMesh, functionValsRecv->get_field_shape(),
                                                                        functionValsRecv->get_num_comp());
    recvData.transferCallback->finish_integration(recvData.quadVals, linearSystemRhs);
    recvData.transferCallback->solve_linear_system(linearSystemRhs, functionValsRecv);
  }

  m_functionValsRecv   = nullptr;
  m_transferInProgress = false;
}

ConservativeTransfer::MeshData ConservativeTransfer::get_mesh_data(mesh::FieldPtr<double> functionVals)
{
  MeshData data;
  if (functionVals)
  {
    bool onMesh1 = functionVals->get_mesh() == m_mesh1Data.inputMesh;
    data = onMesh1 ? m_mesh1Data : m_mesh2Data;   
    auto& quadVals = data.quadVals;

    if (!quadVals || quadVals->get_num_comp() != functionVals->get_num_comp())
    {
      quadVals = mesh::create_field<double>(data.middleMesh, data.middleMeshQuadPointsOnInputMesh->get_field_shape(),
                                            functionVals->get_num_comp());
    }
  }

  return data;
}


void ConservativeTransfer::check_xi_coords_same_on_all_procs_debug_only(
                            std::shared_ptr<ConservativeTransferUser> transferCallback1,
                            std::shared_ptr<ConservativeTransferUser> transferCallback2)
{
  std::shared_ptr<XiCoordinates> xiCoords;
  if (transferCallback1)
    xiCoords = transferCallback1->create_xi_coords();
  else if (transferCallback2)
    xiCoords = transferCallback2->create_xi_coords();

  std::array<mesh::MeshEntityType, 2> types = {mesh::MeshEntityType::Triangle, mesh::MeshEntityType::Quad};
  for (mesh::MeshEntityType type : types)
  {
    int npts = check_number_of_points(xiCoords, type);
    check_point_coords(xiCoords, type, npts);
  }
}

int ConservativeTransfer::check_number_of_points(std::shared_ptr<stk::middle_mesh::XiCoordinates> xiCoords,
                                                 middle_mesh::mesh::MeshEntityType type)
{
  const int rootRank = 0;
  const int commRank = utils::impl::comm_rank(m_unionComm);
  const int commSize = utils::impl::comm_size(m_unionComm);

  int npts = xiCoords ? static_cast<int>(xiCoords->get_xi_coords(type).size()) : -1;
  std::vector<int> nptsGathered;
  if (commRank == rootRank)
  {
    nptsGathered.resize(commSize);
  }

  MPI_Gather(&npts, 1, MPI_INT, nptsGathered.data(), 1, MPI_INT, rootRank, m_unionComm);

  int nptsExpected = -1;

  if (commRank == rootRank)
  {
    for (int i=0; i < commSize; ++i)
    {
      if (nptsGathered[i] != -1)
      {
        if (nptsExpected == -1)
        {
          nptsExpected = nptsGathered[i];
        }

        if (nptsExpected != nptsGathered[i])
        {
          std::stringstream ss;
          ss << "number of quadrature points on " << type << " is inconsistent: some processes have "
              << nptsExpected << " while others have " << nptsGathered[i] << std::endl;
          throw std::runtime_error(ss.str());
        }
      }   
    }

    if (nptsExpected == -1)
    {
      throw std::runtime_error("the XiCoords object was not provided by any process");
    }
  }

  MPI_Bcast(&nptsExpected, 1, MPI_INT, rootRank, m_unionComm);
  return nptsExpected;
}

void ConservativeTransfer::check_point_coords(std::shared_ptr<stk::middle_mesh::XiCoordinates> xiCoords,
                                              middle_mesh::mesh::MeshEntityType type, int npts)
{
  const int rootRank = 0;
  const int commRank = utils::impl::comm_rank(m_unionComm);
  const int commSize = utils::impl::comm_size(m_unionComm);
  const double tol = 1e-12;

  std::vector<middle_mesh::utils::Point> pts(npts);
  int doIHaveXiCoords = -1;
  if (xiCoords)
  {
    assert(xiCoords->get_xi_coords(type).size() == size_t(npts));
    pts = xiCoords->get_xi_coords(type);
    doIHaveXiCoords = 1;
  }

  std::vector<middle_mesh::utils::Point> ptsGathered;
  std::vector<int> mask;
  if (commRank == rootRank)
  {
    ptsGathered.resize(npts * commSize);
    mask.resize(commSize);
  }

  MPI_Gather(pts.data(), 3*npts, MPI_DOUBLE, ptsGathered.data(), 3*npts, MPI_DOUBLE, rootRank, m_unionComm);
  MPI_Gather(&doIHaveXiCoords, 1, MPI_INT, mask.data(), 1, MPI_INT, rootRank, m_unionComm);

  if (commRank == rootRank)
  {
    bool foundPtsExpected = false;
    std::vector<middle_mesh::utils::Point> ptsExpected(npts);
    for (int i=0; i < commSize; ++i)
      if (mask[i] == 1)
      {
        if (!foundPtsExpected)
        {
          for (int j=0; j < npts; ++j)
            ptsExpected[j] = ptsGathered[npts * i + j];

          foundPtsExpected = true;
        }

        for (int j=0; j < npts; ++j)
        {
          middle_mesh::utils::Point disp = ptsGathered[npts * i + j] - ptsExpected[j];
          double dist = std::sqrt(dot(disp, disp));
          if (dist > tol)
          {
            std::stringstream ss;
            ss << "quadrature point location on " << type << " point " << j << " is inconsistent: some processes have "
                << ptsGathered[npts * i + j] << " while others have " << ptsExpected[j] << ", distance between them is "
                << dist << " > " << tol << std::endl;
            throw std::runtime_error(ss.str());
          }
        }
      }
  }
}


}
}