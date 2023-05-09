#include "stk_middle_mesh/application_interface.hpp"
#include "stk_middle_mesh/communication_api.hpp"
#include "stk_middle_mesh/matrix.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>

#include "stk_middle_mesh/mesh_io.hpp"

#include <array>
#include <vector>
#include <utility>

using namespace stk::middle_mesh;


class XiCoordinatesFE : public XiCoordinates
{
  public:
    XiCoordinatesFE(const std::vector<utils::Point>& triangleQuadPoints) :
      m_triangleIntegrationPoints(triangleQuadPoints)
    {}

    const std::vector<utils::Point>& get_xi_coords(mesh::MeshEntityType type) override
    {
      if (type == mesh::MeshEntityType::Triangle)
        return m_triangleIntegrationPoints;
      else
        return m_quadIntegrationPoints;
    }

    std::pair<double, double> get_xi_coord_range(mesh::MeshEntityType type) override { return std::make_pair(0.0, 1.0); }

  private:
    std::vector<utils::Point> m_triangleIntegrationPoints;
    std::vector<utils::Point> m_quadIntegrationPoints;
};


class FiniteElement
{
  public:
    void get_basis_vals(mesh::MeshEntityType type, const utils::Point& ptXi, std::array<double, mesh::MAX_DOWN>& vals)
    {
      if (type == mesh::MeshEntityType::Quad)
      {
        std::array<double, 2> xiVals, etaVals;
        mesh::compute_lagrange_vals(ptXi.x, xiVals.data());
        mesh::compute_lagrange_vals(ptXi.y, etaVals.data());

        vals[0] = xiVals[0]*etaVals[0];
        vals[1] = xiVals[1]*etaVals[0];
        vals[2] = xiVals[1]*etaVals[1];
        vals[3] = xiVals[0]*etaVals[1];
      } else
        throw std::runtime_error("unsupported element type");
    };

    const std::vector<utils::Point>& get_quad_points(mesh::MeshEntityType type) const
    {
      assert(type == mesh::MeshEntityType::Triangle);
      return m_quadPts;
    }

    const std::vector<double>& get_quad_weights(mesh::MeshEntityType type) const
    {
      assert(type == mesh::MeshEntityType::Triangle);
      return m_quadWeights;
    }

    std::shared_ptr<XiCoordinatesFE> create_xi_coords() { return std::make_shared<XiCoordinatesFE>(m_quadPts); }

  private:
    std::vector<utils::Point> m_quadPts = { utils::Point(1.0/6, 1.0/6), utils::Point(2.0/3.0, 1.0/6), utils::Point(1.0/6, 2.0/3.0)};
    std::vector<double> m_quadWeights   = {1.0/3, 1.0/3, 1.0/3};
};

class ConservativeTransferUser
{
  public:
    ConservativeTransferUser(std::shared_ptr<mesh::Mesh> inputMesh) :
      m_inputMesh(inputMesh),
      m_mat(mesh::count_valid(inputMesh->get_vertices()), mesh::count_valid(inputMesh->get_vertices()))
    {
      do_dof_numbering();
    }

    std::shared_ptr<XiCoordinatesFE> create_xi_coords() { return m_finiteElement.create_xi_coords(); }

    void set_middle_mesh(std::shared_ptr<mesh::Mesh> middleMesh,
                         mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMesh    ,
                         mesh::FieldPtr<utils::Point> middleMeshQuadPointOnInputMesh)
    {
      m_middleMesh = middleMesh;
      m_inputMeshToMiddleMesh = inputMeshToMiddleMesh;
      m_middleMeshQuadPointsOnInputMesh = middleMeshQuadPointOnInputMesh;
    }

    void interpolate_to_quad_pts(mesh::FieldPtr<double> functionSolValsPtr, mesh::FieldPtr<double> functionQuadValsPtr)
    {
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);
      assert(functionQuadValsPtr->get_mesh() == m_middleMesh);

      int numQuadPointsPerElement          = functionQuadValsPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *m_inputMeshToMiddleMesh;
      auto& middleMeshQuadPointOnInputMesh = *m_middleMeshQuadPointsOnInputMesh;
      auto& solVals  =  *functionSolValsPtr;
      auto& quadVals = *functionQuadValsPtr;
      std::array<double, mesh::MAX_DOWN> elInputVals, basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      for (auto& elInput : m_inputMesh->get_elements())
        if (elInput)
        {
          int nverts = mesh::get_downward(elInput, 0, elInputVerts.data());
          for (int k=0; k < nverts; ++k)
            elInputVals[k] = solVals(elInputVerts[k], 0, 0);

          for (int i=0; i < inputMeshToMiddleMesh.get_num_comp(elInput, 0); ++i)
          {
            mesh::MeshEntityPtr elMiddle = inputMeshToMiddleMesh(elInput, 0, i);
            for (int j=0; j < numQuadPointsPerElement; ++j)
            {
              utils::Point quadPtXi = middleMeshQuadPointOnInputMesh(elMiddle, j, 0);
              m_finiteElement.get_basis_vals(elInput->get_type(), quadPtXi, basisVals);
              double val_j = 0;
              for (int k=0; k < nverts; ++k)
                val_j += basisVals[k] * elInputVals[k];

              quadVals(elMiddle, j, 0) = val_j;
            }
          }
        }
    }

    void finish_integration(mesh::FieldPtr<double> functionQuadValsPtr, mesh::FieldPtr<double> functionSolValsPtr)
    {
      assert(functionQuadValsPtr->get_mesh() == m_middleMesh);
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);

      int numQuadPointsPerElement          = functionQuadValsPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *m_inputMeshToMiddleMesh;
      auto& middleMeshQuadPointOnInputMesh = *m_middleMeshQuadPointsOnInputMesh;
      auto& quadVals = *functionQuadValsPtr;
      auto& solVals  =  *functionSolValsPtr;
      solVals.set(0);
      std::array<double, mesh::MAX_DOWN> basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      for (auto& elInput : m_inputMesh->get_elements())
        if (elInput)
        {
          int nverts = mesh::get_downward(elInput, 0, elInputVerts.data());
          for (int i=0; i < inputMeshToMiddleMesh.get_num_comp(elInput, 0); ++i)
          {
            mesh::MeshEntityPtr elMiddle = inputMeshToMiddleMesh(elInput, 0, i);
            double detJacobian = compute_triangle_area(elMiddle);

            for (int j=0; j < numQuadPointsPerElement; ++j)
            {
              utils::Point quadPtXiOnInputMesh = middleMeshQuadPointOnInputMesh(elMiddle, j, 0);
              m_finiteElement.get_basis_vals(elInput->get_type(), quadPtXiOnInputMesh, basisVals);              
              double weight = m_finiteElement.get_quad_weights(elMiddle->get_type())[j];

              for (int k=0; k < nverts; ++k)
                solVals(elInputVerts[k], 0, 0) += quadVals(elMiddle, j, 0) * basisVals[k] * weight * detJacobian;
            }
          }
        }   
    }

    void solve_linear_system(mesh::FieldPtr<double> linearSystemRhsPtr, mesh::FieldPtr<double> functionSolValsPtr)
    {                             
      assert(linearSystemRhsPtr->get_mesh() == m_inputMesh);
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);

      utils::impl::Matrix<double> massMatrix(m_numDofs, m_numDofs);
      compute_mass_matrix(massMatrix);

      std::vector<double> rhs(m_numDofs);
      auto& dofNums = *m_dofNums;
      auto& linearSystemRhs = *linearSystemRhsPtr;
      for (auto& vert : m_inputMesh->get_vertices())
        if (vert)
        {
          int dof = dofNums(vert, 0, 0);
          rhs[dof] = linearSystemRhs(vert, 0, 0);
        }

      std::vector<int> ipiv(m_numDofs);
      utils::impl::solve_linear_system(massMatrix, ipiv.data(), rhs.data());

      auto& solVals = *functionSolValsPtr;
      for (auto& vert : m_inputMesh->get_vertices())
        if (vert)
        {
          int dof = dofNums(vert, 0, 0);
          solVals(vert, 0, 0) = rhs[dof];
        }
    }

    // only needed for testing
    double integrate_function(mesh::FieldPtr<double> functionSolValsPtr)
    {
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);

      int numQuadPointsPerElement          = m_middleMeshQuadPointsOnInputMesh->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *m_inputMeshToMiddleMesh;
      auto& middleMeshQuadPointOnInputMesh = *m_middleMeshQuadPointsOnInputMesh;
      auto& solVals  =  *functionSolValsPtr;
      std::array<double, mesh::MAX_DOWN> elInputVals, basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      double integralVal = 0;
      for (auto& elInput : m_inputMesh->get_elements())
        if (elInput)
        {
          int nverts = mesh::get_downward(elInput, 0, elInputVerts.data());
          for (int k=0; k < nverts; ++k)
            elInputVals[k] = solVals(elInputVerts[k], 0, 0);

          for (int i=0; i < inputMeshToMiddleMesh.get_num_comp(elInput, 0); ++i)
          {
            mesh::MeshEntityPtr elMiddle = inputMeshToMiddleMesh(elInput, 0, i);
            double detJacobian = compute_triangle_area(elMiddle);

            for (int j=0; j < numQuadPointsPerElement; ++j)
            {
              utils::Point quadPtXi = middleMeshQuadPointOnInputMesh(elMiddle, j, 0);
              m_finiteElement.get_basis_vals(elInput->get_type(), quadPtXi, basisVals);
              double weight = m_finiteElement.get_quad_weights(elMiddle->get_type())[j];

              double val_j = 0;
              for (int k=0; k < nverts; ++k)
                val_j += basisVals[k] * elInputVals[k];

              integralVal += val_j * weight * detJacobian;
            }
          }
        }

      return integralVal;
    }


  private:

    void compute_mass_matrix(utils::impl::Matrix<double>& massMatrix)
    {
      //TODO: you could also do the integration on the inputMesh rather than the middle mesh, but that would make
      //      it harder to experiment with different geometric terms

      massMatrix.fill(0);

      int numQuadPointsPerElement          = m_middleMeshQuadPointsOnInputMesh->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *m_inputMeshToMiddleMesh;
      auto& middleMeshQuadPointOnInputMesh = *m_middleMeshQuadPointsOnInputMesh;
      std::array<double, mesh::MAX_DOWN> basisVals;
      utils::impl::Matrix<double> elementMatrix(mesh::MAX_DOWN, mesh::MAX_DOWN);
      for (auto& elInput : m_inputMesh->get_elements())
        if (elInput)
        {
          int nverts = elInput->get_type() == mesh::MeshEntityType::Quad ? 4 : 3;
          elementMatrix.fill(0);

          for (int i=0; i < inputMeshToMiddleMesh.get_num_comp(elInput, 0); ++i)
          {
            mesh::MeshEntityPtr elMiddle = inputMeshToMiddleMesh(elInput, 0, i);
            for (int j=0; j < numQuadPointsPerElement; ++j)
            {
              utils::Point quadPtXiOnInputMesh = middleMeshQuadPointOnInputMesh(elMiddle, j, 0);
              m_finiteElement.get_basis_vals(elInput->get_type(), quadPtXiOnInputMesh, basisVals);              
              double weight = m_finiteElement.get_quad_weights(elMiddle->get_type())[j];
              double detJacobian = compute_triangle_area(elMiddle);

              for (int k=0; k < nverts; ++k)
                for (int p=0; p < nverts; ++p)
                  elementMatrix(k, p) += basisVals[k] * basisVals[p] * weight * detJacobian;
            }
          }

          std::array<int, mesh::MAX_DOWN> dofs = get_element_dofs(elInput);
          for (int k=0; k < nverts; ++k)
            for (int p=0; p < nverts; ++p)
              massMatrix(dofs[k], dofs[p]) += elementMatrix(k, p);
        } 
    } 

    std::array<int, mesh::MAX_DOWN> get_element_dofs(mesh::MeshEntityPtr el)
    {
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
      int nverts = mesh::get_downward(el, 0, verts.data());

      std::array<int, mesh::MAX_DOWN> dofs;
      auto& dofNums = *m_dofNums;
      for (int i=0; i < nverts; ++i)
        dofs[i] = dofNums(verts[i], 0, 0);

      return dofs;        
    }

    double compute_triangle_area(mesh::MeshEntityPtr tri)
    {
      assert(tri->get_type() == mesh::MeshEntityType::Triangle);
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
      mesh::get_downward(tri, 0, verts.data());

      utils::Point b1 = verts[1]->get_point_orig(0) - verts[0]->get_point_orig(0);
      utils::Point b2 = verts[2]->get_point_orig(0) - verts[0]->get_point_orig(0);
      utils::Point normal = cross(b1, b2);
      return std::sqrt(dot(normal, normal))/2;
    }

    void do_dof_numbering()
    {
      m_dofNums = mesh::create_field<int>(m_inputMesh, mesh::FieldShape(1, 0, 0), 1);

      // who cares about matrix bandwith anyways
      auto& dofNums = *m_dofNums;
      int dof = 0;
      for (auto v : m_inputMesh->get_vertices())
        if (v)
        {
          dofNums(v, 0, 0) = dof++;
        }

      m_numDofs = dof;
    }


    std::shared_ptr<mesh::Mesh> m_inputMesh;
    mesh::FieldPtr<int> m_dofNums;
    int m_numDofs = 0;
    utils::impl::Matrix<double> m_mat;
    FiniteElement m_finiteElement;


    std::shared_ptr<mesh::Mesh> m_middleMesh;
    mesh::FieldPtr<utils::Point> m_middleMeshQuadPointsOnInputMesh;
	  mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_inputMeshToMiddleMesh;   
};


class ConservativeTransferMPMDStk
{
  public:
    ConservativeTransferMPMDStk(MPI_Comm unionComm, std::shared_ptr<mesh::Mesh> inputMesh1, std::shared_ptr<mesh::Mesh> inputMesh2,
                                std::shared_ptr<ConservativeTransferUser> transferCallback1, std::shared_ptr<ConservativeTransferUser> transferCallback2) :
      m_unionComm(unionComm),
      m_inputMesh1(inputMesh1),
      m_inputMesh2(inputMesh2),
      m_transferCallback1(transferCallback1),
      m_transferCallback2(transferCallback2)
    {
      setup_for_transfer();
    }

    void setup_for_transfer()
    {
      std::shared_ptr<XiCoordinates> xiCoords;
      if (m_transferCallback1)
        xiCoords = m_transferCallback1->create_xi_coords();
      else
        xiCoords = m_transferCallback2->create_xi_coords();

	    std::shared_ptr<ApplicationInterface> interface = application_interface_factory(
	      ApplicationInterfaceType::FakeParallel, m_inputMesh1, m_inputMesh2, m_unionComm, xiCoords);

      interface->create_middle_grid();
      if (m_inputMesh1)
      {
        m_middleMesh1                      = interface->get_middle_grid_for_mesh1();
        m_middleMeshToInputMesh1Info       = interface->get_mesh1_classification();
        m_inputMesh1ToMiddleMesh           = interface->compute_mesh1_inverse_classification();
        m_remoteInfo1To2                   = interface->get_remote_info_mesh_one_to_two();
        m_middleMeshQuadPointsOnInputMesh1 = interface->get_xi_points_on_mesh1();

        m_transferCallback1->set_middle_mesh(m_middleMesh1, m_inputMesh1ToMiddleMesh, m_middleMeshQuadPointsOnInputMesh1);
      }

      if (m_inputMesh2)
      {
        m_middleMesh2                      = interface->get_middle_grid_for_mesh2();
        m_middleMeshToInputMesh2Info       = interface->get_mesh2_classification();
        m_inputMesh2ToMiddleMesh           = interface->compute_mesh2_inverse_classification();
        m_remoteInfo2To1                   = interface->get_remote_info_mesh_two_to_one();
        m_middleMeshQuadPointsOnInputMesh2 = interface->get_xi_points_on_mesh2();

        m_transferCallback2->set_middle_mesh(m_middleMesh2, m_inputMesh2ToMiddleMesh, m_middleMeshQuadPointsOnInputMesh2);
      }

      m_exchanger = std::make_shared<MiddleMeshFieldCommunication<double>>(m_unionComm, m_middleMesh1, m_middleMesh2, 
                                                                           m_remoteInfo1To2, m_remoteInfo2To1);   
    }

    void start_transfer(mesh::FieldPtr<double> functionValsSend, mesh::FieldPtr<double> functionValsRecv)
    {
      MeshData sendData = get_mesh_data_send(functionValsSend);
      MeshData recvData = get_mesh_data_recv(functionValsRecv);

      if (functionValsSend)      
      {
        sendData.transferCallback->interpolate_to_quad_pts(functionValsSend, sendData.quadVals);
      }

	    m_exchanger->start_exchange(sendData.quadVals, recvData.quadVals);
    }

    void finish_transfer(mesh::FieldPtr<double> functionValsSend, mesh::FieldPtr<double> functionValsRecv)
    {
      MeshData sendData = get_mesh_data_send(functionValsSend);
      MeshData recvData = get_mesh_data_recv(functionValsRecv);

	    m_exchanger->finish_exchange(recvData.quadVals);

      if (functionValsRecv)
      {
        mesh::FieldPtr<double> linearSystemRhs = mesh::create_field<double>(recvData.inputMesh, mesh::FieldShape(1, 0, 0), 1);
        recvData.transferCallback->finish_integration(recvData.quadVals, linearSystemRhs);
        recvData.transferCallback->solve_linear_system(linearSystemRhs, functionValsRecv);
      }

      // for testing only
      check_conservation(functionValsSend, functionValsRecv);
    }

  private:

    struct MeshData
    {
      std::shared_ptr<mesh::Mesh> inputMesh;
      std::shared_ptr<ConservativeTransferUser> transferCallback;
      std::shared_ptr<mesh::Mesh> middleMesh;
      mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMesh;
      mesh::FieldPtr<utils::Point> middleMeshQuadPointsOnInputMesh;
      mesh::FieldPtr<double> quadVals;
    };


    MeshData get_mesh_data_send(mesh::FieldPtr<double> functionValsSend)
    {
      MeshData data;
      if (functionValsSend)
      {
        bool sendOneToTwo = functionValsSend->get_mesh() == m_inputMesh1;

        data.inputMesh  = sendOneToTwo ? m_inputMesh1 : m_inputMesh2;
        data.middleMesh = sendOneToTwo ? m_middleMesh1 : m_middleMesh2;
        data.transferCallback       = sendOneToTwo ? m_transferCallback1 : m_transferCallback2;
        data.middleMeshQuadPointsOnInputMesh = sendOneToTwo ? m_middleMeshQuadPointsOnInputMesh1 
                                                            : m_middleMeshQuadPointsOnInputMesh2;
        data.inputMeshToMiddleMesh = sendOneToTwo ? m_inputMesh1ToMiddleMesh : m_inputMesh2ToMiddleMesh;
        auto& quadVals = sendOneToTwo ? m_quadValsMiddleMesh1 : m_quadValsMiddleMesh2;

        if (!quadVals || quadVals->get_num_comp() != functionValsSend->get_num_comp())
        {
          quadVals = mesh::create_field<double>(data.middleMesh, data.middleMeshQuadPointsOnInputMesh->get_field_shape(),
                                                functionValsSend->get_num_comp());
        }

        data.quadVals = quadVals;
      }

      return data;
    }


    MeshData get_mesh_data_recv(mesh::FieldPtr<double> functionValsRecv)
    {
      MeshData data;
      if (functionValsRecv)
      {
        bool sendOneToTwo = functionValsRecv->get_mesh() == m_inputMesh2;

        data.inputMesh  = sendOneToTwo ? m_inputMesh2 : m_inputMesh1;
        data.middleMesh = sendOneToTwo ? m_middleMesh2 : m_middleMesh1;
        data.transferCallback       = sendOneToTwo ? m_transferCallback2 : m_transferCallback1;
        data.middleMeshQuadPointsOnInputMesh = sendOneToTwo ? m_middleMeshQuadPointsOnInputMesh2 
                                                            : m_middleMeshQuadPointsOnInputMesh1;
        data.inputMeshToMiddleMesh = sendOneToTwo ? m_inputMesh2ToMiddleMesh : m_inputMesh1ToMiddleMesh;

        auto& quadVals = sendOneToTwo ? m_quadValsMiddleMesh2 : m_quadValsMiddleMesh1;
        quadVals = mesh::create_field<double>(data.middleMesh, data.middleMeshQuadPointsOnInputMesh->get_field_shape(),
                                              functionValsRecv->get_num_comp());

        if (!quadVals || quadVals->get_num_comp() != functionValsRecv->get_num_comp())
        {
          quadVals = mesh::create_field<double>(data.middleMesh, data.middleMeshQuadPointsOnInputMesh->get_field_shape(),
                                                functionValsRecv->get_num_comp());
        }

        data.quadVals = quadVals;        
      }

      return data;
    }   

    void check_conservation(mesh::FieldPtr<double> functionValsSend, mesh::FieldPtr<double> functionValsRecv)
    {
      MeshData sendData = get_mesh_data_send(functionValsSend);
      MeshData recvData = get_mesh_data_recv(functionValsRecv);
  
      int myRank = utils::impl::comm_rank(m_unionComm);
      if (functionValsSend)
      {
        double integralVal = sendData.transferCallback->integrate_function(functionValsSend);
        MPI_Send(&integralVal, 1, MPI_DOUBLE, 1 - myRank, 666, m_unionComm);
      } else
      {
        double integralVal = recvData.transferCallback->integrate_function(functionValsRecv);        
        double integralValSender;
        MPI_Recv(&integralValSender, 1, MPI_DOUBLE, 1 - myRank, 666, m_unionComm, MPI_STATUS_IGNORE);

        std::cout << "integralValSender = " << integralValSender << ", integralValRecver = " << integralVal << ", diff = " << std::abs(integralValSender - integralVal) << std::endl;
        if (std::abs(integralValSender - integralVal) > 1e-12)
          throw std::runtime_error("transfer was not conservative");
      }
    }

    MPI_Comm m_unionComm;
    std::shared_ptr<mesh::Mesh> m_inputMesh1;
    std::shared_ptr<mesh::Mesh> m_inputMesh2;
    std::shared_ptr<ConservativeTransferUser> m_transferCallback1;
    std::shared_ptr<ConservativeTransferUser> m_transferCallback2;
    std::shared_ptr<MiddleMeshFieldCommunication<double>> m_exchanger;


    std::shared_ptr<mesh::Mesh> m_middleMesh1;
    std::shared_ptr<mesh::Mesh> m_middleMesh2;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_middleMeshToInputMesh1Info;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_middleMeshToInputMesh2Info;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_inputMesh1ToMiddleMesh;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_inputMesh2ToMiddleMesh;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo1To2;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo2To1;
    mesh::FieldPtr<utils::Point> m_middleMeshQuadPointsOnInputMesh1;
    mesh::FieldPtr<utils::Point> m_middleMeshQuadPointsOnInputMesh2;
    mesh::FieldPtr<double> m_quadValsMiddleMesh1;
    mesh::FieldPtr<double> m_quadValsMiddleMesh2;    
};

template <typename Tfunc>
mesh::FieldPtr<double> create_field(std::shared_ptr<mesh::Mesh> mesh, Tfunc func)
{
  auto fieldPtr = mesh::create_field<double>(mesh, mesh::FieldShape(1, 0, 0), 1);

  auto& field = *fieldPtr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      field(vert, 0, 0) = func(pt);
    }

  return fieldPtr;
}


template <typename Tfunc>
void check_field(mesh::FieldPtr<double> fieldPtr, Tfunc func)
{
  auto& field = *fieldPtr;
  for (auto& vert : field.get_mesh()->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      double valExpected = func(pt);
      //std::cout << "field = " << field(vert, 0, 0) << ", val expected = " << valExpected << ", error = " << std::abs(field(vert, 0, 0) - valExpected) << std::endl;
      if (std::abs(field(vert, 0, 0) - valExpected) > 1e-12)
        throw std::runtime_error("field transfer was not exact");
    }
}


std::function<double(const utils::Point&)> function_factory(const std::string& functionName)
{
  if (functionName == "linear")
  {
    return [](const utils::Point& pt) { return pt.x + 2*pt.y + 3*pt.z; };
  } else if (functionName == "quadratic")
  {
    return [](const utils::Point& pt) { return pt.x*pt.x + 2*pt.y*pt.y + 3*pt.z; };
  } else if (functionName == "exponential")
  {
    return [](const utils::Point& pt) { return std::exp(pt.x + pt.y + pt.z ); };
  } else
    throw std::runtime_error("unrecognized function name: " + functionName);
}

int main(int argc, char* argv[])
{
  stk::parallel_machine_init(&argc, &argv);

  int defaultColor = -1;
  int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
  STK_ThrowRequireMsg(color == 0 || color == 1, "app-color command line argument must be provided, and must either 0 or 1");

  int sendColor = stk::get_command_line_option(argc, argv, "send-color", defaultColor);
  STK_ThrowRequireMsg(sendColor == 0 || sendColor == 1, "send-color command line argument must be provided, and must either 0 or 1");

  std::string defaultFileName;
  if (color == 0)
    defaultFileName = "generated:3x3x1|sideset:Z|bbox:0,0,0,1,1,1";
  else
    defaultFileName = "generated:4x4x1|sideset:z|bbox:0,0,1,1,1,2";

  std::string meshFileName = stk::get_command_line_option(argc, argv, "mesh", defaultFileName);

  std::string defaultPartName = "surface_1";
  std::string partName = stk::get_command_line_option(argc, argv, "part-name", defaultPartName);

  std::string defaultFunctionName = "linear";
  std::string functionName = stk::get_command_line_option(argc, argv, "function-name", defaultFunctionName);

  {
    MPI_Comm meshComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

    stk_interface::impl::StkMeshCreator creator(meshFileName, meshComm);
    std::shared_ptr<mesh::Mesh> inputMesh = creator.create_mesh_from_part(partName).mesh;

    mesh::impl::print_vert_edges("mesh" + std::to_string(color), inputMesh);
    std::shared_ptr<mesh::Mesh> inputMesh1 = color == 0 ? inputMesh : nullptr;
    std::shared_ptr<mesh::Mesh> inputMesh2 = color == 0 ? nullptr   : inputMesh;

    auto transferCallback = std::make_shared<ConservativeTransferUser>(inputMesh);
    std::shared_ptr<ConservativeTransferUser> transferCallback1 = color == 0 ? transferCallback : nullptr;
    std::shared_ptr<ConservativeTransferUser> transferCallback2 = color == 0 ? nullptr : transferCallback;
    if (color == 0)
      std::cout << "mesh1 numel = " << inputMesh1->get_elements().size() << std::endl;

    if (color == 1)
      std::cout << "mesh2 numel = " << inputMesh2->get_elements().size() << std::endl;    
      
    MPI_Comm unionComm = MPI_COMM_WORLD;
    ConservativeTransferMPMDStk transfer(unionComm, inputMesh1, inputMesh2, transferCallback1, transferCallback2);

    bool amISender = color == sendColor;
    auto func = function_factory(functionName); // [](const utils::Point& pt) { return pt.x*pt.x + 2*pt.y*pt.y + 3*pt.z; };
    mesh::FieldPtr<double> functionValsSend, functionValsRecv;
    if (amISender) {
      functionValsSend = create_field(inputMesh, func);
      functionValsRecv = nullptr;
    } else {
      functionValsSend = nullptr;
      functionValsRecv = mesh::create_field<double>(inputMesh, mesh::FieldShape(1, 0, 0), 1);
    }

    
    transfer.start_transfer(functionValsSend, functionValsRecv);
    transfer.finish_transfer(functionValsSend, functionValsRecv);

    if (!amISender && functionName == "linear")
      check_field(functionValsRecv, func);

    MPI_Comm_free(&meshComm);
  }

  stk::parallel_machine_finalize();
}