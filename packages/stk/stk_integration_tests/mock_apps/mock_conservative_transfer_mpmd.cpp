#include "stk_middle_mesh/application_interface.hpp"
#include "stk_middle_mesh/communication_api_mpmd.hpp"
#include "stk_middle_mesh/matrix.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>


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

    void interpolate_to_quad_pts(std::shared_ptr<mesh::Mesh> inputMesh, std::shared_ptr<mesh::Mesh> middleMesh,
                               mesh::FieldPtr<utils::Point> middleMeshQuadPointOnInputMeshPtr,
	                             mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMeshPtr,
							                 mesh::FieldPtr<double> functionSolValsPtr, mesh::FieldPtr<double> functionQuadValsPtr)
    {
      assert(inputMesh == m_inputMesh);
      assert(functionSolValsPtr->get_mesh() == inputMesh);
      assert(functionQuadValsPtr->get_mesh() == middleMesh);

      int numQuadPointsPerElement          = functionQuadValsPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *inputMeshToMiddleMeshPtr;
      auto& middleMeshQuadPointOnInputMesh = *middleMeshQuadPointOnInputMeshPtr;
      auto& solVals  =  *functionSolValsPtr;
      auto& quadVals = *functionQuadValsPtr;
      std::array<double, mesh::MAX_DOWN> elInputVals, basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      for (auto& elInput : inputMesh->get_elements())
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

    void finish_integration(std::shared_ptr<mesh::Mesh> inputMesh, std::shared_ptr<mesh::Mesh> middleMesh,
                            mesh::FieldPtr<utils::Point> middleMeshQuadPointsOnInputMeshPtr,
	                          mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMeshPtr,
					                  mesh::FieldPtr<double> functionQuadValsPtr, mesh::FieldPtr<double> functionSolValsPtr)
    {
      assert(inputMesh == m_inputMesh);
      assert(functionQuadValsPtr->get_mesh() == middleMesh);
      assert(functionSolValsPtr->get_mesh() == inputMesh);

      int numQuadPointsPerElement          = functionQuadValsPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *inputMeshToMiddleMeshPtr;
      auto& middleMeshQuadPointOnInputMesh = *middleMeshQuadPointsOnInputMeshPtr;
      auto& quadVals = *functionQuadValsPtr;
      auto& solVals  =  *functionSolValsPtr;
      solVals.set(0);
      std::array<double, mesh::MAX_DOWN> basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      for (auto& elInput : inputMesh->get_elements())
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

    void solve_linear_system(std::shared_ptr<mesh::Mesh> inputMesh, std::shared_ptr<mesh::Mesh> middleMesh,
                             mesh::FieldPtr<utils::Point> middleMeshQuadPointsOnInputMeshPtr,
	                           mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMeshPtr,
							               mesh::FieldPtr<double> linearSystemRhsPtr, mesh::FieldPtr<double> functionSolValsPtr)
    {                             
      assert(inputMesh == m_inputMesh);

      utils::impl::Matrix<double> massMatrix(m_numDofs, m_numDofs);
      compute_mass_matrix(inputMesh, middleMesh, middleMeshQuadPointsOnInputMeshPtr, inputMeshToMiddleMeshPtr, massMatrix);

      std::vector<double> rhs(m_numDofs);
      auto& dofNums = *m_dofNums;
      auto& linearSystemRhs = *linearSystemRhsPtr;
      for (auto& vert : inputMesh->get_vertices())
        if (vert)
        {
          int dof = dofNums(vert, 0, 0);
          rhs[dof] = linearSystemRhs(vert, 0, 0);
        }

      std::vector<int> ipiv(m_numDofs);
      utils::impl::solve_linear_system(massMatrix, ipiv.data(), rhs.data());

      auto& solVals = *functionSolValsPtr;
      for (auto& vert : inputMesh->get_vertices())
        if (vert)
        {
          int dof = dofNums(vert, 0, 0);
          solVals(vert, 0, 0) = rhs[dof];
        }
    }

    // only needed for testing
    double integrate_function(std::shared_ptr<mesh::Mesh> inputMesh, std::shared_ptr<mesh::Mesh> middleMesh,
                              mesh::FieldPtr<utils::Point> middleMeshQuadPointOnInputMeshPtr,
	                            mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMeshPtr,
							                mesh::FieldPtr<double> functionSolValsPtr)
    {
      assert(inputMesh == m_inputMesh);
      assert(functionSolValsPtr->get_mesh() == inputMesh);

      int numQuadPointsPerElement          = middleMeshQuadPointOnInputMeshPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *inputMeshToMiddleMeshPtr;
      auto& middleMeshQuadPointOnInputMesh = *middleMeshQuadPointOnInputMeshPtr;
      auto& solVals  =  *functionSolValsPtr;
      std::array<double, mesh::MAX_DOWN> elInputVals, basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      double integralVal = 0;
      for (auto& elInput : inputMesh->get_elements())
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

    void compute_mass_matrix(std::shared_ptr<mesh::Mesh> inputMesh, std::shared_ptr<mesh::Mesh> middleMesh,
                             mesh::FieldPtr<utils::Point> middleMeshQuadPointsOnInputMeshPtr,
	                           mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMeshPtr,
							               utils::impl::Matrix<double>& massMatrix)
    {
      //TODO: you could also do the integration on the inputMesh rather than the middle mesh, but that would make
      //      it harder to experiment with different geometric terms

      massMatrix.fill(0);

      int numQuadPointsPerElement          = middleMeshQuadPointsOnInputMeshPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *inputMeshToMiddleMeshPtr;
      auto& middleMeshQuadPointOnInputMesh = *middleMeshQuadPointsOnInputMeshPtr;
      std::array<double, mesh::MAX_DOWN> basisVals;
      utils::impl::Matrix<double> elementMatrix(mesh::MAX_DOWN, mesh::MAX_DOWN);
      for (auto& elInput : inputMesh->get_elements())
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
};


class ConservativeTransferMPMDStk
{
  public:
    ConservativeTransferMPMDStk(std::shared_ptr<mesh::Mesh> inputMesh, bool isProjectToMesh,
                                std::shared_ptr<ConservativeTransferUser> user) :
      m_inputMesh(inputMesh),
      m_isProjectToMesh(isProjectToMesh),
      m_user(user)
    {
      setup_for_transfer();
    }

    void setup_for_transfer()
    {
      m_unionComm = MPI_COMM_WORLD;
      std::shared_ptr<XiCoordinates> xiCoords = m_user->create_xi_coords();   
	    std::shared_ptr<ApplicationInterfaceMPMD> interface = application_interface_mpmd_factory(
	      ApplicationInterfaceType::FakeParallel, m_inputMesh, m_isProjectToMesh, m_unionComm, xiCoords);

      m_middleMesh                      = interface->create_middle_grid();
      m_middleMeshToInputMeshInfo       = interface->get_mesh_classification();
      m_inputMeshToMiddleMesh           = interface->compute_mesh_inverse_classification();
      m_remoteInfo                      = interface->get_remote_info();
      m_middleMeshQuadPointsOnInputMesh = interface->get_xi_points_on_input_mesh();
    }

    void do_transfer(mesh::FieldPtr<double> functionVals, bool amISender)
    {
      mesh::FieldPtr<double> quadVals = mesh::create_field<double>(m_middleMesh, m_middleMeshQuadPointsOnInputMesh->get_field_shape(), 1);
      if (amISender)
      {
        m_user->interpolate_to_quad_pts(m_inputMesh, m_middleMesh,
                                        m_middleMeshQuadPointsOnInputMesh, m_inputMeshToMiddleMesh,
                                        functionVals, quadVals);
      }

      MiddleMeshFieldCommunicationMPMD<double> exchanger(m_unionComm, m_middleMesh, m_remoteInfo);	  
	    exchanger.start_exchange(quadVals, amISender);
	    exchanger.finish_exchange(quadVals, amISender);

      if (!amISender)
      {
        mesh::FieldPtr<double> linearSystemRhs = mesh::create_field<double>(m_inputMesh, mesh::FieldShape(1, 0, 0), 1);
        m_user->finish_integration(m_inputMesh, m_middleMesh, m_middleMeshQuadPointsOnInputMesh, m_inputMeshToMiddleMesh,
                                   quadVals, linearSystemRhs);
        m_user->solve_linear_system(m_inputMesh, m_middleMesh, m_middleMeshQuadPointsOnInputMesh, m_inputMeshToMiddleMesh,
                                    linearSystemRhs, functionVals);
      }

      // for testing only
      check_conservation(functionVals, amISender);
    }

  private:

    void check_conservation(mesh::FieldPtr<double> functionSolValsPtr, bool amISender)
    {
      double integralVal = m_user->integrate_function(m_inputMesh, m_middleMesh, m_middleMeshQuadPointsOnInputMesh,
                                                      m_inputMeshToMiddleMesh, functionSolValsPtr);
      
      int myRank = utils::impl::comm_rank(m_unionComm);
      if (amISender)
      {
        MPI_Send(&integralVal, 1, MPI_DOUBLE, 1 - myRank, 666, m_unionComm);
      } else
      {
        double integralValSender;
        MPI_Recv(&integralValSender, 1, MPI_DOUBLE, 1 - myRank, 666, m_unionComm, MPI_STATUS_IGNORE);

        if (std::abs(integralValSender - integralVal) > 1e-12)
          throw std::runtime_error("transfer was not conservative");
      }
    }

    std::shared_ptr<mesh::Mesh> m_inputMesh;
    bool m_isProjectToMesh;
    std::shared_ptr<ConservativeTransferUser> m_user;

    MPI_Comm m_unionComm;
    std::shared_ptr<mesh::Mesh> m_middleMesh;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_middleMeshToInputMeshInfo;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_inputMeshToMiddleMesh;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfo;
    mesh::FieldPtr<utils::Point> m_middleMeshQuadPointsOnInputMesh;   
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

/*
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
*/

int main(int argc, char* argv[])
{
  stk::parallel_machine_init(&argc, &argv);

  int defaultColor = -1;
  int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
  ThrowRequireMsg(color == 0 || color == 1, "app-color command line argument must be provided, and must either 0 or 1");

  int sendColor = stk::get_command_line_option(argc, argv, "send-color", defaultColor);
  ThrowRequireMsg(sendColor == 0 || sendColor == 1, "send-color command line argument must be provided, and must either 0 or 1");


  {
    MPI_Comm meshComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

    mesh::impl::MeshSpec spec;
    if (color == 0)
    {
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;
      spec.numelX = 3;
      spec.numelY = 3;
    } else {
      spec.xmin   = 0;
      spec.xmax   = 1;
      spec.ymin   = 0;
      spec.ymax   = 1;
      spec.numelX = 4;
      spec.numelY = 4;    
    }

    auto f = [](const stk::middle_mesh::utils::Point& pt) { return pt; };
    std::shared_ptr<mesh::Mesh> inputMesh = stk::middle_mesh::mesh::impl::create_mesh(spec, f, meshComm);
    bool isProjectToMesh = color == 0;

    auto user = std::make_shared<ConservativeTransferUser>(inputMesh);
    ConservativeTransferMPMDStk transfer(inputMesh, isProjectToMesh, user);

    bool amISender = color == sendColor;
    auto func = [](const utils::Point& pt) { return pt.x*pt.x + 2*pt.y*pt.y + 3*pt.z; };
    mesh::FieldPtr<double> functionVals;
    if (amISender)
      functionVals = create_field(inputMesh, func);
    else
      functionVals = mesh::create_field<double>(inputMesh, mesh::FieldShape(1, 0, 0), 1);
    
    transfer.do_transfer(functionVals, amISender);

    //if (!amISender)
    //  check_field(functionVals, func);

    MPI_Comm_free(&meshComm);
  }

  stk::parallel_machine_finalize();
}