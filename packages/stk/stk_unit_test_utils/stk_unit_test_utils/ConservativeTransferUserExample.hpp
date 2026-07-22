#ifndef STK_INTEGRATION_TESTS_CONSERVATIVE_TRANSFER_USER_EXAMPLE
#define STK_INTEGRATION_TESTS_CONSERVATIVE_TRANSFER_USER_EXAMPLE

#include <array>
#include <vector>
#include <utility>

#include "stk_middle_mesh/matrix.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_transfer/ConservativeTransfer.hpp"

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

    std::pair<double, double> get_xi_coord_range(mesh::MeshEntityType /*type*/) override { return std::make_pair(0.0, 1.0); }

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

    const std::vector<utils::Point>& get_quad_points([[maybe_unused]] mesh::MeshEntityType type) const
    {
      assert(type == mesh::MeshEntityType::Triangle);
      return m_quadPts;
    }

    const std::vector<double>& get_quad_weights([[maybe_unused]] mesh::MeshEntityType type) const
    {
      assert(type == mesh::MeshEntityType::Triangle);
      return m_quadWeights;
    }

    std::shared_ptr<XiCoordinatesFE> create_xi_coords() { return std::make_shared<XiCoordinatesFE>(m_quadPts); }

  private:
    std::vector<utils::Point> m_quadPts = { utils::Point(1.0/6, 1.0/6), utils::Point(2.0/3.0, 1.0/6), utils::Point(1.0/6, 2.0/3.0)};
    std::vector<double> m_quadWeights   = {1.0/3, 1.0/3, 1.0/3};
};

class ConservativeTransferUserForTest : public stk::transfer::ConservativeTransferUser
{
  public:
    ConservativeTransferUserForTest(std::shared_ptr<mesh::Mesh> inputMesh) :
      m_inputMesh(inputMesh),
      m_mat(mesh::count_valid(inputMesh->get_vertices()), mesh::count_valid(inputMesh->get_vertices()))
    {
      if (utils::impl::comm_size(inputMesh->get_comm()) != 1)
        throw std::runtime_error("mock ConservativeTransferUser only works for a 1 process mesh");
        
      do_dof_numbering();
    }

    std::shared_ptr<XiCoordinates> create_xi_coords() override { return m_finiteElement.create_xi_coords(); }

    void set_middle_mesh(std::shared_ptr<mesh::Mesh> middleMesh,
                         mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inputMeshToMiddleMesh    ,
                         mesh::FieldPtr<utils::Point> middleMeshQuadPointOnInputMesh) override
    {
      m_middleMesh = middleMesh;
      m_inputMeshToMiddleMesh = inputMeshToMiddleMesh;
      m_middleMeshQuadPointsOnInputMesh = middleMeshQuadPointOnInputMesh;
    }

    void interpolate_to_quad_pts(const mesh::FieldPtr<double> functionSolValsPtr, mesh::FieldPtr<double> functionQuadValsPtr) override
    {
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);
      assert(functionQuadValsPtr->get_mesh() == m_middleMesh);
      assert(functionSolValsPtr->get_num_comp() == functionQuadValsPtr->get_num_comp());

      int numQuadPointsPerElement          = functionQuadValsPtr->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *m_inputMeshToMiddleMesh;
      auto& middleMeshQuadPointOnInputMesh = *m_middleMeshQuadPointsOnInputMesh;
      auto& solVals  = *functionSolValsPtr;
      auto& quadVals = *functionQuadValsPtr;
      quadVals.set(0);
      std::array<double, mesh::MAX_DOWN> basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      for (auto& elInput : m_inputMesh->get_elements())
        if (elInput)
        {
          int nverts = mesh::get_downward(elInput, 0, elInputVerts.data());
          for (int i=0; i < inputMeshToMiddleMesh.get_num_comp(elInput, 0); ++i)
          {
            mesh::MeshEntityPtr elMiddle = inputMeshToMiddleMesh(elInput, 0, i);
            for (int j=0; j < numQuadPointsPerElement; ++j)
            {
              utils::Point quadPtXi = middleMeshQuadPointOnInputMesh(elMiddle, j, 0);
              m_finiteElement.get_basis_vals(elInput->get_type(), quadPtXi, basisVals);
              for (int k=0; k < nverts; ++k)
              {
                for (int d=0; d < solVals.get_num_comp(); ++d)
                  quadVals(elMiddle, j, d) += basisVals[k] * solVals(elInputVerts[k], 0, d);
              }
            }
          }
        }
    }

    void finish_integration(const mesh::FieldPtr<double> functionQuadValsPtr, mesh::FieldPtr<double> functionSolValsPtr) override
    {
      assert(functionQuadValsPtr->get_mesh() == m_middleMesh);
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);
      assert(functionQuadValsPtr->get_num_comp() == functionSolValsPtr->get_num_comp());

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
                for (int d=0; d < quadVals.get_num_comp(); ++d)
                  solVals(elInputVerts[k], 0, d) += quadVals(elMiddle, j, d) * basisVals[k] * weight * detJacobian;
            }
          }
        }   
    }

    void solve_linear_system(const mesh::FieldPtr<double> linearSystemRhsPtr, mesh::FieldPtr<double> functionSolValsPtr) override
    {                             
      assert(linearSystemRhsPtr->get_mesh() == m_inputMesh);
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);
      assert(linearSystemRhsPtr->get_num_comp() == functionSolValsPtr->get_num_comp());

      utils::impl::Matrix<double> massMatrix(m_numDofs, m_numDofs);
      compute_mass_matrix(massMatrix);

      std::vector<double> rhs(m_numDofs * functionSolValsPtr->get_num_comp());
      auto& dofNums = *m_dofNums;
      auto& linearSystemRhs = *linearSystemRhsPtr;
      for (auto& vert : m_inputMesh->get_vertices())
        if (vert)
        {
          int dof = dofNums(vert, 0, 0);
          for (int d=0; d < linearSystemRhs.get_num_comp(); ++d)
            rhs[dof + d * m_numDofs] = linearSystemRhs(vert, 0, d);
        }

      std::vector<int> ipiv(m_numDofs);
      utils::impl::solve_linear_system(massMatrix, ipiv.data(), rhs.data(), linearSystemRhs.get_num_comp());

      auto& solVals = *functionSolValsPtr;
      for (auto& vert : m_inputMesh->get_vertices())
        if (vert)
        {
          int dof = dofNums(vert, 0, 0);
          for (int d=0; d < solVals.get_num_comp(); ++d)
          solVals(vert, 0, d) = rhs[dof + d * m_numDofs];
        }
    }

    // only needed for testing
    std::vector<double> integrate_function(const mesh::FieldPtr<double> functionSolValsPtr) override
    {
      assert(functionSolValsPtr->get_mesh() == m_inputMesh);

      int numQuadPointsPerElement          = m_middleMeshQuadPointsOnInputMesh->get_field_shape().count[2];
      auto& inputMeshToMiddleMesh          = *m_inputMeshToMiddleMesh;
      auto& middleMeshQuadPointOnInputMesh = *m_middleMeshQuadPointsOnInputMesh;
      auto& solVals  =  *functionSolValsPtr;
      std::array<double, mesh::MAX_DOWN> elInputVals, basisVals;
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elInputVerts;
      std::vector<double> integralVals(solVals.get_num_comp(), 0);
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

              for (int k=0; k < nverts; ++k)
                for (int d=0; d < solVals.get_num_comp(); ++d)
                  integralVals[d] += basisVals[k] * weight * detJacobian * solVals(elInputVerts[k], 0, d);
            }
          }
        }

      return integralVals;
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

#endif