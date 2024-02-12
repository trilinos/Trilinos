#ifndef STK_TRANSFER_CONSERVATIVE_TRANSFER_USER
#define STK_TRANSFER_CONSERVATIVE_TRANSFER_USER

#include "stk_middle_mesh/application_interface.hpp"

namespace stk {
namespace transfer{

class ConservativeTransferUser
{
  public:
    virtual ~ConservativeTransferUser() = default;

    virtual std::shared_ptr<stk::middle_mesh::XiCoordinates> create_xi_coords() = 0;

    virtual void set_middle_mesh(std::shared_ptr<stk::middle_mesh::mesh::Mesh> middleMesh,
                                 stk::middle_mesh::mesh::VariableSizeFieldPtr<stk::middle_mesh::mesh::MeshEntityPtr> inputMeshToMiddleMesh,
                                 stk::middle_mesh::mesh::FieldPtr<stk::middle_mesh::utils::Point> middleMeshQuadPointOnInputMesh) = 0;

    virtual void interpolate_to_quad_pts(const stk::middle_mesh::mesh::FieldPtr<double> functionSolValsPtr,
                                         stk::middle_mesh::mesh::FieldPtr<double> functionQuadValsPtr) = 0;

    virtual void finish_integration(const stk::middle_mesh::mesh::FieldPtr<double> functionQuadValsPtr,
                                    stk::middle_mesh::mesh::FieldPtr<double> functionSolValsPtr) = 0;

    virtual void solve_linear_system(const stk::middle_mesh::mesh::FieldPtr<double> linearSystemRhsPtr,
                                     stk::middle_mesh::mesh::FieldPtr<double> functionSolValsPtr) = 0;

    // only needed for testing
    virtual std::vector<double> integrate_function(stk::middle_mesh::mesh::FieldPtr<double> functionSolValsPtr) = 0;
};

}
}

#endif