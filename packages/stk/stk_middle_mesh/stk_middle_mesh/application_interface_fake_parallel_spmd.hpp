#ifndef APPLICATION_INTERFACE_FAKE_PARALLEL_SPMD_H
#define APPLICATION_INTERFACE_FAKE_PARALLEL_SPMD_H

#include "application_interface_fake_parallel_impl.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

class ApplicationInterfaceFakeParallelSPMD : public ApplicationInterfaceSPMD
{
  public:
    ApplicationInterfaceFakeParallelSPMD(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                         const ParallelSearchOpts& parallelSearchOpts,
                                         const VolumeSnapOpts& volumeSnapOpts,
                                         const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                                         const MiddleGridOpts& middleGridOpts,
                                         std::shared_ptr<XiCoordinates> xiPts);

    ~ApplicationInterfaceFakeParallelSPMD();

    void create_middle_grid() override;

    std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh1() override;

    std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh2() override;

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh1_classification() override;

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh2_classification() override;

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_inverse_classification() override;

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_inverse_classification() override;

    mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_one_to_two() override;

    mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_two_to_one() override;

    mesh::FieldPtr<utils::Point> get_xi_points_on_mesh1() override;
    
    mesh::FieldPtr<utils::Point> get_xi_points_on_mesh2() override;

  private:
    MPI_Comm m_comm;
    impl::ApplicationInterfaceFakeParallelImpl m_interface;
    bool m_middleGridCreated = false;
};

} // namespace impl
} // namespace middle_mesh
} // namespace stk

#endif