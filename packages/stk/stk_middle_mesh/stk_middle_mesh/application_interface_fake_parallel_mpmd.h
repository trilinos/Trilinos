#ifndef APPLICATION_INTERFACE_FAKE_PARALLEL_MPMD
#define APPLICATION_INTERFACE_FAKE_PARALLEL_MPMD

#include "application_interface.h"
#include "application_interface_fake_parallel_impl.h"

namespace stk {
namespace middle_mesh {
namespace impl {

class ApplicationInterfaceFakeParallelMPMD : public ApplicationInterfaceMPMD
{
  public:
    ApplicationInterfaceFakeParallelMPMD(std::shared_ptr<mesh::Mesh> mesh, bool isMesh1, MPI_Comm unionComm,
                                         const ParallelSearchOpts& parallelSearchOpts,
                                         const VolumeSnapOpts& volumeSnapOpts,
                                         const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                                         const MiddleGridOpts& middleGridOpts);

    ~ApplicationInterfaceFakeParallelMPMD();

    std::shared_ptr<mesh::Mesh> create_middle_grid() override;

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh_classification() override;

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh_inverse_classification() override;

  private:
    void check_union_comm_size();

    std::shared_ptr<mesh::Mesh> m_mesh1Parallel;
    std::shared_ptr<mesh::Mesh> m_mesh2Parallel;
    MPI_Comm m_unionComm;
    impl::ApplicationInterfaceFakeParallelImpl m_interface;
};

} // namespace impl
} // namespace middle_mesh
} // namespace stk

#endif