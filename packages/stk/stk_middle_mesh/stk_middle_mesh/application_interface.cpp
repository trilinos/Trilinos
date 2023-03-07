#include "application_interface.hpp"
#include "application_interface_fake_parallel_mpmd.hpp"
#include "application_interface_fake_parallel_spmd.hpp"

#include <memory>
#include <stdexcept>

namespace stk {
namespace middle_mesh {

std::shared_ptr<ApplicationInterfaceMPMD> application_interface_mpmd_factory(
    ApplicationInterfaceType type, std::shared_ptr<mesh::Mesh> mesh, bool isMesh1, MPI_Comm unionComm,
    std::shared_ptr<XiCoordinates> xiPts,
    const ParallelSearchOpts& parallelSearchOpts, const VolumeSnapOpts& volumeSnapOpts,
    const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts, const MiddleGridOpts& middleGridOpts)
{
  switch (type)
  {
    case ApplicationInterfaceType::FakeParallel: {
      return std::make_shared<impl::ApplicationInterfaceFakeParallelMPMD>(
          mesh, isMesh1, unionComm, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts, xiPts);
    }
    default:
      throw std::runtime_error("unhandled ApplicationInterfaceType enum");
  }
}

std::shared_ptr<ApplicationInterfaceSPMD> application_interface_spmd_factory(
    ApplicationInterfaceType type, std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
    std::shared_ptr<XiCoordinates> xiPts,   
    const ParallelSearchOpts& parallelSearchOpts, const VolumeSnapOpts& volumeSnapOpts,
    const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts, const MiddleGridOpts& middleGridOpts)
{
  switch (type)
  {
    case ApplicationInterfaceType::FakeParallel: {
      return std::make_shared<impl::ApplicationInterfaceFakeParallelSPMD>(
          mesh1, mesh2, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts, xiPts);
    }
    default:
      throw std::runtime_error("unhandled ApplicationInterfaceType enum");
  }
}

} // namespace middle_mesh
} // namespace stk
