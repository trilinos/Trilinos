#include "application_interface.hpp"
#include "application_interface_fake_parallel_impl.hpp"
#include "application_interface_parallel.hpp"

#include <memory>
#include <stdexcept>

namespace stk {
namespace middle_mesh {


std::shared_ptr<ApplicationInterface> application_interface_factory(
    ApplicationInterfaceType type, std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
    MPI_Comm unionComm,
    std::shared_ptr<XiCoordinates> xiPts,   
    const ParallelSearchOpts& parallelSearchOpts, const VolumeSnapOpts& volumeSnapOpts,
    const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts, const MiddleGridOpts& middleGridOpts)
{
  switch (type)
  {
    case ApplicationInterfaceType::FakeParallel: {
      return std::make_shared<impl::ApplicationInterfaceFakeParallelImpl>(
          mesh1, mesh2, unionComm, parallelSearchOpts, 
          volumeSnapOpts, boundarySnapOpts, middleGridOpts, xiPts);
    }

    case ApplicationInterfaceType::Parallel:
    {
      return std::make_shared<impl::ApplicationInterfaceParallel>(
          mesh1, mesh2, unionComm, parallelSearchOpts, 
          volumeSnapOpts, boundarySnapOpts, middleGridOpts, xiPts);      
    }
    default:
      throw std::runtime_error("unhandled ApplicationInterfaceType enum");
  }
}

} // namespace middle_mesh
} // namespace stk
