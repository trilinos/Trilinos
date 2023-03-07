#include "application_interface_fake_parallel_mpmd.hpp"
#include "stk_util/parallel/MPITagManager.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

ApplicationInterfaceFakeParallelMPMD::ApplicationInterfaceFakeParallelMPMD(
    std::shared_ptr<mesh::Mesh> mesh, bool isMesh1, MPI_Comm unionComm, const ParallelSearchOpts& parallelSearchOpts,
    const VolumeSnapOpts& volumeSnapOpts, const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
    const MiddleGridOpts& middleGridOpts,
    std::shared_ptr<XiCoordinates> xiPts)
  : m_mesh1Parallel(isMesh1 ? mesh : nullptr)
  , m_mesh2Parallel(isMesh1 ? nullptr : mesh)
  , m_unionComm(utils::impl::comm_dup(unionComm))
  , m_interface(m_mesh1Parallel, m_mesh2Parallel, m_unionComm, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts,
                middleGridOpts, xiPts)
{
  if (!mesh)
    throw std::runtime_error("mesh cannot be nullptr");

  check_union_comm_size();
}

ApplicationInterfaceFakeParallelMPMD::~ApplicationInterfaceFakeParallelMPMD()
{
  MPI_Comm_free(&m_unionComm);
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelMPMD::create_middle_grid()
{
  m_interface.create_middle_grid();

  if (m_mesh1Parallel)
    return m_interface.get_middle_grid_for_mesh1();
  else
    return m_interface.get_middle_grid_for_mesh2();
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelMPMD::get_mesh_classification()
{
  if (m_mesh1Parallel)
    return m_interface.get_mesh1_classification();
  else
    return m_interface.get_mesh2_classification();
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
ApplicationInterfaceFakeParallelMPMD::compute_mesh_inverse_classification()
{
  if (m_mesh1Parallel)
    return m_interface.compute_mesh1_inverse_classification();
  else
    return m_interface.compute_mesh2_inverse_classification();
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceFakeParallelMPMD::get_remote_info()
{
  if (m_mesh1Parallel)
    return m_interface.get_remote_info_mesh_one_to_two();
  else
    return m_interface.get_remote_info_mesh_two_to_one();
}

mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_two_to_one();

void ApplicationInterfaceFakeParallelMPMD::check_union_comm_size()
{
  int rootRank   = 0;
  MPI_Comm comm1 = m_mesh1Parallel ? m_mesh1Parallel->get_comm() : MPI_COMM_NULL;
  MPI_Comm comm2 = m_mesh2Parallel ? m_mesh2Parallel->get_comm() : MPI_COMM_NULL;
  auto commSizes = impl::get_comm_sizes_on_root(comm1, comm2, m_unionComm, rootRank);

  if (utils::impl::comm_rank(m_unionComm) == rootRank)
  {
    if (utils::impl::comm_size(m_unionComm) < (commSizes.first + commSizes.second))
      throw std::runtime_error("union comm is too small to be the union of the two mesh communicators");
  }
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceFakeParallelMPMD::get_xi_points_on_input_mesh()
{
  if (m_mesh1Parallel)
    return m_interface.get_xi_points_on_mesh1();
  else
    return m_interface.get_xi_points_on_mesh2();
}


} // namespace impl
} // namespace middle_mesh
} // namespace stk
