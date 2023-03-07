#include "application_interface_fake_parallel_spmd.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

ApplicationInterfaceFakeParallelSPMD::ApplicationInterfaceFakeParallelSPMD(
    std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2, const ParallelSearchOpts& parallelSearchOpts,
    const VolumeSnapOpts& volumeSnapOpts, const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
    const MiddleGridOpts& middleGridOpts,
    std::shared_ptr<XiCoordinates> xiPts)
  : m_comm(utils::impl::comm_dup(mesh1->get_comm()))
  , m_interface(mesh1, mesh2, m_comm, parallelSearchOpts, volumeSnapOpts, boundarySnapOpts, middleGridOpts, xiPts)
{
  if (!mesh1 && !mesh2)
    throw std::runtime_error("for ApplicationInterfaceSPMD, both meshes must be non-null on all processes");

  if (!(mesh1->get_comm() == mesh2->get_comm()))
    throw std::runtime_error("for ApplicationInterfaceSPMD, meshes must be defined on same communicator");
}

ApplicationInterfaceFakeParallelSPMD::~ApplicationInterfaceFakeParallelSPMD()
{
  MPI_Comm_free(&m_comm);
}

void ApplicationInterfaceFakeParallelSPMD::create_middle_grid()
{
  m_interface.create_middle_grid();
  m_middleGridCreated = true;
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelSPMD::get_middle_grid_for_mesh1()
{
  if (m_middleGridCreated)
    return m_interface.get_middle_grid_for_mesh1();
  else
    throw std::runtime_error("must create middle grid before retrieving it");
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelSPMD::get_middle_grid_for_mesh2()
{
  if (m_middleGridCreated)
    return m_interface.get_middle_grid_for_mesh2();
  else
    throw std::runtime_error("must create middle grid before retrieving it");
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelSPMD::get_mesh1_classification()
{
  if (m_middleGridCreated)
    return m_interface.get_mesh1_classification();
  else
    throw std::runtime_error("must create middle grid before retrieving information about it");
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelSPMD::get_mesh2_classification()
{
  if (m_middleGridCreated)
    return m_interface.get_mesh2_classification();
  else
    throw std::runtime_error("must create middle grid before retrieving information about it");
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
ApplicationInterfaceFakeParallelSPMD::compute_mesh1_inverse_classification()
{
  if (m_middleGridCreated)
    return m_interface.compute_mesh1_inverse_classification();
  else
    throw std::runtime_error("must create middle grid before retrieving information about it");
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
ApplicationInterfaceFakeParallelSPMD::compute_mesh2_inverse_classification()
{
  if (m_middleGridCreated)
    return m_interface.compute_mesh2_inverse_classification();
  else
    throw std::runtime_error("must create middle grid before retrieving information about it");
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceFakeParallelSPMD::get_remote_info_mesh_one_to_two()
{
  if (m_middleGridCreated)
  {
    return m_interface.get_remote_info_mesh_one_to_two();
  } else
    throw std::runtime_error("must create middle grid before retrieving information about it");  
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceFakeParallelSPMD::get_remote_info_mesh_two_to_one()
{
  if (m_middleGridCreated)
  {
    return m_interface.get_remote_info_mesh_two_to_one();
  } else
    throw std::runtime_error("must create middle grid before retrieving information about it");  
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceFakeParallelSPMD::get_xi_points_on_mesh1()
{
  if (m_middleGridCreated)
  {
    return m_interface.get_xi_points_on_mesh1();
  } else
    throw std::runtime_error("must create middle grid before retrieving information about it");  
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceFakeParallelSPMD::get_xi_points_on_mesh2()
{
  if (m_middleGridCreated)
  {
    return m_interface.get_xi_points_on_mesh2();
  } else
    throw std::runtime_error("must create middle grid before retrieving information about it");  
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
