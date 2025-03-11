#include "application_interface_fake_parallel_impl.hpp"

#include "gathered_mesh_coordinate_update.hpp"
#include "incremental_mesh_boundary_snapper.hpp"
#include "mesh_classification_scatter.hpp"
#include "mesh_gather_to_root.hpp"
#include "mesh_scatter_from_root.hpp"
#include "field_scatter_from_root.hpp"
#include "mesh_snapper.hpp"
#include "nonconformal4.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace impl {

void ApplicationInterfaceFakeParallelImpl::check_union_comm_size()
{
  MPI_Comm comm1 = m_mesh1Parallel ? m_mesh1Parallel->get_comm() : MPI_COMM_NULL;
  MPI_Comm comm2 = m_mesh2Parallel ? m_mesh2Parallel->get_comm() : MPI_COMM_NULL;
  auto commSizes = impl::get_comm_sizes_on_root(comm1, comm2, m_unionComm, m_rootRankOnUnionComm);

  if (utils::impl::comm_rank(m_unionComm) == m_rootRankOnUnionComm &&
      utils::impl::comm_size(m_unionComm) < std::max(commSizes.first, commSizes.second))
    throw std::runtime_error("union comm is too small to be the union of the two mesh communicators");
}

void ApplicationInterfaceFakeParallelImpl::check_both_meshes_exist_somewhere()
{
  int mesh1CountLocal  = m_mesh1Parallel ? 1 : 0;
  int mesh2CountLocal  = m_mesh2Parallel ? 1 : 0;
  int mesh1CountGlobal = 0, mesh2CountGlobal = 0;

  MPI_Allreduce(&mesh1CountLocal, &mesh1CountGlobal, 1, MPI_INT, MPI_SUM, m_unionComm);
  MPI_Allreduce(&mesh2CountLocal, &mesh2CountGlobal, 1, MPI_INT, MPI_SUM, m_unionComm);

  if (mesh1CountGlobal == 0)
    throw std::runtime_error("no process provided mesh1");

  if (mesh2CountGlobal == 0)
    throw std::runtime_error("no process provided mesh2");
}

void ApplicationInterfaceFakeParallelImpl::create_middle_grid()
{
  gather_meshes_to_root();

  do_volume_snap();

  do_boundary_snap();

  create_serial_middle_grid();

  scatter_meshes();

  m_middleGridCreated = true;
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelImpl::get_middle_grid_for_mesh1()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh1Parallel)
    return m_middleGridParallel1;
  else
    throw std::runtime_error(
        "attempting to retrieve middle grid for mesh1 even though mesh 1 is the nullptr on this process");
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelImpl::get_middle_grid_for_mesh2()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh2Parallel)
    return m_middleGridParallel2;
  else
    throw std::runtime_error(
        "attempting to retrieve middle grid for mesh 2 even though mesh 2 is the nullptr on this process");
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelImpl::get_mesh1_classification()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh1Parallel)
    return m_mesh1ClassificationParallel;
  else
    throw std::runtime_error(
        "attempting to retrieve mesh 1 classification even though mesh 1 is the nullptr on this process");
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelImpl::get_mesh2_classification()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh2Parallel)
    return m_mesh2ClassificationParallel;
  else
    throw std::runtime_error(
        "attempting to retrieve mesh 2 classification even though mesh 2 is the nullptr on this process");
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
ApplicationInterfaceFakeParallelImpl::compute_mesh1_inverse_classification()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh1Parallel)
    return nonconformal4::impl::invert_classification_field(m_middleGridParallel1, m_mesh1Parallel,
                                                            m_mesh1ClassificationParallel);
  else
    throw std::runtime_error(
        "attempting to compute mesh 1 inverse classification even though mesh 1 is the nullptr on this process");
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
ApplicationInterfaceFakeParallelImpl::compute_mesh2_inverse_classification()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh2Parallel)
    return nonconformal4::impl::invert_classification_field(m_middleGridParallel2, m_mesh2Parallel,
                                                            m_mesh2ClassificationParallel);
  else
    throw std::runtime_error(
        "attempting to compute mesh 2 inverse classification even though mesh 2 is the nullptr on this process");
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceFakeParallelImpl::get_remote_info_mesh_one_to_two()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh1Parallel)
  {
    return m_remoteInfoMeshOneToTwo;
  } else
      throw std::runtime_error(
        "attempting to retrieve middle mesh 1 to 2 remote info even though mesh 1 is the nullptr on this process");
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceFakeParallelImpl::get_remote_info_mesh_two_to_one()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh2Parallel)
  {
    return m_remoteInfoMeshTwoToOne;
  } else
      throw std::runtime_error(
        "attempting to retrieve middle mesh 2 to 1 remote info even though mesh 2 is the nullptr on this process");
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceFakeParallelImpl::get_xi_points_on_mesh1()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh1Parallel && m_xiPts)
  {
    return m_mesh1XiPointsParallel;
  } else
      throw std::runtime_error(
        "attempting to retrieve xi points on mesh 1 info even though mesh 1 is the nullptr on this process or no xi points were given");
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceFakeParallelImpl::get_xi_points_on_mesh2()
{
  STK_ThrowRequireMsg(m_middleGridCreated, "must call create_middle_grid() before accessing information about the middle grid");

  if (m_mesh2Parallel && m_xiPts)
  {
    return m_mesh2XiPointsParallel;
  } else
      throw std::runtime_error(
        "attempting to retrieve xi points on mesh 2 info even though mesh 2 is the nullptr on this process or no xi points were given");
}

void ApplicationInterfaceFakeParallelImpl::gather_meshes_to_root()
{
  mesh::impl::MeshGatherToRoot gather1(m_unionComm, m_rootRankOnUnionComm, m_mesh1Parallel);
  m_mesh1Serial         = gather1.gather();
  m_mesh1ElementOrigins = gather1.get_element_origins();

  mesh::impl::MeshGatherToRoot gather2(m_unionComm, m_rootRankOnUnionComm, m_mesh2Parallel);
  m_mesh2Serial         = gather2.gather();
  m_mesh2ElementOrigins = gather2.get_element_origins();

  if (utils::impl::comm_rank(m_unionComm) != m_rootRankOnUnionComm && m_mesh1Serial != nullptr)
    throw std::runtime_error("mesh1 serial should be null");

  if (utils::impl::comm_rank(m_unionComm) != m_rootRankOnUnionComm && m_mesh1Serial != nullptr)
    throw std::runtime_error("mesh1 serial should be null");
}

void ApplicationInterfaceFakeParallelImpl::do_volume_snap()
{
  if (!(m_mesh1Serial && m_mesh2Serial))
    return;

  if (m_volumeSnapOpts.type == VolumeSnapType::Standard)
  {
    mesh::impl::MeshSnapper snapper(m_volumeSnapOpts.standardOpts);
    snapper.snap(m_mesh1Serial, m_mesh2Serial);
  } else if (m_volumeSnapOpts.type != VolumeSnapType::None)
  {
    throw std::runtime_error("unhandled volume snap enum");
  }
}

void ApplicationInterfaceFakeParallelImpl::do_boundary_snap()
{
  if (!(m_mesh1Serial && m_mesh2Serial))
    return;

  if (m_boundarySnapOpts.type == BoundarySnapAndQualityImprovementType::SnapThenQuality)
  {
    mesh::impl::BoundaryFixture fixture1(m_mesh1Serial);
    auto qualityImprover1 =
        mesh::impl::make_standard_improver(m_mesh1Serial, fixture1, m_boundarySnapOpts.snapThenQualityOpts);

    mesh::impl::BoundaryFixture fixture2(m_mesh2Serial);
    auto qualityImprover2 =
        mesh::impl::make_standard_improver(m_mesh2Serial, fixture2, m_boundarySnapOpts.snapThenQualityOpts);

    mesh::impl::MeshBoundarySnapper snapper;
    snapper.snap(m_mesh1Serial, m_mesh2Serial, MPI_COMM_SELF);

    qualityImprover1->run();
    qualityImprover2->run();
  } else if (m_boundarySnapOpts.type == BoundarySnapAndQualityImprovementType::IncrementalBoundarySnap)
  {
    auto snapper = mesh::impl::make_incremental_boundary_snapper(m_mesh1Serial, m_mesh2Serial, MPI_COMM_SELF,
                                                                 m_boundarySnapOpts.incrementalMeshBoundarySnapOpts);
    snapper->snap();
  } else if (m_boundarySnapOpts.type != BoundarySnapAndQualityImprovementType::None)
  {
    throw std::runtime_error("unhandled boundary snap enum");
  }
}

void ApplicationInterfaceFakeParallelImpl::create_serial_middle_grid()
{
  if (!(m_mesh1Serial && m_mesh2Serial))
    return;

  if (m_middleGridOpts.type == MiddleGridType::NormalProjection)
  {
    nonconformal4::impl::Nonconformal4 nonconformal(m_mesh1Serial, m_mesh2Serial,
                                                    m_middleGridOpts.normalProjectionOpts);
    m_middleGridSerial          = nonconformal.create();
    m_mesh1ClassificationSerial = nonconformal.get_mesh1_classification();
    m_mesh2ClassificationSerial = nonconformal.get_mesh2_classification();
    if (m_xiPts)
    {
      m_mesh1XiPointsSerial = nonconformal.compute_points_on_mesh1(m_xiPts);
      m_mesh2XiPointsSerial = nonconformal.compute_points_on_mesh2(m_xiPts);
    }
  } else
  {
    throw std::runtime_error("unhandled middle grid enum");
  }
}

void ApplicationInterfaceFakeParallelImpl::scatter_meshes()
{
  scatter_mesh(m_mesh1ClassificationSerial, m_mesh1ElementOrigins, m_mesh1Parallel, m_middleGridSerialEntityDestinations1,
               m_middleGridParallel1, m_mesh1ClassificationParallel);
  scatter_mesh(m_mesh2ClassificationSerial, m_mesh2ElementOrigins, m_mesh2Parallel, m_middleGridSerialEntityDestinations2,
               m_middleGridParallel2, m_mesh2ClassificationParallel);

  mesh::impl::GatheredMeshCoordinateUpdate updator1(m_unionComm, m_mesh1Serial, m_mesh1ElementOrigins, m_mesh1Parallel);
  updator1.update();

  mesh::impl::GatheredMeshCoordinateUpdate updator2(m_unionComm, m_mesh2Serial, m_mesh2ElementOrigins, m_mesh2Parallel);
  updator2.update();

  scatter_remote_info();

  scatter_xi_points();

}

void ApplicationInterfaceFakeParallelImpl::scatter_mesh(mesh::FieldPtr<mesh::MeshEntityPtr> meshClassificationSerial,
                                                        mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins,
                                                        std::shared_ptr<mesh::Mesh> inputMeshParallel,
                                                        mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity>& middleGridSerialEntityDestinations,
                                                        std::shared_ptr<mesh::Mesh>& middleGridParallel,
                                                        mesh::FieldPtr<mesh::MeshEntityPtr>& meshClassificationParallel)
{
  auto elementDestinations = get_element_destinations(m_middleGridSerial, meshClassificationSerial, elementOrigins);
  MPI_Comm meshComm        = inputMeshParallel ? inputMeshParallel->get_comm() : MPI_COMM_NULL;
  mesh::impl::MeshScatterFromRoot scatter1(m_unionComm, m_middleGridSerial, meshComm, elementDestinations);
  middleGridParallel = scatter1.scatter();
  middleGridSerialEntityDestinations = scatter1.get_entity_destinations();

  mesh::impl::MeshClassificationScatter classificationScatter(
      m_unionComm, m_middleGridSerial, middleGridSerialEntityDestinations, middleGridParallel, 
      inputMeshParallel, elementOrigins, meshClassificationSerial);
  meshClassificationParallel = classificationScatter.scatter();
}


void ApplicationInterfaceFakeParallelImpl::scatter_remote_info()
{
  mesh::FieldPtr<mesh::RemoteSharedEntity> remoteInfoSerialOneToTwo, remoteInfoSerialTwoToOne;
  if (m_middleGridSerial)
  {
    remoteInfoSerialOneToTwo = mesh::create_field<mesh::RemoteSharedEntity>(m_middleGridSerial, mesh::impl::FieldShape(0, 0, 1), 1);
    remoteInfoSerialTwoToOne = mesh::create_field<mesh::RemoteSharedEntity>(m_middleGridSerial, mesh::impl::FieldShape(0, 0, 1), 1);

    auto& middleGridSerialEntityDestinations1 = *m_middleGridSerialEntityDestinations1;
    auto& middleGridSerialEntityDestinations2 = *m_middleGridSerialEntityDestinations2;
    for (auto& el : m_middleGridSerial->get_elements())
    {
      (*remoteInfoSerialOneToTwo)(el, 0, 0) = middleGridSerialEntityDestinations2(el, 0, 0);
      (*remoteInfoSerialTwoToOne)(el, 0, 0) = middleGridSerialEntityDestinations1(el, 0, 0);
    }
  }

  if (m_middleGridParallel1)
  {
    m_remoteInfoMeshOneToTwo = mesh::create_field<mesh::RemoteSharedEntity>(m_middleGridParallel1, mesh::impl::FieldShape(0, 0, 1), 1);
  }

  mesh::impl::FieldScatterFromRoot<mesh::RemoteSharedEntity> scatterer1(m_unionComm, m_rootRankOnUnionComm, 
                            m_middleGridSerialEntityDestinations1, remoteInfoSerialOneToTwo, m_remoteInfoMeshOneToTwo);
  scatterer1.scatter();

  if (m_middleGridParallel2)
  {
    m_remoteInfoMeshTwoToOne = mesh::create_field<mesh::RemoteSharedEntity>(m_middleGridParallel2, mesh::impl::FieldShape(0, 0, 1), 1);
  }

  mesh::impl::FieldScatterFromRoot<mesh::RemoteSharedEntity> scatterer2(m_unionComm, m_rootRankOnUnionComm, 
                            m_middleGridSerialEntityDestinations2, remoteInfoSerialTwoToOne, m_remoteInfoMeshTwoToOne);
  scatterer2.scatter();    
}

void ApplicationInterfaceFakeParallelImpl::scatter_xi_points()
{

  scatter_xi_points(m_middleGridParallel1, m_middleGridSerialEntityDestinations1, m_mesh1XiPointsSerial, m_mesh1XiPointsParallel);
  scatter_xi_points(m_middleGridParallel2, m_middleGridSerialEntityDestinations2, m_mesh2XiPointsSerial, m_mesh2XiPointsParallel);

}

void ApplicationInterfaceFakeParallelImpl::scatter_xi_points(std::shared_ptr<mesh::Mesh> middleMeshParallel,
                                                             mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityDestinations,
                                                             mesh::FieldPtr<utils::Point> xiPointsSerial,
                                                             mesh::FieldPtr<utils::Point>& xiPointsParallel)
{

  if (!stk::is_true_on_any_proc(m_unionComm, m_xiPts != nullptr))
    return;

  if (middleMeshParallel)
  {
    int npts = std::max(m_xiPts->get_xi_coords(mesh::MeshEntityType::Triangle).size(),
                        m_xiPts->get_xi_coords(mesh::MeshEntityType::Quad).size());
    xiPointsParallel = mesh::create_field<utils::Point>(middleMeshParallel, mesh::impl::FieldShape(0, 0, npts), 1);
  }

  mesh::impl::FieldScatterFromRoot<utils::Point> scatterer(m_unionComm, m_rootRankOnUnionComm, 
                            entityDestinations, xiPointsSerial, xiPointsParallel);
  scatterer.scatter();
}

mesh::FieldPtr<int>
ApplicationInterfaceFakeParallelImpl::get_element_destinations(std::shared_ptr<mesh::Mesh> middleGridSerial,
                                                               mesh::FieldPtr<mesh::MeshEntityPtr> meshClassification,
                                                               mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins)
{
  if (!middleGridSerial)
    return nullptr;

  auto elementDestinationRanks =
      mesh::create_field<int>(middleGridSerial, ::stk::middle_mesh::mesh::FieldShape(0, 0, 1), 1);
  auto& destRanksField          = *elementDestinationRanks;
  auto& meshClassificationField = *meshClassification;
  auto& elementOriginsField     = *elementOrigins;
  for (auto& middleGridEl : middleGridSerial->get_elements())
    if (middleGridEl)
    {
      mesh::MeshEntityPtr inputMeshEl    = meshClassificationField(middleGridEl, 0, 0);
      destRanksField(middleGridEl, 0, 0) = elementOriginsField(inputMeshEl, 0, 0).remoteRank;
    }

  return elementDestinationRanks;
}

std::pair<int, int> get_comm_sizes_on_root(MPI_Comm comm1, MPI_Comm comm2, MPI_Comm unionComm, int rootRank)
{
  auto tag1 = get_mpi_tag_manager().get_tag(unionComm);
  auto tag2 = get_mpi_tag_manager().get_tag(unionComm);

  MPI_Request sendReq1 = MPI_REQUEST_NULL, sendReq2 = MPI_REQUEST_NULL, recvReq1 = MPI_REQUEST_NULL,
              recvReq2 = MPI_REQUEST_NULL;

  int mesh1CommSizeSend = 0, mesh2CommSizeSend = 0;
  if (comm1 != MPI_COMM_NULL && utils::impl::comm_rank(comm1) == 0)
  {
    mesh1CommSizeSend = utils::impl::comm_size(comm1);
    MPI_Isend(&mesh1CommSizeSend, 1, MPI_INT, rootRank, tag1, unionComm, &sendReq1);
  }

  if (comm2 != MPI_COMM_NULL && utils::impl::comm_rank(comm2) == 0)
  {
    mesh2CommSizeSend = utils::impl::comm_size(comm2);
    MPI_Isend(&mesh2CommSizeSend, 1, MPI_INT, rootRank, tag2, unionComm, &sendReq2);
  }

  int mesh1CommSize = 0, mesh2CommSize = 0;
  if (utils::impl::comm_rank(unionComm) == rootRank)
  {
    MPI_Irecv(&mesh1CommSize, 1, MPI_INT, MPI_ANY_SOURCE, tag1, unionComm, &recvReq1);
    MPI_Irecv(&mesh2CommSize, 1, MPI_INT, MPI_ANY_SOURCE, tag2, unionComm, &recvReq2);

    MPI_Wait(&recvReq1, MPI_STATUS_IGNORE);
    MPI_Wait(&recvReq2, MPI_STATUS_IGNORE);
  }

  MPI_Wait(&sendReq1, MPI_STATUS_IGNORE);
  MPI_Wait(&sendReq2, MPI_STATUS_IGNORE);

  return std::make_pair(mesh1CommSize, mesh2CommSize);
}

int ApplicationInterfaceFakeParallelImpl::decide_root_rank(MPI_Comm unionComm, std::shared_ptr<mesh::Mesh> mesh1,
                                                           std::shared_ptr<mesh::Mesh> /*mesh2*/)
{
  int val = -1;
  if (mesh1 && utils::impl::comm_rank(mesh1->get_comm()) == 0)
    val = utils::impl::comm_rank(unionComm);

  int rootRank;
  MPI_Allreduce(&val, &rootRank, 1, MPI_INT, MPI_MAX, unionComm);
  
  if (rootRank == -1)
    throw std::runtime_error("failed to determine root rank, was mesh1 not provided on all procs?");

  return rootRank;
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
