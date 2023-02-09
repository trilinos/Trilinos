#include "application_interface_fake_parallel_impl.h"

#include "gathered_mesh_coordinate_update.h"
#include "incremental_mesh_boundary_snapper.h"
#include "mesh_classification_scatter.h"
#include "mesh_gather_to_root.h"
#include "mesh_scatter_from_root.h"
#include "mesh_snapper.h"
#include "nonconformal4.h"

#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace impl {

void ApplicationInterfaceFakeParallelImpl::check_union_comm_size()
{
  int rootRank   = 0;
  MPI_Comm comm1 = m_mesh1Parallel ? m_mesh1Parallel->get_comm() : MPI_COMM_NULL;
  MPI_Comm comm2 = m_mesh2Parallel ? m_mesh2Parallel->get_comm() : MPI_COMM_NULL;
  auto commSizes = impl::get_comm_sizes_on_root(comm1, comm2, m_unionComm, rootRank);

  if (utils::impl::comm_rank(m_unionComm) == rootRank &&
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
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelImpl::get_middle_grid_for_mesh1()
{
  if (m_mesh1Parallel)
    return m_middleGridParallel1;
  else
    throw std::runtime_error(
        "attempting to retrieve middle grid for mesh1 even though mesh 1 is the nullptr on this process");
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceFakeParallelImpl::get_middle_grid_for_mesh2()
{
  if (m_mesh2Parallel)
    return m_middleGridParallel2;
  else
    throw std::runtime_error(
        "attempting to retrieve middle grid for mesh 2 even though mesh 2 is the nullptr on this process");
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelImpl::get_mesh1_classification()
{
  if (m_mesh1Parallel)
    return m_mesh1ClassificationParallel;
  else
    throw std::runtime_error(
        "attempting to retrieve mesh 1 classification even though mesh 1 is the nullptr on this process");
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceFakeParallelImpl::get_mesh2_classification()
{
  if (m_mesh2Parallel)
    return m_mesh2ClassificationParallel;
  else
    throw std::runtime_error(
        "attempting to retrieve mesh 2 classification even though mesh 2 is the nullptr on this process");
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
ApplicationInterfaceFakeParallelImpl::compute_mesh1_inverse_classification()
{
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
  if (m_mesh2Parallel)
    return nonconformal4::impl::invert_classification_field(m_middleGridParallel2, m_mesh2Parallel,
                                                            m_mesh2ClassificationParallel);
  else
    throw std::runtime_error(
        "attempting to compute mesh 2 inverse classification even though mesh 2 is the nullptr on this process");
}

void ApplicationInterfaceFakeParallelImpl::gather_meshes_to_root()
{
  const int root = 0;

  mesh::impl::MeshGatherToRoot gather1(m_unionComm, root, m_mesh1Parallel);
  m_mesh1Serial         = gather1.gather();
  m_mesh1ElementOrigins = gather1.get_element_origins();

  mesh::impl::MeshGatherToRoot gather2(m_unionComm, root, m_mesh2Parallel);
  m_mesh2Serial         = gather2.gather();
  m_mesh2ElementOrigins = gather2.get_element_origins();

  if (utils::impl::comm_rank(m_unionComm) != root && m_mesh1Serial != nullptr)
    throw std::runtime_error("mesh1 serial should be null");

  if (utils::impl::comm_rank(m_unionComm) != root && m_mesh1Serial != nullptr)
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
    snapper.snap(m_mesh1Serial, m_mesh2Serial);

    qualityImprover1->run();
    qualityImprover2->run();
  } else if (m_boundarySnapOpts.type == BoundarySnapAndQualityImprovementType::IncrementalBoundarySnap)
  {
    auto snapper = mesh::impl::make_incremental_boundary_snapper(m_mesh1Serial, m_mesh2Serial,
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

  } else
  {
    throw std::runtime_error("unhandled middle grid enum");
  }
}

void ApplicationInterfaceFakeParallelImpl::scatter_meshes()
{
  scatter_mesh(m_mesh1ClassificationSerial, m_mesh1ElementOrigins, m_mesh1Parallel, m_middleGridParallel1,
               m_mesh1ClassificationParallel);
  scatter_mesh(m_mesh2ClassificationSerial, m_mesh2ElementOrigins, m_mesh2Parallel, m_middleGridParallel2,
               m_mesh2ClassificationParallel);

  mesh::impl::GatheredMeshCoordinateUpdate updator1(m_unionComm, m_mesh1Serial, m_mesh1ElementOrigins, m_mesh1Parallel);
  updator1.update();

  mesh::impl::GatheredMeshCoordinateUpdate updator2(m_unionComm, m_mesh2Serial, m_mesh2ElementOrigins, m_mesh2Parallel);
  updator2.update();
}

void ApplicationInterfaceFakeParallelImpl::scatter_mesh(mesh::FieldPtr<mesh::MeshEntityPtr> meshClassificationSerial,
                                                        mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins,
                                                        std::shared_ptr<mesh::Mesh> inputMeshParallel,
                                                        std::shared_ptr<mesh::Mesh>& middleGridParallel,
                                                        mesh::FieldPtr<mesh::MeshEntityPtr>& meshClassificationParallel)
{
  auto elementDestinations = get_element_destinations(m_middleGridSerial, meshClassificationSerial, elementOrigins);
  MPI_Comm meshComm        = inputMeshParallel ? inputMeshParallel->get_comm() : MPI_COMM_NULL;
  mesh::impl::MeshScatterFromRoot scatter1(m_unionComm, m_middleGridSerial, meshComm, elementDestinations);
  middleGridParallel = scatter1.scatter();

  mesh::impl::MeshClassificationScatter classificationScatter(
      m_unionComm, m_middleGridSerial, middleGridParallel, inputMeshParallel, elementOrigins, meshClassificationSerial);
  meshClassificationParallel = classificationScatter.scatter();
}

mesh::FieldPtr<int>
ApplicationInterfaceFakeParallelImpl::get_element_destinations(std::shared_ptr<mesh::Mesh> middleGridSerial,
                                                               mesh::FieldPtr<mesh::MeshEntityPtr> meshClassification,
                                                               mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins)
{
  if (!middleGridSerial)
    return nullptr;

  auto elementDestinationRanks =
      mesh::create_field<int>(middleGridSerial, ::stk::middle_mesh::mesh::impl::FieldShape(0, 0, 1), 1);
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

  if (comm1 != MPI_COMM_NULL && utils::impl::comm_rank(comm1) == 0)
  {
    int mesh1CommSizeSend = utils::impl::comm_size(comm1);
    MPI_Isend(&mesh1CommSizeSend, 1, MPI_INT, rootRank, tag1, unionComm, &sendReq1);
  }

  if (comm2 != MPI_COMM_NULL && utils::impl::comm_rank(comm2) == 0)
  {
    int mesh2CommSizeSend = utils::impl::comm_size(comm2);
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

} // namespace impl
} // namespace middle_mesh
} // namespace stk
