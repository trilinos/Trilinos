#include "application_interface_parallel.hpp"

#include "boundary_fixture.hpp"
//#include "field_base.hpp"
#include "incremental_mesh_boundary_snapper.hpp"
#include "mesh_quality_improver.hpp"
#include "middle_mesh_field_scatter.hpp"
#include "middle_mesh_point_projection.hpp"
#include "stk_middle_mesh/mesh_relational_data_scatter.hpp"
#include "stk_middle_mesh/predicates/point_classifier_normal_wrapper.hpp"
#include "stk_middle_mesh/mesh_scatter.hpp"
#include "stk_middle_mesh/mesh_scatter_spec.hpp"
#include "stk_middle_mesh/bounding_box_search.hpp"
#include "stk_middle_mesh/middle_grid_constraint_generator.hpp"
#include "stk_middle_mesh/mesh_projection_calculator.hpp"
#include "stk_middle_mesh/middle_grid_triangulator.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlocking.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "stk_middle_mesh/edge_tracer_opts.hpp"
#include "stk_middle_mesh/create_mesh_quality_improver.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

ApplicationInterfaceParallel::ApplicationInterfaceParallel(std::shared_ptr<mesh::Mesh> mesh1,
                              std::shared_ptr<mesh::Mesh> mesh2,
                              MPI_Comm unionComm, const ParallelSearchOpts& parallelSearchOpts,
                              const VolumeSnapOpts& volumeSnapOpts,
                              const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                              const MiddleGridOpts& middleGridOpts,
                              std::shared_ptr<XiCoordinates> xiPts)
  : m_mesh1(mesh1)
  , m_mesh2(mesh2)
  , m_unionComm(unionComm)
  , m_parallelSearchOpts(parallelSearchOpts)
  , m_volumeSnapOpts(volumeSnapOpts)
  , m_boundarySnapOpts(boundarySnapOpts)
  , m_middleGridOpts(middleGridOpts)
  , m_xiPts(xiPts)   
{
  check_inputs();
}

// creates the middle grid
void ApplicationInterfaceParallel::create_middle_grid()
{
  do_volume_snap();
  do_boundary_snap();

  //TODO: don't keep MeshRelationalData around the entire time
  auto scatterSpec1To2 = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_mesh1);
  auto scatterSpecReversed2To1 = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_mesh2);
  create_scatter_spec(m_mesh1, m_mesh2, scatterSpec1To2, scatterSpecReversed2To1);

  scatter_mesh_1to2(scatterSpec1To2);
  if (m_mesh2)
  {
    m_meshRelationalDataOnMesh2 = std::make_shared<nonconformal4::impl::MeshRelationalData>(m_mesh1ScatteredToMesh2, m_mesh2);
  }

  do_mesh_projections();

  scatter_mesh_2to1(scatterSpecReversed2To1);
  scatter_mesh_relational_data_2to1();

  create_middle_mesh_verts_and_edges();
  create_middle_mesh_triangles();

  project_xi_points_onto_input_meshes();
  apply_geometry_improvers();

  if (m_mesh1)
  {
    m_mesh1Classification = m_meshRelationalDataOnMesh1->meshInElementsToMesh1Elements;
  }

  scatter_meshin_1to2();
  scatter_meshin_fields_1to2();

  m_middleMeshCreated = true;
}


std::shared_ptr<mesh::Mesh> ApplicationInterfaceParallel::get_middle_grid_for_mesh1()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh1Procs, "attempting to retrieve middle grid for mesh1 even though mesh 1 is the nullptr on this process");

  return m_meshInOnMesh1Procs;
}

std::shared_ptr<mesh::Mesh> ApplicationInterfaceParallel::get_middle_grid_for_mesh2()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh2Procs, "attempting to retrieve info about middle grid for mesh2 even though mesh 1 is the nullptr on this process");
  return m_meshInOnMesh2Procs;
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceParallel::get_mesh1_classification()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh1Procs, "attempting to retrieve info about middle grid for mesh1 even though mesh 1 is the nullptr on this process");

  return m_mesh1Classification;      
}

mesh::FieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceParallel::get_mesh2_classification()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh2Procs, "attempting to retrieve info about middle grid for mesh2 even though mesh 1 is the nullptr on this process");
  return m_mesh2Classification;
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceParallel::compute_mesh1_inverse_classification()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh1Procs, "attempting to retrieve info about middle grid for mesh1 even though mesh 1 is the nullptr on this process");

  if (!m_mesh1InverseClassification)
  {
    m_mesh1InverseClassification = nonconformal4::impl::invert_classification_field(m_meshInOnMesh1Procs, m_mesh1, m_mesh1Classification);
  }

  return m_mesh1InverseClassification;
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> ApplicationInterfaceParallel::compute_mesh2_inverse_classification()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh2Procs, "attempting to retrieve info about middle grid for mesh2 even though mesh 1 is the nullptr on this process");

  if (!m_mesh2InverseClassification)
  {
    m_mesh2InverseClassification = nonconformal4::impl::invert_classification_field(m_meshInOnMesh2Procs, m_mesh2, m_mesh2Classification);
  }

  return m_mesh2InverseClassification;      
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceParallel::get_remote_info_mesh_one_to_two()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh1Procs, "attempting to retrieve info about middle grid for mesh1 even though mesh 1 is the nullptr on this process");

  return m_meshInRemoteInfo1to2;
}

mesh::FieldPtr<mesh::RemoteSharedEntity> ApplicationInterfaceParallel::get_remote_info_mesh_two_to_one()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh2Procs, "attempting to retrieve info about middle grid for mesh2 even though mesh 1 is the nullptr on this process");

  return m_meshInRemoteInfo2to1;
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceParallel::get_xi_points_on_mesh1()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh1Procs, "attempting to retrieve info about middle grid for mesh1 even though mesh 1 is the nullptr on this process");
  STK_ThrowRequireMsg(m_xiPts, "Must pass in the XiCoordinates object to retrieve information about it");

  return m_xiPtsProjectedOntoMesh1;
}

mesh::FieldPtr<utils::Point> ApplicationInterfaceParallel::get_xi_points_on_mesh2()
{
  STK_ThrowRequireMsg(m_middleMeshCreated, "Must create middle mesh before calling the getter");
  STK_ThrowRequireMsg(m_meshInOnMesh2Procs, "attempting to retrieve info about middle grid for mesh2 even though mesh 1 is the nullptr on this process");
  STK_ThrowRequireMsg(m_xiPts, "Must pass in the XiCoordinates object to retrieve information about it");

  return m_xiPtsProjectedOntoMesh2;
}



void ApplicationInterfaceParallel::check_inputs()
{
  check_meshes_exist();
  check_mesh_comm_size(m_mesh1, "mesh1");
  check_mesh_comm_size(m_mesh2, "mesh2");

  if (m_mesh1)
  {
    check_mesh_element_count(m_mesh1, "mesh1");
  }

  if (m_mesh2)
  {
    check_mesh_element_count(m_mesh2, "mesh2");
  }

#ifndef NDEBUG
  if (m_mesh1)
  {
    mesh::check_topology(m_mesh1);
  }

  if (m_mesh2)
  {
    mesh::check_topology(m_mesh2);
  }
#endif
}

void ApplicationInterfaceParallel::check_meshes_exist()
{
  std::array<int, 2> meshExistsLocally   = { static_cast<bool>(m_mesh1), static_cast<bool>(m_mesh2) };
  std::array<int, 2> meshExistsSomewhere = { false, false };

  MPI_Allreduce(meshExistsLocally.data(), meshExistsSomewhere.data(), 2, MPI_INT, MPI_LOR, m_unionComm);
  STK_ThrowRequireMsg(meshExistsSomewhere[0], 
  "Mesh1 does not exist on any procs.  Maybe this is an MPMD case where both application passed their mesh as mesh2?");

  STK_ThrowRequireMsg(meshExistsSomewhere[1], 
  "Mesh2 does not exist on any procs.  Maybe this is an MPMD case where both application passed their mesh as mesh1?");    
}

void ApplicationInterfaceParallel::check_mesh_comm_size(std::shared_ptr<mesh::Mesh> mesh, const std::string& name)
{
  int globalCount = 0;
  int localCount = mesh ? 1 : 0;
  MPI_Allreduce(&localCount, &globalCount, 1, MPI_INT, MPI_SUM, m_unionComm);
  if (mesh)
  {
    STK_ThrowRequireMsg(globalCount == utils::impl::comm_size(mesh->get_comm()), 
                        std::string("number of ") + name + " procs is not equal to the mesh1 MPI_Comm size");
  }
}

void ApplicationInterfaceParallel::check_mesh_element_count(std::shared_ptr<mesh::Mesh> mesh, const std::string& name)
{
  unsigned long long int numElsLocal = mesh::count_valid(mesh->get_elements());
  unsigned long long int numElsGlobal = 0;
  MPI_Allreduce(&numElsLocal, &numElsGlobal, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mesh->get_comm());

  STK_ThrowRequireMsg(numElsGlobal > 0, name + " has zero elements across all procs");
}

void ApplicationInterfaceParallel::do_volume_snap()
{
  if (m_volumeSnapOpts.type != VolumeSnapType::None)
  {
    std::stringstream ss;
    ss << "Volume snapping is not supported in parallel.  Instead, set:" << std::endl;
    ss << "  MiddleMeshOpts.normalProjectionOpts = NormalProjectionOpts(tol)" << std::endl;
    ss << "where tol is the tolerance you would have set for the volume snapper" << std::endl;
    throw std::runtime_error("Volume snapping is not supported in parallel.  Instead, set MiddleMeshOpts.normalProjectionOpts = NormalProjectionOpts(tol), where tol is the tolerance you would have set for the volume snapper");
  }
}

void ApplicationInterfaceParallel::do_boundary_snap()
{
  if (m_boundarySnapOpts.type == BoundarySnapAndQualityImprovementType::SnapThenQuality)
  {
    std::shared_ptr<mesh::impl::MeshQualityImprover> qualityImprover1, qualityImprover2;
    if (m_mesh1)
    {
      mesh::impl::BoundaryFixture fixture1(m_mesh1);
      qualityImprover1 = mesh::impl::make_standard_improver(m_mesh1, fixture1, m_boundarySnapOpts.snapThenQualityOpts);
    }

    if (m_mesh2)
    {
      mesh::impl::BoundaryFixture fixture2(m_mesh2);
      qualityImprover2 = mesh::impl::make_standard_improver(m_mesh2, fixture2, m_boundarySnapOpts.snapThenQualityOpts);
    }

    mesh::impl::MeshBoundarySnapper snapper;
    snapper.snap(m_mesh1, m_mesh2, m_unionComm);

    if (m_mesh1)
    {
      qualityImprover1->run();
    }

    if (m_mesh2)
    {
      qualityImprover2->run();
    }
  } else if (m_boundarySnapOpts.type == BoundarySnapAndQualityImprovementType::IncrementalBoundarySnap)
  {
    auto snapper = mesh::impl::make_incremental_boundary_snapper(m_mesh1, m_mesh2, m_unionComm,
                                                              m_boundarySnapOpts.incrementalMeshBoundarySnapOpts);
    snapper->snap();
  } else if (m_boundarySnapOpts.type != BoundarySnapAndQualityImprovementType::None)
  {
    throw std::runtime_error("unhandled boundary snap enum");
  }
}

void ApplicationInterfaceParallel::create_scatter_spec(std::shared_ptr<mesh::Mesh> sendMesh,
                          std::shared_ptr<mesh::Mesh> recvMesh,
                          std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpecSendToRecv,
                          std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpecRecvToSend)
{  
  bool isRecvMeshMesh2 = recvMesh == m_mesh2;
  mesh::FieldPtr<utils::Point> mesh2AveragedNormalVectors = nullptr;
  if (m_mesh2)
  {
    predicates::impl::AveragedNormalField avgNormalField(m_mesh2);
    mesh2AveragedNormalVectors = avgNormalField.get_field();
  }


  std::shared_ptr<stk::middle_mesh::mesh::impl::SearchMeshElementBoundingBoxBase> sendMeshSearchAdaptor, recvMeshSearchAdaptor;
  if (sendMesh)
  {
    if (isRecvMeshMesh2)
    {
      sendMeshSearchAdaptor = std::make_shared<mesh::impl::SearchMeshElementBoundingBox>(sendMesh, m_unionComm);
    } else
    {
      sendMeshSearchAdaptor = std::make_shared<mesh::impl::SearchMeshElementBoundingBoxNormal>(sendMesh, m_unionComm, 
                                              mesh2AveragedNormalVectors, m_middleGridOpts.searchOpts.normalDirectionFactor);
    }
  }

  if (recvMesh)
  {
    if (isRecvMeshMesh2)
    {
      recvMeshSearchAdaptor = std::make_shared<mesh::impl::SearchMeshElementBoundingBoxNormal>(recvMesh, m_unionComm, 
                                               mesh2AveragedNormalVectors, m_middleGridOpts.searchOpts.normalDirectionFactor);
    } else
    {
      recvMeshSearchAdaptor = std::make_shared<mesh::impl::SearchMeshElementBoundingBox>(recvMesh, m_unionComm);
    }
  }

  auto search = std::make_shared<BoundingBoxSearch>(sendMeshSearchAdaptor, recvMeshSearchAdaptor, 
                                                    "mesh1To2Search", m_unionComm, true, m_middleGridOpts.searchOpts);

  search->coarse_search();
  BoundingBoxSearch::EntityProcRelationVec& meshRecvToSendRelations = search->get_range_to_domain();

  if (search->get_unpaired_recv_entities().size() > 0)
    throw std::runtime_error("coarse search could not find a destination for some mesh1 entities");

  if (sendMesh)
  {
    int myrankOnUnionComm = utils::impl::comm_rank(m_unionComm);
    for (const BoundingBoxSearch::EntityProcRelation& searchRelation : meshRecvToSendRelations)
    {
      if (searchRelation.second.proc() == size_t(myrankOnUnionComm))
      {
        int destProc = searchRelation.first.proc();
        int meshSendElementLocalId = searchRelation.second.id();
        assert(size_t(meshSendElementLocalId) < sendMesh->get_elements().size());

        destProc = translate_comm_rank(m_unionComm, m_unionComm, destProc);
        mesh::MeshEntityPtr sendMeshEl = sendMesh->get_elements()[meshSendElementLocalId];
        scatterSpecSendToRecv->add_destination(sendMeshEl, destProc);
      }
    }
  }

  if (scatterSpecRecvToSend)
  {
    invert_coarse_search_result(sendMesh, recvMesh, meshRecvToSendRelations, scatterSpecRecvToSend);
  }
}

void ApplicationInterfaceParallel::invert_coarse_search_result(std::shared_ptr<mesh::Mesh> sendMesh,
                                  std::shared_ptr<mesh::Mesh> recvMesh,
                                  const BoundingBoxSearch::EntityProcRelationVec& meshRecvToSendRelations,
                                  std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpecRecvToSend)
{
  stk::DataExchangeUnknownPatternNonBlockingBuffer<int> exchanger(m_unionComm);
  if (sendMesh)
  {
    int myrankOnUnionComm = utils::impl::comm_rank(m_unionComm);
    for (const BoundingBoxSearch::EntityProcRelation& recvToSend : meshRecvToSendRelations)
    {
      auto recvProcAndId = recvToSend.first;
      auto sendProcAndId = recvToSend.second;          
      if (sendProcAndId.proc() == size_t(myrankOnUnionComm))
      {
        auto& sendBuf = exchanger.get_send_buf(recvProcAndId.proc());
        sendBuf.push_back(recvProcAndId.id());
      }
    }
  }

  exchanger.start_nonblocking();
  exchanger.post_nonblocking_receives();


  auto unpacker = [&](int rank, const std::vector<int>& buf)
  {
    for (int recvMeshElId : buf)
    {
      assert(recvMeshElId < int(recvMesh->get_elements().size()));
      mesh::MeshEntityPtr recvMeshEl = recvMesh->get_elements()[recvMeshElId];
      scatterSpecRecvToSend->add_destination(recvMeshEl, rank);
    }
  };

  exchanger.complete_receives(unpacker);
}

int ApplicationInterfaceParallel::translate_comm_rank(MPI_Comm inputComm, MPI_Comm outputComm, int rankOnInputComm)
{
MPI_Group inputGroup, outputGroup;
MPI_Comm_group(inputComm, &inputGroup);
MPI_Comm_group(outputComm, &outputGroup);

int outputRank;
MPI_Group_translate_ranks(inputGroup, 1, &rankOnInputComm, outputGroup, &outputRank);

return outputRank;
}    

void ApplicationInterfaceParallel::scatter_mesh_1to2(std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec)
{
  MPI_Comm scatteredMeshComm = m_mesh2 ? m_mesh2->get_comm() : MPI_COMM_NULL;
  mesh::impl::MeshScatter scatterer(scatterSpec, m_mesh1, scatteredMeshComm, true);
  m_mesh1ScatteredToMesh2   = scatterer.scatter();
  m_mesh1EntityOrigins      = scatterer.get_entity_origins();
  m_mesh1EntityDestinations = scatterer.get_entity_destinations();
}

void ApplicationInterfaceParallel::scatter_mesh_2to1(std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec)
{
  MPI_Comm scatteredMeshComm = m_mesh1 ? m_mesh1->get_comm() : MPI_COMM_NULL;
  mesh::impl::MeshScatter scatterer(scatterSpec, m_mesh2, scatteredMeshComm, true);
  m_mesh2ScatteredToMesh1   = scatterer.scatter();
  m_mesh2EntityOrigins      = scatterer.get_entity_origins();
  m_mesh2EntityDestinations = scatterer.get_entity_destinations();
}

void ApplicationInterfaceParallel::do_mesh_projections()
{
  if (m_mesh2)
  {
    m_pointClassifierOnMesh2 = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(m_mesh2, m_middleGridOpts.normalProjectionOpts.classifierTolerances);
    nonconformal4::impl::MeshProjectionCalculator meshProjection(m_mesh1ScatteredToMesh2, m_mesh2, m_meshRelationalDataOnMesh2, m_pointClassifierOnMesh2,
                                                                  m_middleGridOpts.normalProjectionOpts.edgeTracerTolerances,
                                                                  m_middleGridOpts.searchOpts);
    meshProjection.project();
  }
}

void ApplicationInterfaceParallel::scatter_mesh_relational_data_2to1()
{
  if (m_mesh2ScatteredToMesh1)
  {
    m_pointClassifierOnMesh1 = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(m_mesh2ScatteredToMesh1, m_middleGridOpts.normalProjectionOpts.classifierTolerances);
  }

  nonconformal4::impl::MeshRelationalDataScatterInput scatterInput;
  scatterInput.mesh1 = m_mesh1;
  scatterInput.mesh2 = m_mesh2;
  scatterInput.mesh1ScatteredToMesh2 = m_mesh1ScatteredToMesh2;
  scatterInput.mesh2ScatteredToMesh1 = m_mesh2ScatteredToMesh1;
  scatterInput.mesh1EntityOrigins = m_mesh1EntityOrigins;
  scatterInput.mesh1EntityDestinations = m_mesh1EntityDestinations;
  scatterInput.mesh2EntityOrigins = m_mesh2EntityOrigins;
  scatterInput.mesh2EntityDestinations = m_mesh2EntityDestinations;

  nonconformal4::impl::MeshRelationalDataScatter scatterer(scatterInput, m_meshRelationalDataOnMesh2,
                                                            m_pointClassifierOnMesh2,
                                                            m_pointClassifierOnMesh1,
                                                            MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
  m_meshRelationalDataOnMesh1 = scatterer.scatter();
  m_meshInOnMesh1Procs        = scatterer.get_middle_mesh();
}

void ApplicationInterfaceParallel::create_middle_mesh_verts_and_edges()
{
  if (m_mesh1)
  {
    nonconformal4::impl::MiddleGridConstraintGenerator generator(
      m_mesh1, m_mesh2ScatteredToMesh1, m_meshInOnMesh1Procs, m_meshRelationalDataOnMesh1, m_pointClassifierOnMesh1);

    generator.generate();
  }
}

void ApplicationInterfaceParallel::create_middle_mesh_triangles()
{
  if (m_mesh1)
  {
    nonconformal4::impl::MiddleGridTriangulator triangulator(m_mesh1, m_mesh2ScatteredToMesh1,
                            m_meshInOnMesh1Procs, m_meshRelationalDataOnMesh1, m_pointClassifierOnMesh1);

    triangulator.triangulate();
  }
}

void ApplicationInterfaceParallel::project_xi_points_onto_input_meshes()
{
  if (m_mesh1 && m_xiPts)
  {
    nonconformal4::impl::MiddleMeshPointProjection projector(m_mesh1, m_mesh2ScatteredToMesh1,
                    m_meshInOnMesh1Procs, m_meshRelationalDataOnMesh1, m_pointClassifierOnMesh1);
    m_xiPtsProjectedOntoMesh1                 = projector.projection_onto_mesh1(m_xiPts);
    m_xiPtsProjectedOntoMesh2ScatteredToMesh1 = projector.projection_onto_mesh2(m_xiPts);
  }
}

void ApplicationInterfaceParallel::apply_geometry_improvers()
{
  if (m_middleGridOpts.normalProjectionOpts.geometryImprovers.size() > 0)
  {
    throw std::runtime_error("GeometryImprovers are not available in parallel (yet)");
  }
}


void ApplicationInterfaceParallel::scatter_meshin_1to2()
{
  auto scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(m_unionComm, m_meshInOnMesh1Procs);

  if (m_meshInOnMesh1Procs)
  {
    auto& meshInToMesh2Els = *(m_meshRelationalDataOnMesh1->meshInElementsToMesh2Elements);
    auto& mesh2EntityOrigins = *m_mesh2EntityOrigins;

    for (auto& elIn : m_meshInOnMesh1Procs->get_elements())
      if (elIn)
      {
        mesh::MeshEntityPtr el2ScatteredToMesh1 = meshInToMesh2Els(elIn, 0, 0);
        mesh::RemoteSharedEntity origin = mesh2EntityOrigins(el2ScatteredToMesh1, 0, 0);
        scatterSpec->add_destination(elIn, origin.remoteRank);
      }
  }

  MPI_Comm scatteredMeshComm = m_mesh2 ? m_mesh2->get_comm() : MPI_COMM_NULL;
  mesh::impl::MeshScatter scatterer(scatterSpec, m_meshInOnMesh1Procs, scatteredMeshComm, true);

  m_meshInOnMesh2Procs = scatterer.scatter();
  create_meshin_remote_info(scatterer);
}

void ApplicationInterfaceParallel::create_meshin_remote_info(const mesh::impl::MeshScatter& scatterer)
{
  if (m_meshInOnMesh1Procs)
  {
    m_meshInRemoteInfo1to2 = extract_element_remote_info(m_meshInOnMesh1Procs,
                                                          scatterer.get_entity_destinations());
  }

  if (m_meshInOnMesh2Procs)
  {
    m_meshInRemoteInfo2to1 = extract_element_remote_info(m_meshInOnMesh2Procs, 
                                                          scatterer.get_entity_origins());
  }
}

mesh::FieldPtr<mesh::RemoteSharedEntity>
ApplicationInterfaceParallel::extract_element_remote_info(std::shared_ptr<mesh::Mesh> meshIn,
                                                          mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityRemoteInfoPtr)
{
  auto elementRemoteInfoPtr = mesh::create_field<mesh::RemoteSharedEntity>(meshIn, mesh::FieldShape(0, 0, 1), 1);

  auto& entityRemoteInfo = *entityRemoteInfoPtr;
  auto& elementRemoteInfo = *elementRemoteInfoPtr;
  for (auto& el : meshIn->get_elements())
  {
    if (el)
    {
      assert(entityRemoteInfo.get_num_comp(el, 0) == 1);
      elementRemoteInfo(el, 0, 0) = entityRemoteInfo(el, 0, 0);
    }
  }

  return elementRemoteInfoPtr;
}

void ApplicationInterfaceParallel::scatter_meshin_fields_1to2()
{
  mesh::FieldPtr<mesh::MeshEntityPtr> meshInClassificationOnMesh2ScatteredToMesh1;
  if (m_mesh1)
  {
    meshInClassificationOnMesh2ScatteredToMesh1 = m_meshRelationalDataOnMesh1->meshInElementsToMesh2Elements;
  }
  nonconformal4::impl::MiddleMeshFieldScatter scatterer(
                                                m_meshInOnMesh1Procs, m_meshInOnMesh2Procs, m_mesh2,
                                                m_xiPts, meshInClassificationOnMesh2ScatteredToMesh1,
                                                m_xiPtsProjectedOntoMesh2ScatteredToMesh1,
                                                m_mesh2EntityOrigins,
                                                m_meshInRemoteInfo1to2,
                                                m_meshInRemoteInfo2to1,
                                                m_unionComm);
  scatterer.scatter();
  m_mesh2Classification = scatterer.get_meshin_classification_on_mesh2();
  m_xiPtsProjectedOntoMesh2 = scatterer.get_xi_pts_projected_onto_mesh2();

}

}
}
}
