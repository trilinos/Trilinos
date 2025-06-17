// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "SimpleTransfer.hpp"

#include "stk_transfer_util/spmd/GeometricTransfer.hpp"
#include "stk_transfer_util/spmd/GeometricTransferDispatch.hpp"
#include "stk_transfer_util/spmd/GeometricTransferOptions.hpp"
#include "stk_transfer_util/spmd/GeometricTransferUtils.hpp"
#include <stk_transfer_util/spmd/ElementRecvMesh.hpp>
#include <stk_transfer_util/spmd/ElementSendMesh.hpp>
#include <stk_transfer_util/spmd/NodeRecvMesh.hpp>
#include <stk_transfer_util/spmd/NodeSendMesh.hpp>
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include <stk_mesh/base/Entity.hpp>         // for Entity
#include <stk_mesh/base/FieldBase.hpp>      // for FieldBase
#include "stk_mesh/base/MetaData.hpp"       // for MetaData
#include "stk_mesh/base/Selector.hpp"       // for Selector
#include "stk_mesh/base/Bucket.hpp"         // for Bucket
#include <stk_mesh/base/BulkData.hpp>       // for BulkData
#include "stk_mesh/base/CompositeRank.hpp"                 // for CompositeRank
#include "stk_topology/topology.hpp"        // for topology, operator<<, top...
#include "stk_util/environment/RuntimeWarning.hpp"         // for RuntimeWar...
#include "stk_util/util/ReportHandler.hpp"  // for eval_test_condition, STK_...
#include "stk_util/util/SortAndUnique.hpp"

#include <algorithm>                        // for max, sort, fill
#include <cstddef>                          // for size_t
#include <memory>                           // for shared_ptr, __shared_ptr_...
#include <type_traits>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

namespace impl {

void add_field(std::vector<stk::transfer::FieldSpec> & fieldSpecs, stk::mesh::FieldBase & field)
{
  stk::transfer::FieldSpec fieldSpec(field.name(), field.state());
  fieldSpecs.push_back(fieldSpec);

  // all fields need same entity rank
  const stk::mesh::MetaData& meta = field.mesh_meta_data();

  stk::mesh::FieldBase* firstField = meta.get_field(field.entity_rank(), fieldSpecs.front().name);
  STK_ThrowRequire(nullptr != firstField);
}

void add_part_name(std::vector<std::string>& partNames, const std::string& partName)
{
  stk::util::insert_keep_sorted_and_unique(partName, partNames);
}

void add_part_name(std::vector<std::string>& partNames, const std::vector<std::string>& partNames_)
{
  for(const auto& partName_ : partNames_) {
    stk::util::insert_keep_sorted_and_unique(partName_, partNames);
  }
}

}

stk::mesh::EntityRank SimpleTransfer::recv_type_to_entity_rank(stk::transfer::spmd::RecvMeshType type) const
{
  if(type == stk::transfer::spmd::RecvMeshType::NODE)                { return stk::topology::NODE_RANK; }
  if(type == stk::transfer::spmd::RecvMeshType::EDGE_CENTROID)       { return stk::topology::EDGE_RANK; }
  if(type == stk::transfer::spmd::RecvMeshType::EDGE_GAUSS_POINT)    { return stk::topology::EDGE_RANK; }
  if(type == stk::transfer::spmd::RecvMeshType::FACE_CENTROID)       { return stk::topology::FACE_RANK; }
  if(type == stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT)    { return stk::topology::FACE_RANK; }
  if(type == stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID)    { return stk::topology::ELEM_RANK; }
  if(type == stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT) { return stk::topology::ELEM_RANK; }

  return stk::topology::INVALID_RANK;
}

stk::transfer::spmd::SendMeshType SimpleTransfer::entity_rank_to_send_type(stk::mesh::EntityRank rank) const
{
  if(rank == stk::topology::NODE_RANK) { return stk::transfer::spmd::SendMeshType::NODE; }
  if(rank == stk::topology::EDGE_RANK) { return stk::transfer::spmd::SendMeshType::EDGE; }
  if(rank == stk::topology::FACE_RANK) { return stk::transfer::spmd::SendMeshType::FACE; }
  if(rank == stk::topology::ELEM_RANK) { return stk::transfer::spmd::SendMeshType::ELEMENT; }

  return stk::transfer::spmd::SendMeshType::INVALID;
}

void SimpleTransfer::add_send_field(stk::mesh::FieldBase & field)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  impl::add_field(m_sendFieldSpecs, field);
}

void SimpleTransfer::add_send_field(const stk::transfer::FieldSpec & fieldSpec)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  stk::util::insert_keep_sorted_and_unique(fieldSpec, m_sendFieldSpecs);
}

void SimpleTransfer::add_send_fields(const std::vector<stk::mesh::FieldBase*>& fields)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");

  for(stk::mesh::FieldBase* field : fields) {
    impl::add_field(m_sendFieldSpecs, *field);
  }
}

void SimpleTransfer::add_send_fields(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");

  for(const stk::transfer::FieldSpec& fieldSpec: fieldSpecs) {
    add_send_field(fieldSpec);
  }
}

void SimpleTransfer::add_recv_field(stk::mesh::FieldBase & field)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  impl::add_field(m_recvFieldSpecs, field);
}

void SimpleTransfer::add_recv_field(const stk::transfer::FieldSpec & fieldSpec)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  stk::util::insert_keep_sorted_and_unique(fieldSpec, m_recvFieldSpecs);
}

void SimpleTransfer::add_recv_fields(const std::vector<stk::mesh::FieldBase*>& fields)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");

  for(stk::mesh::FieldBase* field : fields) {
    impl::add_field(m_recvFieldSpecs, *field);
  }
}

void SimpleTransfer::add_recv_fields(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");

  for(const stk::transfer::FieldSpec& fieldSpec: fieldSpecs) {
    add_recv_field(fieldSpec);
  }
}

void SimpleTransfer::add_send_part_name(const std::string& partName)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  impl::add_part_name(m_sendPartNames, partName);
}

void SimpleTransfer::add_send_part_names(const std::vector<std::string>& partNames)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  impl::add_part_name(m_sendPartNames, partNames);
}

void SimpleTransfer::add_recv_part_name(const std::string& partName)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  impl::add_part_name(m_recvPartNames, partName);
}

void SimpleTransfer::add_recv_part_names(const std::vector<std::string>& partNames)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  impl::add_part_name(m_recvPartNames, partNames);
}

void SimpleTransfer::set_send_mesh_transfer_options(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::EntityRank sendRank,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = m_transferOptions->get_send_options();
  stk::transfer::spmd::SendMeshType sendType = entity_rank_to_send_type(sendRank);

  sendOptions.set_mesh_type(sendType);

  stk::mesh::MetaData& sendMeta = sendBulk.mesh_meta_data();
  sendOptions.set_coordinate_field(sendMeta.coordinate_field());

  sendOptions.set_extrapolate_option(extrapolateOption);
  sendOptions.set_element_mesh_entity_rank(sendRank);

  // parts and fields
  sendOptions.add_part(m_sendPartNames);
  sendOptions.add_field_spec(m_sendFieldSpecs);
}

void SimpleTransfer::set_send_mesh_transfer_options(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::EntityRank sendRank,
    std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = m_transferOptions->get_send_options();
  stk::transfer::spmd::SendMeshType sendType = entity_rank_to_send_type(sendRank);

  sendOptions.set_mesh_type(sendType);

  stk::mesh::MetaData& sendMeta = sendBulk.mesh_meta_data();
  sendOptions.set_coordinate_field(sendMeta.coordinate_field());

  sendOptions.set_extrapolate_option(extrapolateOption);
  sendOptions.set_element_mesh_entity_rank(sendRank);

  // Use defaults to create other master element components
  sendOptions.set_master_element_provider(masterElemProvider);

  // parts and fields
  sendOptions.add_part(m_sendPartNames);
  sendOptions.add_field_spec(m_sendFieldSpecs);
}

void SimpleTransfer::set_recv_mesh_transfer_options(stk::mesh::BulkData& recvBulk,
                                                    stk::transfer::spmd::RecvMeshType recvType)
{
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = m_transferOptions->get_recv_options();
  stk::mesh::EntityRank rank = recv_type_to_entity_rank(recvType);

  recvOptions.set_mesh_type(recvType);

  stk::mesh::MetaData& recvMeta = recvBulk.mesh_meta_data();
  recvOptions.set_coordinate_field(recvMeta.coordinate_field());

  recvOptions.set_element_mesh_entity_rank(rank);

  // parts and fields
  recvOptions.add_part(m_recvPartNames);
  recvOptions.add_field_spec(m_recvFieldSpecs);
}

void SimpleTransfer::set_recv_mesh_transfer_options(
    stk::mesh::BulkData& recvBulk,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
{
  set_recv_mesh_transfer_options(recvBulk, recvType);

  // Use defaults to create other master element components
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = m_transferOptions->get_recv_options();
  recvOptions.set_master_element_provider(masterElemProvider);
}

void SimpleTransfer::setup_master_element_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> sendMasterElemProvider,
    std::shared_ptr<stk::search::MasterElementProviderInterface> recvMasterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  STK_ThrowRequireMsg(sendMasterElemProvider, "NULL send master element provider for transfer named: " << m_transferName);
  STK_ThrowRequireMsg(recvMasterElemProvider, "NULL recv master element provider for transfer named: " << m_transferName);

  // Use default values for the rest
  m_transferOptions = std::make_shared<stk::transfer::spmd::GeometricTransferOptions>(m_transferName, sendBulk, recvBulk);
  m_transferOptions->set_interpolation_type(stk::transfer::spmd::InterpolationType::MASTER_ELEMENT);

  if(m_hasParallelMachine) {
    m_transferOptions->set_parallel_machine(m_parallelMachine);
  }

  set_send_mesh_transfer_options(sendBulk, sendRank, sendMasterElemProvider, extrapolateOption);
  set_recv_mesh_transfer_options(recvBulk, recvType, recvMasterElemProvider);

  std::tie(m_transfer, m_dispatch) = stk::transfer::spmd::create_transfer(*m_transferOptions);

  m_committed = true;
}

void SimpleTransfer::setup_master_element_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  setup_master_element_transfer(sendBulk, recvBulk, sendRank, recvType,
                                masterElemProvider, masterElemProvider, extrapolateOption);
}

void SimpleTransfer::setup_patch_recovery_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> sendMasterElemProvider,
    std::shared_ptr<stk::search::MasterElementProviderInterface> recvMasterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  STK_ThrowRequireMsg(sendMasterElemProvider, "NULL send master element provider for transfer named: " << m_transferName);
  STK_ThrowRequireMsg(recvMasterElemProvider, "NULL recv master element provider for transfer named: " << m_transferName);

  // Use default values for the rest
  m_transferOptions = std::make_shared<stk::transfer::spmd::GeometricTransferOptions>(m_transferName, sendBulk, recvBulk);
  m_transferOptions->set_interpolation_type(stk::transfer::spmd::InterpolationType::PATCH_RECOVERY);
  m_transferOptions->get_send_options().set_patch_recovery_type(m_defaultPatchRecovery);

  if(m_hasParallelMachine) {
    m_transferOptions->set_parallel_machine(m_parallelMachine);
  }

  set_send_mesh_transfer_options(sendBulk, sendRank, sendMasterElemProvider, extrapolateOption);
  set_recv_mesh_transfer_options(recvBulk, recvType, recvMasterElemProvider);

  std::tie(m_transfer, m_dispatch) = stk::transfer::spmd::create_transfer(*m_transferOptions);

  m_committed = true;
}

void SimpleTransfer::setup_patch_recovery_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  setup_patch_recovery_transfer(sendBulk, recvBulk, sendRank, recvType,
                                masterElemProvider, masterElemProvider, extrapolateOption);
}

void SimpleTransfer::setup_copy_nearest_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> sendMasterElemProvider,
    std::shared_ptr<stk::search::MasterElementProviderInterface> recvMasterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  STK_ThrowRequireMsg(sendMasterElemProvider, "NULL send master element provider for transfer named: " << m_transferName);
  STK_ThrowRequireMsg(recvMasterElemProvider, "NULL recv master element provider for transfer named: " << m_transferName);

  // Use default values for the rest
  m_transferOptions = std::make_shared<stk::transfer::spmd::GeometricTransferOptions>(m_transferName, sendBulk, recvBulk);
  m_transferOptions->set_interpolation_type(stk::transfer::spmd::InterpolationType::COPY);

  if(m_hasParallelMachine) {
    m_transferOptions->set_parallel_machine(m_parallelMachine);
  }

  set_send_mesh_transfer_options(sendBulk, sendRank, sendMasterElemProvider, extrapolateOption);
  set_recv_mesh_transfer_options(recvBulk, recvType, recvMasterElemProvider);

  std::tie(m_transfer, m_dispatch) = stk::transfer::spmd::create_transfer(*m_transferOptions);

  m_committed = true;
}

void SimpleTransfer::setup_copy_nearest_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  setup_copy_nearest_transfer(sendBulk, recvBulk, sendRank, recvType,
                              masterElemProvider, masterElemProvider, extrapolateOption);
}

void SimpleTransfer::setup_sum_nearest_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> sendMasterElemProvider,
    std::shared_ptr<stk::search::MasterElementProviderInterface> recvMasterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
  STK_ThrowRequireMsg(sendMasterElemProvider, "NULL send master element provider for transfer named: " << m_transferName);
  STK_ThrowRequireMsg(recvMasterElemProvider, "NULL recv master element provider for transfer named: " << m_transferName);

  // Use default values for the rest
  m_transferOptions = std::make_shared<stk::transfer::spmd::GeometricTransferOptions>(m_transferName, sendBulk, recvBulk);
  m_transferOptions->set_interpolation_type(stk::transfer::spmd::InterpolationType::SUM);

  if(m_hasParallelMachine) {
    m_transferOptions->set_parallel_machine(m_parallelMachine);
  }

  set_send_mesh_transfer_options(sendBulk, sendRank, sendMasterElemProvider, extrapolateOption);
  set_recv_mesh_transfer_options(recvBulk, recvType, recvMasterElemProvider);

  std::tie(m_transfer, m_dispatch) = stk::transfer::spmd::create_transfer(*m_transferOptions);

  m_committed = true;
}

void SimpleTransfer::setup_sum_nearest_transfer(
    stk::mesh::BulkData& sendBulk,
    stk::mesh::BulkData& recvBulk,
    stk::mesh::EntityRank sendRank,
    stk::transfer::spmd::RecvMeshType recvType,
    std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
    stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  setup_sum_nearest_transfer(sendBulk, recvBulk, sendRank, recvType,
                             masterElemProvider, masterElemProvider, extrapolateOption);
}

void SimpleTransfer::initialize_meshes()
{
  STK_ThrowRequireMsg(m_dispatch, "Transfer dispatch not created");
  STK_ThrowRequireMsg(!m_dispatch->is_initialized(), "Transfer already initialized");

  m_dispatch->initialize_meshes();
}

void SimpleTransfer::set_transfer_options(bool overrideParametricSearchWithGeometricClosestNodeToCentroid,
                                          bool useCentroidForGeometricSearchProximity, bool printSearchWarnings)
{
  STK_ThrowRequireMsg(m_dispatch, "Transfer dispatch not created");

  m_dispatch->closest_bounding_box_using_nearest_node(overrideParametricSearchWithGeometricClosestNodeToCentroid);
  m_dispatch->use_centroid_for_geometric_proximity(useCentroidForGeometricSearchProximity);
  m_dispatch->print_search_warnings(printSearchWarnings);
}

void SimpleTransfer::execute_search()
{
  STK_ThrowRequireMsg(m_transfer && m_dispatch, "Transfer not created");
  STK_ThrowRequireMsg(!m_dispatch->is_initialized(), "Transfer already initialized");

  m_transfer->coarse_search();
  m_dispatch->ghost_from_elements();
  m_transfer->local_search();
}

void SimpleTransfer::initialize(bool overrideParametricSearchWithGeometricClosestNodeToCentroid,
                                    bool useCentroidForGeometricProximity, bool printSearchWarnings)
{
  STK_ThrowRequireMsg(m_dispatch, "Transfer dispatch not created");

  if(m_dispatch->is_initialized()) return;

  initialize_meshes();

  set_transfer_options(overrideParametricSearchWithGeometricClosestNodeToCentroid,
                       useCentroidForGeometricProximity, printSearchWarnings);

  execute_search();

  m_transferOptions->check_field_validity();

  m_dispatch->set_initialized(true);
}

void SimpleTransfer::apply()
{
  STK_ThrowRequireMsg(m_transfer, "Transfer not created");
  STK_ThrowRequireMsg(m_dispatch, "Transfer dispatch not created");
  STK_ThrowRequireMsg(m_dispatch->is_initialized(), "Transfer not initialized");

  m_transfer->apply();
}

} // namespace spmd
} // namespace transfer
} // namespace stk


