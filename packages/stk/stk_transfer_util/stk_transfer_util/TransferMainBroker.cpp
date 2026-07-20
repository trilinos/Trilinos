/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "mpi.h"

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/CompositeRank.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_search_util/MasterElementProvider.hpp"
#include "stk_search_util/MasterElementProviderIntrepid2.hpp"
#include "stk_search_util/CoordTransform.hpp"
#include "stk_transfer_util/TransferMainBroker.hpp"
#include "stk_transfer_util/TransferMainSettings.hpp"
#include "stk_transfer_util/spmd/SimpleTransfer.hpp"
#include "stk_transfer_util/LogMessage.hpp"

namespace stk {
namespace transfer_util {

std::vector<std::string> get_default_parts_for_transfer(const stk::mesh::PartVector & parts,
                                                        const stk::mesh::EntityRank defaultRank)
{
  std::vector<std::string> partNamesForTransfer;

  for (auto part : parts) {
    if (part->primary_entity_rank() == defaultRank) {
      partNamesForTransfer.push_back(part->name());
    }
  }

  return partNamesForTransfer;
}

TransferMainBroker::TransferMainBroker(MPI_Comm comm,
                                       std::shared_ptr<stk::mesh::BulkData> sendBulk,
                                       std::shared_ptr<stk::mesh::BulkData> recvBulk,
                                       const TransferMainSettings& settings,
                                       std::shared_ptr<stk::search::MasterElementProviderInterface> handlerMEP)
  : m_comm(comm),
    m_sendBulk(sendBulk),
    m_recvBulk(recvBulk),
    m_sendMeta(sendBulk->mesh_meta_data()),
    m_recvMeta(recvBulk->mesh_meta_data()),
    m_settings(settings),
    m_masterElementProvider(handlerMEP)
{

  if (m_masterElementProvider == nullptr) {
    const std::string& masterElemName = m_settings.get_master_elements_name();
    STK_ThrowRequireMsg(masterElemName == "INTREPID2" || masterElemName == "",
                       "Invalid MasterElements name '" << masterElemName << "', Currently 'INTREPID2' is the only choice.");
    m_masterElementProvider = std::make_shared<stk::search::MasterElementProviderIntrepid2>();
  }

  m_sendPartNamesForTransfer = m_settings.get_transfer_send_parts();
  m_recvPartNamesForTransfer = m_settings.get_transfer_recv_parts();

  if (m_sendPartNamesForTransfer.size() == 0) {
    //currently using element-rank send parts as default due to Intrepid2 limitation
    auto sendDefaultRank = stk::topology::ELEMENT_RANK;
    m_sendPartNamesForTransfer = get_default_parts_for_transfer(m_sendMeta.get_mesh_parts(), sendDefaultRank);
  }

  if (m_recvPartNamesForTransfer.size() == 0) {
    auto recvDefaultRank = m_settings.get_default_recv_part_rank();
    m_recvPartNamesForTransfer = get_default_parts_for_transfer(m_recvMeta.get_mesh_parts(), recvDefaultRank);
  }

  STK_ThrowRequireMsg(m_sendPartNamesForTransfer.size() > 0, "Cannot complete transfer. There are no send parts!");
  STK_ThrowRequireMsg(m_recvPartNamesForTransfer.size() > 0, "Cannot complete transfer. There are no recv parts!");
}

stk::mesh::FieldBase* TransferMainBroker::add_created_recvField(const stk::mesh::FieldBase* sendField, const std::string& fieldName)
{
  m_recvMeta.enable_late_fields();

  auto createdRecvField = sendField->clone(m_recvMeta.get_field_repository(), fieldName, m_settings.get_recv_field_rank());
  m_recvMeta.set_mesh_on_field(m_recvBulk.get(), *createdRecvField);

  stk::mesh::PartVector recvParts;
  for (auto partName : m_recvPartNamesForTransfer) {
    recvParts.push_back(m_recvMeta.get_part(partName));
  }

  stk::mesh::Selector recvPartSelector = stk::mesh::selectUnion(recvParts);

  int numCopies = 1;
  int numComponents = sendField->max_size();
  const std::byte* initValsPtr = createdRecvField->get_initial_value_bytes().extent(0) > 0 ? createdRecvField->get_initial_value_bytes().data() : nullptr;

  if (m_settings.get_recv_type_string() == "ELEMENT_GAUSS_POINT") {
    for (auto* part : recvParts) {
       if (part->primary_entity_rank() == stk::topology::ELEM_RANK) {
         stk::search::SearchTopology searchTopology(part->topology());
         numCopies = m_masterElementProvider->num_integration_points(searchTopology);
         numComponents = sendField->max_size()/numCopies;
         stk::mesh::put_field_on_mesh(*createdRecvField, *part, numComponents, numCopies, initValsPtr);
       }
    }
  }
  else if (m_settings.get_recv_type_string().ends_with("GAUSS_POINT")) {
    STK_ThrowErrorMsg("Creation of receive field not handled for "
                          + m_settings.get_recv_type_string() + " receive type");
  }
  else {

    for(const auto& r : sendField->restrictions()) {
      if (r.dimension() != numComponents || r.num_scalars_per_entity() != numComponents*numCopies) {
      STK_ThrowErrorMsg("Only constant-size fields are supported at this time. Cannot create a receive field from "
                              + sendField->name());
      }
    }

    stk::mesh::put_field_on_mesh(*createdRecvField, recvPartSelector, numComponents, numCopies, initValsPtr);
  }

  auto fieldOutputType = stk::io::get_field_output_type(*sendField);
  stk::io::set_field_output_type(*createdRecvField, fieldOutputType);
  stk::io::set_field_role(*createdRecvField, Ioss::Field::TRANSIENT);

  m_recvFieldsForTransfer.push_back(createdRecvField);

  return createdRecvField;

}

void check_send_field_for_interpolation(const stk::mesh::FieldBase* sendField, const std::string& sendFieldName)
{
  STK_ThrowRequireMsg(sendField->entity_rank() == stk::topology::NODE_RANK,
                                         "Specified send field " + sendFieldName + " must have NODE rank for interpolation transfer");
}

void check_send_field_for_copy(const stk::mesh::FieldBase* sendField, const std::string& sendFieldName)
{
  STK_ThrowRequireMsg(sendField->entity_rank() == stk::topology::ELEMENT_RANK ||
                                                  sendField->entity_rank() == stk::topology::NODE_RANK,
                                         "Specified send field " + sendFieldName + " must have NODE or ELEMENT rank for copy transfer");
}

void check_send_field_for_patch_recovery(const stk::mesh::FieldBase* sendField, const std::string& sendFieldName)
{
  STK_ThrowRequireMsg(sendField->entity_rank() == stk::topology::ELEMENT_RANK,
                                         "Specified send field " + sendFieldName + " must have ELEMENT rank for patch recovery transfer");
}

void TransferMainBroker::check_send_field(const stk::mesh::FieldBase* sendField, const std::string& fieldName)
{
  STK_ThrowRequireMsg(sendField != nullptr, "Specified field " + fieldName + " was not found on the send mesh!");

  if(stk::io::get_field_role(*sendField) != nullptr) {
    STK_ThrowRequireMsg(*stk::io::get_field_role(*sendField) == Ioss::Field::TRANSIENT,
                                         "Specified send field " + fieldName + " must be TRANSIENT");
  }

  if (m_settings.get_transfer_type() == "INTERP") {
    check_send_field_for_interpolation(sendField, fieldName);
  }
  else if (m_settings.get_transfer_type() == "COPY") {
    check_send_field_for_copy(sendField, fieldName);
  }
  else if (m_settings.get_transfer_type() == "PATCH") {
    check_send_field_for_patch_recovery(sendField, fieldName);
  }
  else {
    STK_ThrowErrorMsg("Invalid transfer type " + m_settings.get_transfer_type());
  }
}

stk::mesh::EntityRank get_required_send_field_rank(const std::string& transferType, const std::string& recvType)
{
  if (transferType == "INTERP")
  {
    return stk::topology::NODE_RANK;
  }
  if (transferType == "COPY")
  {
      if (recvType == "NODE") {
        return stk::topology::NODE_RANK;
      }
      else if (recvType == "ELEMENT_CENTROID" || recvType == "ELEMENT_GAUSS_POINT") {
        return stk::topology::ELEM_RANK;
      }
      else {
        STK_ThrowErrorMsg("COPY Transfer type is incompatible with receive type " + recvType);
        return stk::topology::INVALID_RANK;
      }
  }
  if (transferType == "PATCH")
  {
    return stk::topology::ELEM_RANK;
  }

  STK_ThrowErrorMsg("Transfer type " + transferType + " is not valid" );
  return stk::topology::INVALID_RANK;
}

stk::mesh::EntityRank get_required_recv_field_rank(const std::string& recvType)
{
  if (recvType == "NODE") {
    return stk::topology::NODE_RANK;
  }
  if (recvType == "ELEMENT_CENTROID" || recvType == "ELEMENT_GAUSS_POINT") {
    return stk::topology::ELEM_RANK;
  }
  if (recvType == "FACE_CENTROID" || recvType == "FACE_GAUSS_POINT") {
    return stk::topology::FACE_RANK;
  }
  if (recvType == "EDGE_CENTROID" || recvType == "EDGE_GAUSS_POINT") {
    return stk::topology::EDGE_RANK;
  }
    return stk::topology::INVALID_RANK;
}

void check_recv_field_for_interpolation(const stk::mesh::FieldBase* recvField, const std::string& recvType)
{
  STK_ThrowRequireMsg(recvField->entity_rank() == get_required_recv_field_rank(recvType),
                      "Specified recv field " + recvField->name() + " must have "
                      << get_required_recv_field_rank(recvType) << " for interpolation with recvType " + recvType);
}

void check_recv_field_for_copy(const stk::mesh::FieldBase* recvField, const std::string& recvType)
{
  STK_ThrowRequireMsg(recvField->entity_rank() == stk::topology::ELEMENT_RANK ||
                                                  recvField->entity_rank() == stk::topology::NODE_RANK,
                                         "Specified recv field " + recvField->name() + " must have NODE or ELEMENT rank for copy transfer");

  STK_ThrowRequireMsg(recvField->entity_rank() == get_required_recv_field_rank(recvType),
                      "Specified recv field " + recvField->name() + " must have "
                      << get_required_recv_field_rank(recvType) << " for copy with recvType " + recvType);
}

void check_recv_field_for_patch_recovery(const stk::mesh::FieldBase* recvField, const std::string& recvType) {
  check_recv_field_for_interpolation(recvField, recvType);
}

void TransferMainBroker::check_recv_field(const stk::mesh::FieldBase* recvField, const std::string& fieldName)
{
  if(stk::io::get_field_role(*recvField) != nullptr) {
    STK_ThrowRequireMsg(*stk::io::get_field_role(*recvField) == Ioss::Field::TRANSIENT,
                                          "Specified recv field " + recvField->name() + " must be TRANSIENT");
  }

  if (m_settings.get_transfer_type() == "INTERP") {
    check_recv_field_for_interpolation(recvField, fieldName);
  }
  else if (m_settings.get_transfer_type() == "COPY") {
    check_recv_field_for_copy(recvField, fieldName);
  }
  else if (m_settings.get_transfer_type() == "PATCH") {
    check_recv_field_for_patch_recovery(recvField, fieldName);
  }
  else {
    STK_ThrowErrorMsg("Invalid transfer type " + m_settings.get_transfer_type());
  }
}

void TransferMainBroker::check_and_create_fields_impl()
{
  std::vector<std::pair<std::string,std::string>> transfer_fields = m_settings.get_transfer_fields();

  if (transfer_fields.size() != 0) {
    for (auto fieldPair : transfer_fields) {
      stk::mesh::FieldBase* sendField = get_field_by_name(fieldPair.first, m_sendMeta);
      check_send_field(sendField, fieldPair.first);
      m_sendFieldsForTransfer.push_back(sendField);

      stk::mesh::FieldBase* recvField = get_field_by_name(fieldPair.second, m_recvMeta);

      if (recvField) {
        m_recvFieldsForTransfer.push_back(recvField);
      }
      else {
        stk::transfer_util::log_message(m_comm, "Specified recv field " + fieldPair.second + " not found; creating.");
        recvField = add_created_recvField(sendField, fieldPair.second);
      }
      check_recv_field(recvField, m_settings.get_recv_type_string());
    }
  }

  else {
    stk::mesh::EntityRank reqSendFieldRank = get_required_send_field_rank(m_settings.get_transfer_type(),
                                                                          m_settings.get_recv_type_string());
    auto sendFields = m_sendMeta.get_fields(reqSendFieldRank);
    for (auto field : sendFields) {
      const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
      if (*fieldRole == Ioss::Field::TRANSIENT) {
        m_sendFieldsForTransfer.push_back(field);
      }
    }

    STK_ThrowRequireMsg(m_sendFieldsForTransfer.size() > 0, "There are no fields on the send mesh eligible for transfer");

    for (auto sendField : m_sendFieldsForTransfer) {
      stk::mesh::EntityRank reqRecvFieldRank = get_required_recv_field_rank(m_settings.get_recv_type_string());
      auto recvField = m_recvMeta.get_field(reqRecvFieldRank, sendField->name());
      if (recvField) {
        STK_ThrowRequireMsg(*stk::io::get_field_role(*recvField) == Ioss::Field::TRANSIENT,
                                          "Recv field " + recvField->name() + " must be TRANSIENT");
        m_recvFieldsForTransfer.push_back(recvField);
      }
      else {
        auto recvFieldWrongRank = stk::mesh::get_field_by_name(sendField->name(), m_recvMeta);
        STK_ThrowRequireMsg(recvFieldWrongRank == nullptr, "Recv field " + recvFieldWrongRank->name() +
                                                           " exists but does not have compatible rank required for specified " + m_settings.get_recv_type_string() + " transfer type");

        stk::transfer_util::log_message(m_comm, "Send field " + sendField->name() +
                                               " not found on receive mesh; creating.");
        recvField = add_created_recvField(sendField);
      }
      check_recv_field(recvField, m_settings.get_recv_type_string());
    }
  }
}
void TransferMainBroker::check_and_create_fields()
{
  if (m_settings.get_transfer_type() != "INTERP" &&
      m_settings.get_transfer_type() != "COPY" &&
      m_settings.get_transfer_type() != "PATCH") {
    STK_ThrowErrorMsg("Invalid transfer type " + m_settings.get_transfer_type());
  }

  check_and_create_fields_impl();
}

void TransferMainBroker::check_part_names()
{
  for (auto sendPartName : m_sendPartNamesForTransfer) {
    STK_ThrowRequireMsg(m_sendMeta.get_part(sendPartName) != nullptr,
                        "Cannot complete transfer. Send part name " + sendPartName + " does not exist!");
  }

  for (auto recvPartName : m_recvPartNamesForTransfer) {
    STK_ThrowRequireMsg(m_recvMeta.get_part(recvPartName) != nullptr,
                        "Cannot complete transfer. Recv part name " + recvPartName + " does not exist!");
  }
}

bool is_compatible_send_part_rank_for_interpolation(const stk::mesh::Part* part)
{
  auto partRank = part->primary_entity_rank();
  return (partRank == stk::topology::ELEMENT_RANK);
}

bool is_compatible_recv_part_rank_for_interpolation(const stk::mesh::Part* part,
                                                    const std::string& recvType)
{
  auto partRank = part->primary_entity_rank();
  if (recvType == "NODE") {
    return (partRank == stk::topology::EDGE_RANK ||
            partRank == stk::topology::FACE_RANK ||
            partRank == stk::topology::ELEMENT_RANK);
  }
  if (recvType == "EDGE_CENTROID" || recvType == "EDGE_GAUSS_POINT") {
    return (partRank == stk::topology::EDGE_RANK);
  }
  if (recvType == "FACE_CENTROID" || recvType == "FACE_GAUSS_POINT") {
    return (partRank == stk::topology::FACE_RANK);
  }
  if (recvType == "ELEMENT_CENTROID" || recvType == "ELEMENT_GAUSS_POINT") {
    return (partRank == stk::topology::ELEMENT_RANK);
  }
  return false;
}

void TransferMainBroker::check_part_ranks_for_interpolation()
{
  for (auto sendPartName : m_sendPartNamesForTransfer) {
    STK_ThrowRequireMsg(is_compatible_send_part_rank_for_interpolation(m_sendMeta.get_part(sendPartName)),
                        "Cannot complete transfer. Send part name " + sendPartName + " does not have element rank!");
  }

  for (auto recvPartName : m_recvPartNamesForTransfer) {
    STK_ThrowRequireMsg(is_compatible_recv_part_rank_for_interpolation(m_recvMeta.get_part(recvPartName),
                        m_settings.get_recv_type_string()),
                        "Cannot complete transfer. Recv part name " + recvPartName + " does not have correct rank for transfer with receive type " + m_settings.get_recv_type_string() + "!");
  }
}

bool is_compatible_part_rank_for_copy(const stk::mesh::Part* part)
{
  auto partRank = part->primary_entity_rank();
  return (partRank == stk::topology::ELEMENT_RANK);
}

bool is_compatible_recv_part_rank_for_copy(const stk::mesh::Part* part,
                                                    const std::string& recvType)
{
  bool supported_recvType = (recvType == "NODE" || recvType == "ELEMENT_CENTROID" || recvType == "ELEMENT_GAUSS_POINT");
  return (is_compatible_part_rank_for_copy(part) && supported_recvType);
}

void TransferMainBroker::check_part_ranks_for_copy()
{
  for (auto sendPartName : m_sendPartNamesForTransfer) {
    STK_ThrowRequireMsg(is_compatible_part_rank_for_copy(m_sendMeta.get_part(sendPartName)),
                        "Cannot complete transfer. Send part name " + sendPartName + " does not have element rank!");
  }

  for (auto recvPartName : m_recvPartNamesForTransfer) {
    STK_ThrowRequireMsg(is_compatible_recv_part_rank_for_copy(m_recvMeta.get_part(recvPartName),
                        m_settings.get_recv_type_string()),
                        "Cannot complete transfer. Recv part name " + recvPartName + " does not have correct rank for transfer with receive type " + m_settings.get_recv_type_string() + "!");
  }
}

bool is_compatible_send_part_rank_for_patch_recovery(const stk::mesh::Part* part)
{
  auto partRank = part->primary_entity_rank();
  return (partRank == stk::topology::ELEMENT_RANK);
}

bool is_compatible_recv_part_rank_for_patch_recovery(const stk::mesh::Part* part,
                                                    const std::string& recvType)
{
  return is_compatible_recv_part_rank_for_interpolation(part, recvType);
}

void TransferMainBroker::check_part_ranks_for_patch_recovery()
{
  for (auto sendPartName : m_sendPartNamesForTransfer) {
    STK_ThrowRequireMsg(is_compatible_send_part_rank_for_patch_recovery(m_sendMeta.get_part(sendPartName)),
                        "Cannot complete transfer. Send part name " + sendPartName + " does not have element rank!");
  }

  for (auto recvPartName : m_recvPartNamesForTransfer) {
    STK_ThrowRequireMsg(is_compatible_recv_part_rank_for_patch_recovery(m_recvMeta.get_part(recvPartName),
                        m_settings.get_recv_type_string()),
                        "Cannot complete transfer. Recv part name " + recvPartName + " does not have correct rank for transfer with receive type " + m_settings.get_recv_type_string() + "!");
  }
}

void TransferMainBroker::check_part_ranks()
{
  if (m_settings.get_transfer_type() == "INTERP") {
    check_part_ranks_for_interpolation();
  }
  else if (m_settings.get_transfer_type() == "COPY") {
    check_part_ranks_for_copy();
  }
  else if (m_settings.get_transfer_type() == "PATCH") {
    check_part_ranks_for_patch_recovery();
  }
  else {
    STK_ThrowErrorMsg("Invalid transfer type " + m_settings.get_transfer_type());
  }
}

bool all_parts_same_rank(const std::vector<std::string>& partNames, const stk::mesh::MetaData& metaData)
{
  auto requiredPartRank = metaData.get_part(partNames[0])->primary_entity_rank();

  for (auto partName : partNames) {
    if (metaData.get_part(partName)->primary_entity_rank() != requiredPartRank) {
      return false;
    }
  }
  return true;
}

void TransferMainBroker::check_part_consistency()
{
  STK_ThrowRequireMsg(all_parts_same_rank(m_sendPartNamesForTransfer, m_sendMeta),
                      "Cannot complete transfer. Send parts are not all same rank!");

  STK_ThrowRequireMsg(all_parts_same_rank(m_recvPartNamesForTransfer, m_recvMeta),
                      "Cannot complete transfer. Recv parts are not all same rank!");
}

void TransferMainBroker::check_parts()
{

  check_part_names();
  check_part_ranks();
  check_part_consistency();
}

std::vector<std::string> TransferMainBroker::get_send_fields_for_transfer() const
{
  std::vector<std::string> sendFieldNamesForTransfer;
  for (auto field : m_sendFieldsForTransfer) {
    sendFieldNamesForTransfer.push_back(field->name());
  }
  return sendFieldNamesForTransfer;
}

std::vector<std::string> TransferMainBroker::get_recv_fields_for_transfer() const
{
  std::vector<std::string> recvFieldNamesForTransfer;
  for (auto field : m_recvFieldsForTransfer) {
    recvFieldNamesForTransfer.push_back(field->name());
  }
  return recvFieldNamesForTransfer;
}

std::vector<std::string> TransferMainBroker::get_send_parts_for_transfer() const
{
  return m_sendPartNamesForTransfer;
}

std::vector<std::string> TransferMainBroker::get_recv_parts_for_transfer() const
{
  return m_recvPartNamesForTransfer;
}

std::shared_ptr<stk::mesh::BulkData> TransferMainBroker::get_send_bulk()
{
  return m_sendBulk;
}

std::shared_ptr<stk::mesh::BulkData> TransferMainBroker::get_recv_bulk()
{
  return m_recvBulk;
}

void TransferMainBroker::set_master_element_provider(std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
{
  m_masterElementProvider = masterElemProvider;
}

void TransferMainBroker::transfer_initialize()
{
  m_simpleTransfer = std::make_shared<stk::transfer::spmd::SimpleTransfer>("TransferMainSimpleTransfer", m_comm);

  m_simpleTransfer->add_send_fields(m_sendFieldsForTransfer);
  m_simpleTransfer->add_recv_fields(m_recvFieldsForTransfer);

  m_simpleTransfer->add_send_part_names(m_sendPartNamesForTransfer);
  m_simpleTransfer->add_recv_part_names(m_recvPartNamesForTransfer);

  const unsigned recvMeshDim = m_recvBulk->mesh_meta_data().spatial_dimension();
  const unsigned sendMeshDim = m_sendBulk->mesh_meta_data().spatial_dimension();

  if (m_settings.get_enable_2d_to_3d_axisymmetric_transfer())
  {
    STK_ThrowRequireMsg(m_sendBulk->mesh_meta_data().spatial_dimension() == 2 &&
                        m_recvBulk->mesh_meta_data().spatial_dimension() == 3,
                        "2D to 3D axisymmetric transfer requires a 2D send mesh and a 3D receive mesh");
    const TwoDTo3DAxisymmetricParams& params = m_settings.get_2d_to_3d_axisymmetric_transfer_params();
    auto transform = std::make_shared<stk::search::CoordTransformAxisymmetric2D>(params.axis, params.trans);
    m_simpleTransfer->add_recv_mesh_coord_transform(transform);
  }
  else if (m_settings.get_enable_3d_to_3d_axisymmetric_transfer())
  {
    STK_ThrowRequireMsg(m_sendBulk->mesh_meta_data().spatial_dimension() == 3 &&
                        m_recvBulk->mesh_meta_data().spatial_dimension() == 3,
                        "3D to 3D axisymmetric transfer requires a 3D send mesh and a 3D receive mesh");
    const ThreeDTo3DAxisymmetricParams& params = m_settings.get_3d_to_3d_axisymmetric_transfer_params();
    auto transform = std::make_shared<stk::search::CoordTransformAxisymmetric3D>(params.theta_min, params.theta_max, params.axis, params.trans);
    m_simpleTransfer->add_recv_mesh_coord_transform(transform);
  }
  else if (recvMeshDim == 2 && sendMeshDim == 3) {
    auto coordTransform = std::make_shared<stk::search::CoordTransformAddZ>(m_settings.get_coord_transf_z_expr());
    m_simpleTransfer->add_recv_mesh_coord_transform(coordTransform);
  }
  else if (sendMeshDim == 2 && recvMeshDim == 3) {
    std::string mappingType = m_settings.get_2d_to_3d_mapping_type();
    if (mappingType == "EXTRUDE") {
      auto coordTransform = std::make_shared<stk::search::CoordTransformRemoveZ>();
      m_simpleTransfer->add_recv_mesh_coord_transform(coordTransform);
    }
    if (mappingType == "ZPLANE") {
      auto coordTransform = std::make_shared<stk::search::CoordTransformAddZ>(m_settings.get_coord_transf_z_expr());
      m_simpleTransfer->add_send_mesh_coord_transform(coordTransform);
    }
  }

  if (m_settings.get_transfer_type() == "INTERP") {
    m_simpleTransfer->setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                                    stk::topology::ELEM_RANK, m_settings.get_recv_type(),
                                                    m_masterElementProvider, m_settings.get_extrapolate_option());
  }
  else if (m_settings.get_transfer_type() == "COPY") {
    if (m_settings.get_recv_type_string() == "NODE") {
      m_simpleTransfer->setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                                    stk::topology::NODE_RANK, m_settings.get_recv_type(),
                                                    m_masterElementProvider, m_settings.get_extrapolate_option());
    }
    else if (m_settings.get_recv_type_string() == "ELEMENT_CENTROID" ||
             m_settings.get_recv_type_string() == "ELEMENT_GAUSS_POINT")
    {
      m_simpleTransfer->setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                                    stk::topology::ELEM_RANK, m_settings.get_recv_type(),
                                                    m_masterElementProvider, m_settings.get_extrapolate_option());
    }
    else {
      STK_ThrowErrorMsg("Unsupported receive type for COPY transfer type.");
    }
  }
  else if (m_settings.get_transfer_type() == "PATCH") {
    m_simpleTransfer->setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                                    stk::topology::ELEM_RANK, m_settings.get_recv_type(),
                                                    m_masterElementProvider, m_settings.get_extrapolate_option());
  }
  else {
    STK_ThrowErrorMsg("Invalid transfer type " + m_settings.get_transfer_type());
  }


  m_simpleTransfer->initialize();
}

void TransferMainBroker::transfer_apply()
{
  m_simpleTransfer->apply();
}

void TransferMainBroker::transfer()
{
  transfer_initialize();
  transfer_apply();
}

} // namespace transfer_util
} // namespace stk
