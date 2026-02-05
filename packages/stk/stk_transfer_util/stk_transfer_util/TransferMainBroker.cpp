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
#include "stk_transfer_util/TransferMainBroker.hpp"
#include "stk_transfer_util/spmd/SimpleTransfer.hpp"
#include "stk_transfer_util/MockMasterElementProvider.hpp"
#include "stk_transfer_util/LogMessage.hpp"

namespace stk {
namespace transfer_util {

TransferMainBroker::TransferMainBroker(MPI_Comm comm, std::shared_ptr<stk::mesh::BulkData> sendBulk,
                                                      std::shared_ptr<stk::mesh::BulkData> recvBulk,
                                                      const TransferMainSettings& settings)
  : m_comm(comm),
    m_sendBulk(sendBulk),
    m_recvBulk(recvBulk),
    m_sendMeta(sendBulk->mesh_meta_data()),
    m_recvMeta(recvBulk->mesh_meta_data()),
    m_settings(settings)
{
  //extract all the element rank parts
  auto sendParts = m_sendMeta.get_parts();
  auto recvParts = m_recvMeta.get_parts();

  for (auto sendPart : sendParts) {
    if (stk::mesh::is_element_block(*sendPart)) {
      m_sendPartNamesForTransfer.push_back(sendPart->name());
    }
  }

  for (auto recvPart : recvParts) {
    if (stk::mesh::is_element_block(*recvPart)) {
      m_recvPartNamesForTransfer.push_back(recvPart->name());
    }
  }
}

void TransferMainBroker::add_created_recvField(const stk::mesh::FieldBase* sendField, const std::string& fieldName)
{
  m_recvMeta.enable_late_fields();
  auto createdRecvField = sendField->clone(m_recvMeta.get_field_repository(), fieldName);

  auto fieldOutputType = stk::io::get_field_output_type(*sendField);
  stk::io::set_field_output_type(*createdRecvField, fieldOutputType);
  stk::io::set_field_role(*createdRecvField, Ioss::Field::TRANSIENT);

  const stk::mesh::FieldRestrictionVector& oldRestrictions = sendField->restrictions();
  for(const stk::mesh::FieldRestriction& res : oldRestrictions)
  {
    stk::mesh::Selector selectNewParts = res.selector().clone_for_different_mesh(m_recvMeta);
    const std::byte* initValsPtr = createdRecvField->get_initial_value_bytes().extent(0) > 0 ? createdRecvField->get_initial_value_bytes().data() : nullptr;
    m_recvMeta.declare_field_restriction(*createdRecvField,
                                      selectNewParts,
                                      res.num_scalars_per_entity(),
                                      res.dimension(),
                                      initValsPtr);
  }
  m_recvFieldsForTransfer.push_back(createdRecvField); 

}

void TransferMainBroker::check_fields()
{
  std::vector<std::pair<std::string,std::string>> transfer_fields = m_settings.get_transfer_fields();
   
  if (transfer_fields.size() != 0) {
    for (auto fieldPair : transfer_fields) {
      stk::mesh::FieldBase* sendField = get_field_by_name(fieldPair.first, m_sendMeta);

      STK_ThrowRequireMsg(sendField != nullptr, "Specified field " + fieldPair.first + " was not found on the send mesh!");
      STK_ThrowRequireMsg(sendField->entity_rank() == stk::topology::NODE_RANK,
                                         "Specified send field " + fieldPair.first + " must have  NODE rank");
      STK_ThrowRequireMsg(*stk::io::get_field_role(*sendField) == Ioss::Field::TRANSIENT,
                                         "Specified send field " + fieldPair.first + " must be TRANSIENT");

      m_sendFieldsForTransfer.push_back(sendField);

      stk::mesh::FieldBase* recvField = get_field_by_name(fieldPair.second, m_recvMeta);

      if (recvField) {
        if (recvField->entity_rank() == stk::topology::NODE_RANK) {
          STK_ThrowRequireMsg(*stk::io::get_field_role(*recvField) == Ioss::Field::TRANSIENT,
                                          "Specified recv field " + fieldPair.second + " must be TRANSIENT");
          m_recvFieldsForTransfer.push_back(recvField);
        }
        else {
          stk::transfer_util::log_message(m_comm, "Found specified field " + fieldPair.second + 
                                                 " with incompatible (non-NODE) rank. Creating field " + fieldPair.second + 
                                                " of NODE rank to be transferred from the send mesh.");
          add_created_recvField(sendField, fieldPair.second);
        }
      }
      else {
        stk::transfer_util::log_message(m_comm, "Unable to find specified field " + fieldPair.second + 
                                                " on receive mesh. Creating field " + fieldPair.first + 
                                                ", to be transferred from the send mesh.");
        add_created_recvField(sendField);
      }
    }
  }

  else {
    auto sendFields = m_sendMeta.get_fields(stk::topology::NODE_RANK);
    for (auto field : sendFields) {
      const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
      if (*fieldRole == Ioss::Field::TRANSIENT) {
        m_sendFieldsForTransfer.push_back(field);
      } 
    }
    STK_ThrowRequireMsg(m_sendFieldsForTransfer.size() > 0, "There are no fields on the send mesh eligible for transfer");
    for (auto sendField : m_sendFieldsForTransfer) {
      auto recvField = m_recvMeta.get_field(stk::topology::NODE_RANK, sendField->name());
      if (recvField) {
        STK_ThrowRequireMsg(*stk::io::get_field_role(*recvField) == Ioss::Field::TRANSIENT,
                                          "Recv field " + recvField->name() + " must be TRANSIENT");
        m_recvFieldsForTransfer.push_back(recvField);
      }
      else {
        stk::transfer_util::log_message(m_comm, "Unable to find field " + sendField->name() + 
                                               " with NODE rank on receive mesh. Creating field of this name to be transferred from the send mesh.");
        add_created_recvField(sendField);
      }
    }
  }

  m_recvMeta.set_mesh_on_fields(m_recvBulk.get()); 
}

void TransferMainBroker::check_parts()
{
  for (auto sendPartName : m_sendPartNamesForTransfer) {
    STK_ThrowRequireMsg(stk::mesh::is_element_block(*(m_sendMeta.get_part(sendPartName))),
                        "Cannot complete transfer. Send part name " + sendPartName + " is not an element block!");
  }

  for (auto recvPartName : m_recvPartNamesForTransfer) {
    STK_ThrowRequireMsg(stk::mesh::is_element_block(*(m_recvMeta.get_part(recvPartName))),
                        "Cannot complete transfer. Recv part name " + recvPartName + " is not an element block!");
  }
}

std::vector<std::string> TransferMainBroker::get_recv_fields_for_transfer()
{ 
  std::vector<std::string> recvFieldNamesForTransfer;
  for (auto field : m_recvFieldsForTransfer) {
    recvFieldNamesForTransfer.push_back(field->name());
  }
  return recvFieldNamesForTransfer;
}

void TransferMainBroker::transfer()
{
  stk::transfer::spmd::SimpleTransfer simpleTransfer("TransferMainSimpleTransfer", m_comm);

  simpleTransfer.add_send_fields(m_sendFieldsForTransfer);
  simpleTransfer.add_recv_fields(m_recvFieldsForTransfer);

  simpleTransfer.add_send_part_names(m_sendPartNamesForTransfer);
  simpleTransfer.add_recv_part_names(m_recvPartNamesForTransfer);

  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider;
  masterElemProvider = std::make_shared<MasterElementProvider>();

  simpleTransfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                               stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                               masterElemProvider, m_settings.get_extrapolate_option());

  simpleTransfer.initialize();
  simpleTransfer.apply();
}

} // namespace transfer_util
} // namespace stk
