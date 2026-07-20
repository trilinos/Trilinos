/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAINBROKER_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAINBROKER_HPP

#include "stk_util/parallel/Parallel.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_search_util/MasterElementProvider.hpp"
#include "stk_transfer_util/TransferMainSettings.hpp"
#include "stk_transfer_util/spmd/SimpleTransfer.hpp"

namespace stk {
namespace transfer_util {

class TransferMainBroker
{
public:

  TransferMainBroker(MPI_Comm comm,
                     std::shared_ptr<stk::mesh::BulkData> sendBulk,
                     std::shared_ptr<stk::mesh::BulkData> recvBulk,
                     const TransferMainSettings& settings,
                     std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider = nullptr);

  std::vector<std::string> get_send_fields_for_transfer() const;
  std::vector<std::string> get_recv_fields_for_transfer() const;
  std::vector<std::string> get_send_parts_for_transfer() const;
  std::vector<std::string> get_recv_parts_for_transfer() const;

  void check_send_field(const stk::mesh::FieldBase* sendField, const std::string& fieldName);
  void check_recv_field(const stk::mesh::FieldBase* recvField, const std::string& fieldName);
  void check_and_create_fields_impl();
  void check_and_create_fields();
  
  void check_part_names();
  void check_part_types();
  void check_part_ranks_for_interpolation();
  void check_part_ranks_for_copy();
  void check_part_ranks_for_patch_recovery();
  void check_part_ranks();
  void check_part_consistency();
  void check_parts();

  void set_master_element_provider(std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider);

  void transfer_initialize();
  void transfer_apply();
  void transfer();

  std::shared_ptr<stk::mesh::BulkData> get_send_bulk();
  std::shared_ptr<stk::mesh::BulkData> get_recv_bulk();

private:

  stk::mesh::FieldBase* add_created_recvField(const stk::mesh::FieldBase* sendField, const std::string& fieldName = "");

  const MPI_Comm m_comm;
  std::shared_ptr<stk::mesh::BulkData> m_sendBulk;
  std::shared_ptr<stk::mesh::BulkData> m_recvBulk;
  stk::mesh::MetaData& m_sendMeta;
  stk::mesh::MetaData& m_recvMeta;
  const TransferMainSettings& m_settings;
  std::vector<stk::mesh::FieldBase*> m_sendFieldsForTransfer;
  std::vector<stk::mesh::FieldBase*> m_recvFieldsForTransfer;
  std::vector<std::string> m_sendPartNamesForTransfer;
  std::vector<std::string> m_recvPartNamesForTransfer;
  std::shared_ptr<stk::search::MasterElementProviderInterface> m_masterElementProvider;
  std::shared_ptr<stk::transfer::spmd::SimpleTransfer> m_simpleTransfer;
};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINBROKER_HPP
