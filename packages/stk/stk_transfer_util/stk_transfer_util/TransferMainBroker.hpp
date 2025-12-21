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

#include "mpi.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_transfer_util/TransferMainSettings.hpp"


namespace stk {
namespace transfer_util {

class TransferMainBroker
{
public:

  TransferMainBroker(MPI_Comm comm, std::shared_ptr<stk::mesh::BulkData> sendBulk,
                                    std::shared_ptr<stk::mesh::BulkData> recvBulk,
                                    const TransferMainSettings& settings);

  void check_fields();
  void check_parts();
  void transfer();
  std::vector<std::string> get_recv_fields_for_transfer();

private:

  void add_created_recvField(const stk::mesh::FieldBase* sendField, const std::string& fieldName = "");

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

};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINBROKER_HPP
