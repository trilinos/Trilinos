/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAINIO_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAINIO_HPP

#include "mpi.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

#include "stk_transfer_util/TransferMainSettings.hpp"


namespace stk {
namespace transfer_util {

class TransferMainIO
{
public:

  TransferMainIO(MPI_Comm comm, const std::string& sendMesh, const std::string& recvMesh);

  void load_meshes();

  std::shared_ptr<stk::mesh::BulkData> get_sendBulkData() { return m_sendBulk; }
  std::shared_ptr<stk::mesh::BulkData> get_recvBulkData() { return m_recvBulk; }

  void write_transfer_output(std::vector<std::string> transferredFields, std::string outputFile = "transferredReceive.exo");

private:
  const MPI_Comm m_comm;
  const std::string m_sendMeshFilename;
  const std::string m_recvMeshFilename;

  std::shared_ptr<stk::mesh::BulkData> m_sendBulk;
  stk::mesh::MetaData& m_sendMeta;

  std::shared_ptr<stk::mesh::BulkData> m_recvBulk;
  stk::mesh::MetaData& m_recvMeta;

  stk::io::StkMeshIoBroker m_sendBroker;
  stk::io::StkMeshIoBroker m_recvBroker;

};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINIO_HPP
