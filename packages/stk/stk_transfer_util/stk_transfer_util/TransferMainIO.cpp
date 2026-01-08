/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "stk_transfer_util/TransferMainIO.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_transfer_util/TransientFieldTransferById.hpp"
#include "stk_transfer_util/LogMessage.hpp"
#include "stk_io/WriteMesh.hpp"

namespace stk {
namespace transfer_util {

//only mesh names instead of settings
TransferMainIO::TransferMainIO(MPI_Comm comm, const std::string& sendMesh, const std::string& recvMesh)
  : m_comm(comm),
    m_sendMeshFilename(sendMesh),
    m_recvMeshFilename(recvMesh),
    m_sendBulk(stk::mesh::MeshBuilder(m_comm).create()),
    m_sendMeta(m_sendBulk->mesh_meta_data()),
    m_recvBulk(stk::mesh::MeshBuilder(m_comm).create()),
    m_recvMeta(m_recvBulk->mesh_meta_data())
{
}

void TransferMainIO::load_meshes()
{
  m_sendBroker.set_bulk_data(m_sendBulk);
  m_sendBroker.add_mesh_database(m_sendMeshFilename, stk::io::READ_MESH);
  m_sendBroker.create_input_mesh();
  m_sendBroker.add_all_mesh_fields_as_input_fields();
  m_sendBroker.populate_bulk_data();
  m_sendBroker.populate_field_data();

  const std::vector<double> sendTimeSteps = m_sendBroker.get_time_steps();
  if (!sendTimeSteps.empty()) {
    const double lastSendTimeStep = sendTimeSteps.back();
    m_sendBroker.read_defined_input_fields(lastSendTimeStep);
  }

  m_recvBroker.set_bulk_data(m_recvBulk);
  m_recvBroker.add_mesh_database(m_recvMeshFilename, stk::io::READ_MESH);
  m_recvBroker.create_input_mesh();
  m_recvBroker.add_all_mesh_fields_as_input_fields();
  m_recvBroker.populate_bulk_data();
  m_recvBroker.populate_field_data();

  const std::vector<double> recvTimeSteps = m_recvBroker.get_time_steps();
  if (!recvTimeSteps.empty()) {
    const double lastRecvTimeStep = recvTimeSteps.back();
    m_recvBroker.read_defined_input_fields(lastRecvTimeStep);
  }
}

void TransferMainIO::write_transfer_output( std::vector<std::string> transferredFields, std::string outputFile)
{
  const std::vector<double> sendTimeSteps = m_sendBroker.get_time_steps();
  stk::io::write_mesh_with_specified_fields(outputFile, *m_recvBulk, transferredFields, 1, sendTimeSteps.back());
  stk::transfer_util::log_message(m_comm, "Finished writing output mesh");
}

} // namespace transfer_util
} // namespace stk
