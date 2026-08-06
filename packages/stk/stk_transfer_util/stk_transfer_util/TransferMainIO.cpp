/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "stk_transfer_util/TransferMainIO.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/util/ParseCsv.hpp"
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
    m_recvMeta(m_recvBulk->mesh_meta_data()),
    m_outputBroker(),
    m_outputIndex(0)
{
}

void TransferMainIO::load_send_fields_at_time(double time)
{
  STK_ThrowAssertMsg(!m_sendBroker.get_time_steps().empty(),"Error, send mesh has no time steps");
  m_sendBroker.read_defined_input_fields(time);
}

void TransferMainIO::load_recv_fields_at_time(double time)
{
  STK_ThrowAssertMsg(!m_recvBroker.get_time_steps().empty(),"Error, recv mesh has no time steps");
  m_recvBroker.read_defined_input_fields(time);
}

void TransferMainIO::load_meshes()
{
  m_sendBroker.set_bulk_data(m_sendBulk);
  m_sendBroker.add_mesh_database(m_sendMeshFilename, stk::io::READ_MESH);
  m_sendBroker.create_input_mesh();
  m_sendBroker.add_all_mesh_fields_as_input_fields();
  m_sendBroker.populate_bulk_data();
  m_sendBroker.populate_field_data();
  m_sendBroker.cache_entity_list_for_transient_steps(true);

  m_recvBroker.set_bulk_data(m_recvBulk);
  m_recvBroker.add_mesh_database(m_recvMeshFilename, stk::io::READ_MESH);
  m_recvBroker.create_input_mesh();
  m_recvBroker.add_all_mesh_fields_as_input_fields();
  m_recvBroker.populate_bulk_data();
  m_recvBroker.populate_field_data();
  m_recvBroker.cache_entity_list_for_transient_steps(true);
}

std::vector<double> TransferMainIO::get_send_time_steps() const
{
  return m_sendBroker.get_time_steps();
}

std::vector<double> TransferMainIO::get_recv_time_steps() const
{
  return m_recvBroker.get_time_steps();
}

void TransferMainIO::initialize_transfer_output(const std::string& outputFile)
{
  m_outputBroker = std::make_shared<stk::io::StkMeshIoBroker>(m_comm);

  m_outputBroker->set_bulk_data(*m_recvBulk);
  m_outputBroker->set_attribute_field_ordering_stored_by_part_ordinal(m_recvBroker.get_attribute_field_ordering_stored_by_part_ordinal());

  m_outputIndex = m_outputBroker->create_output_mesh(outputFile, stk::io::WRITE_RESULTS);

  const stk::mesh::FieldVector& fields = m_recvBulk->mesh_meta_data().get_fields();
  for(stk::mesh::FieldBase* field : fields) {
    const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*field);
    if(fieldRole == nullptr || *fieldRole == Ioss::Field::TRANSIENT) {
      m_outputBroker->add_field(m_outputIndex, *field);
    }
  }

  m_outputBroker->cache_entity_list_for_transient_steps(true);
  m_outputBroker->write_output_mesh(m_outputIndex);
}

void TransferMainIO::write_transfer_output(int step, double time)
{
  STK_ThrowRequireMsg(step > 0, "Only write output steps greater than 0.");

  m_outputBroker->begin_output_step(m_outputIndex, time);
  m_outputBroker->write_defined_output_fields(m_outputIndex);
  m_outputBroker->end_output_step(m_outputIndex);
}

void TransferMainIO::write_transfer_output(const std::string& outputFile)
{
  initialize_transfer_output(outputFile);

  write_transfer_output(1, get_send_time_steps().back());
  m_outputBroker->flush_output();

  stk::transfer_util::log_message(m_comm, "Finished writing transfer output");
}

std::vector<double> get_times(const std::string& option, const std::vector<double>& times)
{
  if (sierra::make_upper(option) == "ALL") {
    return times;
  }
  if (sierra::make_upper(option) == "LAST") {
    return {times.back()};
  }
  if (sierra::make_upper(option) == "FIRST") {
    return {times.front()};
  }

  try {
    std::vector<int> steps = stk::util::get_ids_from_strings({option});
    stk::util::sort_and_unique(steps);
    std::vector<double> timeSteps;
    for(int step : steps) {
      int stepOffset = step-1;
      if (stepOffset >= 0 && stepOffset < static_cast<int>(times.size())) {
        timeSteps.push_back(times[stepOffset]);
      }
    }
    STK_ThrowRequireMsg(stk::util::is_sorted_and_unique(timeSteps),"Error, specified time-steps not sorted and unique.");

    return timeSteps;
  }
  catch(std::exception& e) {
    STK_ThrowErrorMsg(e.what());
  }

  return std::vector<double>();
}

} // namespace transfer_util
} // namespace stk
