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

#include "stk_transfer_util/TransferMainHandler.hpp"
#include "stk_transfer_util/TransferMainIO.hpp"
#include "stk_transfer_util/TransferMainBroker.hpp"
#include "stk_transfer_util/LogMessage.hpp"

namespace stk {
namespace transfer_util {

TransferMainHandler::TransferMainHandler(MPI_Comm c, int argc, const char** argv)
  : m_comm(c),
    m_argc(argc),
    m_argv(argv),
    m_exitCode(TransferMainStatus::SUCCESS),
    m_isProc0(stk::parallel_machine_rank(m_comm) == 0),
    m_validator(m_comm),
    m_parser(m_comm)
{
  parse();
}

void TransferMainHandler::run()
{
  if (m_exitCode != TransferMainStatus::SUCCESS) {
    return;
  }

  print_running_message();

  if (is_no_op()) {
    print_no_op_message();
    return;
  }

  try {
    transfer();
  }
  catch(std::exception& e) {
    print_transfer_error(e.what());
    m_exitCode = TransferMainStatus::EXECUTION_ERROR;
  }
}

TransferMainStatus TransferMainHandler::exit_code() const
{
  return m_exitCode;
}

void TransferMainHandler::set_master_element_provider(std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
{
  m_masterElementProvider = masterElemProvider;
}

void TransferMainHandler::parse()
{
  m_parser.parse_command_line_options(m_argc, m_argv);

  if (m_parser.get_parser_status() == TransferParserStatus::SUCCESS) {
    m_settings = m_parser.generate_transfer_settings();
    m_validator.require_file_exists(m_settings.get_sendMesh_filename(), m_settings.get_num_input_processors());
    m_validator.require_file_exists(m_settings.get_recvMesh_filename(), m_settings.get_num_output_processors());
  }
  else if (m_parser.get_parser_status() == TransferParserStatus::PARSE_ONLY) {
    m_exitCode = TransferMainStatus::PARSE_ONLY;
  }
  else {
    m_exitCode = TransferMainStatus::PARSE_ERROR;
  }
}

bool TransferMainHandler::is_no_op() const
{
  return m_settings.get_sendMesh_filename() == m_settings.get_recvMesh_filename();
}

void TransferMainHandler::print_no_op_message() const
{
  stk::outputP0() << "Aborting transfer! Input and output mesh names are the same" << std::endl;
}

void TransferMainHandler::transfer()
{
  TransferMainIO io(m_comm, m_settings.get_sendMesh_filename(), m_settings.get_recvMesh_filename());
  io.load_meshes();

  TransferMainBroker broker(m_comm, io.get_sendBulkData(), io.get_recvBulkData(), m_settings, m_masterElementProvider);
  
  broker.check_parts();
  broker.check_and_create_fields();
 
  broker.transfer_initialize();

  const std::vector<double> sendTimeSteps = io.get_send_time_steps();
  const std::vector<double> recvTimeSteps = io.get_recv_time_steps();
  STK_ThrowRequireMsg(recvTimeSteps.empty() || sendTimeSteps.size() == recvTimeSteps.size(),
      "The recv mesh is required to either have no time steps or have the same number of time steps as the send mesh.");

  const bool recvMeshHasTimeSteps = sendTimeSteps.size() == recvTimeSteps.size();

  io.initialize_transfer_output(m_settings.get_outputMesh_filename());

  const std::vector<double> xferTimeSteps = get_times(m_settings.get_time_steps_spec(), sendTimeSteps);

  for(unsigned i=0; i<xferTimeSteps.size(); ++i) {
    int step = i + 1;
    double time = xferTimeSteps[i];

    io.load_send_fields_at_time(time);

    if (recvMeshHasTimeSteps) {
      io.load_recv_fields_at_time(time);
    }

    broker.transfer_apply();

    io.write_transfer_output(step, time);
  }

  stk::transfer_util::log_message(m_comm, "Finished writing output mesh.");
}

void TransferMainHandler::print_parse_error(const char* what) const
{
  if (m_isProc0) {
    std::cerr << what << std::endl;
  }
}

void TransferMainHandler::print_transfer_error(const char* what) const
{
  if (m_isProc0) {
    std::cerr << what << std::endl;
  }
}

void TransferMainHandler::print_running_message() const
{
  stk::transfer_util::log_message(m_comm, "STK Transfer: STK Version "+stk::version_string());
  stk::transfer_util::log_message(m_comm, "Transfer Information:");
  stk::transfer_util::log_message(m_comm, "Send-mesh: "+ m_settings.get_sendMesh_filename());
  stk::transfer_util::log_message(m_comm, "Recv-mesh: "+ m_settings.get_recvMesh_filename());
  stk::transfer_util::log_message(m_comm, "Output-mesh: "+ m_settings.get_outputMesh_filename());
  stk::transfer_util::log_message(m_comm, "Transfer type: "+ m_settings.get_transfer_type());
  stk::transfer_util::log_message(m_comm, "Receive type: "+ m_settings.get_recv_type_string());
  stk::transfer_util::log_message(m_comm, "Fields being transferred: "
                                          + (m_settings.get_transfer_fields().size() == 0 ? "all fields" : 
                                             m_settings.get_field_list_string()));
  stk::transfer_util::log_message(m_comm, "Send parts for transfer: "
                                          + (m_settings.get_transfer_send_parts().size() == 0 ? "all send parts" : 
                                             m_settings.get_send_part_list_string()));
  stk::transfer_util::log_message(m_comm, "Recv parts for transfer: "
                                          + (m_settings.get_transfer_recv_parts().size() == 0 ? "all recv parts" : 
                                             m_settings.get_recv_part_list_string()));
  stk::transfer_util::log_message(m_comm, "Time-steps: "+ m_settings.get_time_steps_spec());
  stk::transfer_util::log_message(m_comm, "Extrapolate option: "+ m_settings.get_extrapolate_option_string());
  stk::transfer_util::log_message(m_comm, "MasterElements implementation: "+ m_settings.get_master_elements_name());
}

} }
