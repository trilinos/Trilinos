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
#include "stk_util/util/ReportHandler.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

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
  //During parsing of command-line and/or input-file syntax, we will just output
  //"clean" error messages, meaning no stack-traces.
  //(These are the handlers used by macros such as STK_ThrowRequireMsg etc.)
  stk::ErrorHandler origAssertHandler = stk::set_assert_handler(stk::clean_error_handler);
  stk::ErrorHandler origErrorHandler = stk::set_error_handler(stk::clean_error_handler);
  stk::ErrorHandler origInvArgHandler = stk::set_invalid_arg_handler(stk::clean_error_handler);

  m_parser.parse_command_line_options(m_argc, m_argv);

  if (m_parser.get_parser_status() == TransferParserStatus::SUCCESS) {
    try {
      m_settings = m_parser.generate_transfer_settings();
      m_validator.require_file_exists(m_settings.get_sendMesh_filename(), m_settings.get_num_input_processors());
      m_validator.require_file_exists(m_settings.get_recvMesh_filename(), m_settings.get_num_output_processors());
    }
    catch(std::exception& e) {
      std::cerr << e.what() << std::endl;
      m_exitCode = TransferMainStatus::PARSE_ERROR;
    }
  }
  else if (m_parser.get_parser_status() == TransferParserStatus::PARSE_ONLY) {
    m_exitCode = TransferMainStatus::PARSE_ONLY;
  }
  else {
    m_exitCode = TransferMainStatus::PARSE_ERROR;
  }

  //Now restore the error handlers so that the STK_Throw* macros resume showing
  //stack-traces.
  stk::set_assert_handler(origAssertHandler);
  stk::set_error_handler(origErrorHandler);
  stk::set_invalid_arg_handler(origInvArgHandler);
}

bool TransferMainHandler::is_no_op() const
{
  return m_settings.get_sendMesh_filename() == m_settings.get_recvMesh_filename();
}

void TransferMainHandler::print_no_op_message() const
{
  stk::outputP0() << "Aborting transfer! Input and output mesh names are the same" << std::endl;
}

namespace {

void require_unique_recv_fields(const std::vector<std::unique_ptr<TransferMainBroker>>& brokers)
{
  std::map<std::string, unsigned> recvFieldToTransferIndex;
  for (unsigned i = 0; i < brokers.size(); ++i) {
    for (const std::string& recvFieldName : brokers[i]->get_recv_fields_for_transfer()) {
      auto [iter, inserted] = recvFieldToTransferIndex.emplace(recvFieldName, i);
      STK_ThrowRequireMsg(inserted, "Recv field " + recvFieldName + " is the target of more than one field group "
                                    "(groups " + std::to_string(iter->second + 1) + " and " + std::to_string(i + 1) + ").");
    }
  }
}

}

void TransferMainHandler::transfer()
{
  TransferMainIO io(m_comm, m_settings.get_sendMesh_filename(), m_settings.get_recvMesh_filename());
  io.load_meshes();

  // One broker (and therefore one SimpleTransfer) per field group, since each group can
  // have a different transfer type, receive type, etc. The brokers hold references to
  // their settings, so 'transferSettings' must outlive them and not be resized.
  const std::vector<TransferMainSettings> transferSettings = m_settings.get_settings_for_each_transfer();
  std::vector<std::unique_ptr<TransferMainBroker>> brokers;

  // All recv fields must be created before any transfer is initialized and before the
  // output mesh is set up.
  for (const TransferMainSettings& settings : transferSettings) {
    brokers.push_back(std::make_unique<TransferMainBroker>(m_comm, io.get_sendBulkData(), io.get_recvBulkData(),
                                                            settings, m_masterElementProvider));
    brokers.back()->check_parts();
    brokers.back()->check_and_create_fields();
  }

  require_unique_recv_fields(brokers);

  for (auto& broker : brokers) {
    broker->transfer_initialize();
  }

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

    for (auto& broker : brokers) {
      broker->transfer_apply();
    }

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
  if (m_parser.input_file_specified()) {
    stk::transfer_util::log_message(m_comm, "Input-file: "+m_parser.get_input_file_name());
  }
  stk::transfer_util::log_message(m_comm, "Send-mesh: "+ m_settings.get_sendMesh_filename());
  stk::transfer_util::log_message(m_comm, "Recv-mesh: "+ m_settings.get_recvMesh_filename());
  stk::transfer_util::log_message(m_comm, "Output-mesh: "+ m_settings.get_outputMesh_filename());
  const std::vector<TransferMainSettings> transferSettings = m_settings.get_settings_for_each_transfer();
  for (unsigned i = 0; i < transferSettings.size(); ++i) {
    const TransferMainSettings& settings = transferSettings[i];
    const std::string indent = m_settings.has_field_groups() ? "  " : "";
    if (m_settings.has_field_groups()) {
      stk::transfer_util::log_message(m_comm, "Field group " + std::to_string(i + 1) + ":");
    }
    stk::transfer_util::log_message(m_comm, indent + "Transfer type: "+ settings.get_transfer_type());
    stk::transfer_util::log_message(m_comm, indent + "Receive type: "+ settings.get_recv_type_string());
    stk::transfer_util::log_message(m_comm, indent + "Fields being transferred: "
                                            + (settings.get_transfer_fields().size() == 0 ? "all fields" :
                                               settings.get_field_list_string()));
    stk::transfer_util::log_message(m_comm, indent + "Send parts for transfer: "
                                            + (settings.get_transfer_send_parts().size() == 0 ? "all send parts" :
                                               settings.get_send_part_list_string()));
    stk::transfer_util::log_message(m_comm, indent + "Recv parts for transfer: "
                                            + (settings.get_transfer_recv_parts().size() == 0 ? "all recv parts" :
                                               settings.get_recv_part_list_string()));
    stk::transfer_util::log_message(m_comm, indent + "Extrapolate option: "+ settings.get_extrapolate_option_string());
  }
  stk::transfer_util::log_message(m_comm, "Time-steps: "+ m_settings.get_time_steps_spec());
  stk::transfer_util::log_message(m_comm, "MasterElements implementation: "+ m_settings.get_master_elements_name());
}

} }
