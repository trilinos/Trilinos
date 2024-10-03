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

#include "LifeCycle.hpp"

#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/EnvData.hpp"
#include "stk_util/environment/OutputLog.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_balance/mesh/BalanceMesh.hpp"
#include "stk_balance/io/BalanceIO.hpp"
#include "stk_balance/internal/Balancer.hpp"
#include "stk_balance/internal/LogUtils.hpp"
#include "stk_balance/rebalance.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_balance/internal/Diagnostics.hpp"
#include "stk_balance/internal/DiagnosticsPrinter.hpp"

#include <stk_io/FillMesh.hpp>

namespace stk {
namespace balance {

LifeCycle::LifeCycle(MPI_Comm c, int argc, const char** argv)
  : m_comm(c),
    m_argc(argc),
    m_argv(argv),
    m_exitCode(LifeCycleStatus::SUCCESS),
    m_isProc0(stk::parallel_machine_rank(m_comm) == 0),
    m_validator(m_comm),
    m_parser(m_comm)
{
  initialize_environment(m_comm, argv);

  try {
    parse();
  }
  catch(const std::exception &e) {
    print_parse_error(e.what());
    m_exitCode = LifeCycleStatus::PARSE_ERROR;
    return;
  }

  set_output_streams();
}

void LifeCycle::run()
{
  if (m_exitCode != 0) {
    return;
  }

  print_running_message();
  print_banner(sierra::Env::outputP0());

  if (is_serial_no_op()) {
    print_serial_no_op_message();
    return;
  }

  try {
    if (m_settings.get_is_rebalancing()) {
      rebalance();
    }
    else {
      balance();
    }
  }
  catch(std::exception& e) {
    print_balance_error(e.what());
    m_exitCode = LifeCycleStatus::EXECUTION_ERROR;
  }
}

LifeCycleStatus LifeCycle::exit_code() const
{
  return m_exitCode;
}

void LifeCycle::parse()
{
  m_parser.parse_command_line_options(m_argc, m_argv, m_settings);
  m_validator.require_file_exists(m_settings.get_input_filename(), m_settings.get_num_input_processors());
}

bool LifeCycle::is_serial_no_op() const
{
  return m_validator.serial_input_equals_output(m_settings.get_input_filename(), m_settings.get_output_filename()) &&
         (not m_settings.get_is_rebalancing());
}

bool LifeCycle::rebalance_will_corrupt_data(const stk::io::StkMeshIoBroker & ioBroker,
                                            const stk::mesh::MetaData & meta) const
{
  const bool isRebalancing =  m_settings.get_is_rebalancing();
  const bool sameInputAndOutputFile = m_validator.input_equals_output(m_settings.get_input_filename(), m_settings.get_output_filename());
  const bool sameInputAndOutputProcCount = m_settings.get_num_input_processors() == m_settings.get_num_output_processors();
  const bool hasTransientFieldData = (not stk::io::get_transient_fields(meta).empty());
  const bool has64BitIds = (ioBroker.check_integer_size_requirements_serial() > 4);

  bool willCorruptData = false;
  if (isRebalancing && sameInputAndOutputFile && sameInputAndOutputProcCount) {
    if (hasTransientFieldData) {
      willCorruptData = true;
      sierra::Env::outputP0()
          << "Aborting rebalance: Overwriting input files that contain transient fields will" << std::endl
          << "lead to data corruption.  Please specify a different outputDirectory." << std::endl;
    }

    if (has64BitIds) {
      willCorruptData = true;
      sierra::Env::outputP0()
          << "Aborting rebalance: Overwriting input files that contain 64-bit IDs will" << std::endl
          << "lead to data corruption.  Please specify a different outputDirectory." << std::endl;
    }
  }

  return willCorruptData;
}

void LifeCycle::balance()
{
  stk::balance::BalanceIO io(m_comm, m_settings);
  const stk::balance::Balancer balancer(m_settings);

  if (m_settings.shouldPrintDiagnostics()) {
    set_up_diagnostics(m_settings);
  }

  stk::balance::BalanceMesh& mesh = io.initial_decomp();
  balancer.balance(mesh);
  io.write(mesh);

  DiagnosticsPrinter diagPrinter(m_comm, m_settings.get_num_output_processors());
  diagPrinter.print(sierra::Env::outputP0());
}

void LifeCycle::rebalance()
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(m_comm).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  stk::io::StkMeshIoBroker ioBroker;

  meta.set_coordinate_field_name(m_settings.getCoordinateFieldName());
  stk::balance::internal::register_internal_fields_and_parts(*bulk, m_settings);
  stk::io::fill_mesh_preexisting(ioBroker, m_settings.get_input_filename(), *bulk);

  if (rebalance_will_corrupt_data(ioBroker, meta)) {
    m_exitCode = LifeCycleStatus::REBALANCE_CORRUPTION_ERROR;
    return;
  }

  if (m_settings.shouldPrintDiagnostics()) {
    set_up_diagnostics(m_settings);
  }

  stk::balance::rebalance(ioBroker, m_settings);

  DiagnosticsPrinter diagPrinter(m_comm, m_settings.get_num_output_processors());
  diagPrinter.print(sierra::Env::outputP0());
}

void LifeCycle::set_output_streams()
{
  if (m_isProc0) {
    const std::string & logName = m_settings.get_log_filename();
    if (logName == "cout"  || logName == "cerr") {
      stk::EnvData::instance().m_outputP0 = stk::get_log_ostream(logName);
    }
    else {
      stk::bind_output_streams("log=\"" + logName + "\"");
      stk::EnvData::instance().m_outputP0 = stk::get_log_ostream("log");
    }
  }
  else {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
  }

  Ioss::Utils::set_output_stream(sierra::Env::outputP0());
}

void LifeCycle::print_parse_error(const char* what) const
{
  if (m_isProc0) {
    std::cerr << what << std::endl;
  }
}

void LifeCycle::print_balance_error(const char* what) const
{
  if (m_isProc0) {
    std::cerr << what << std::endl;
  }
}

void LifeCycle::print_serial_no_op_message() const
{
  sierra::Env::outputP0()
    << "Aborting balance: Rewriting the same input serial mesh file.  Please specify" << std::endl
    << "a different outputDirectory if you wish to copy the input file to an output" << std::endl
    << "file of the same name." << std::endl;
}

void LifeCycle::print_running_message() const
{
  if (m_isProc0) {
    std::ostream diag_stream(std::cout.rdbuf());
    stk::register_ostream(diag_stream, "diag_stream");

    const std::string & logName = m_settings.get_log_filename();
    const bool usingLogFile = not (logName == "cout" || logName == "cerr");
    if (usingLogFile) {
      stk::bind_output_streams("diag_stream>log");
      stk::bind_output_streams("diag_stream>+cout");
    }
    else {
      stk::bind_output_streams("diag_stream>" + logName);
    }

    const unsigned inputRanks = m_settings.get_num_input_processors();
    const unsigned outputRanks = m_settings.get_num_output_processors();
    diag_stream << "stk_balance converting from " << inputRanks << " to " << outputRanks << " MPI ranks" << std::endl;

    if (usingLogFile) {
      diag_stream << "        Log file: " << logName << std::endl;
    }

    diag_stream << "     Input files:  "
                << construct_generic_parallel_file_name(m_settings.get_input_filename(), inputRanks) << std::endl;

    diag_stream << "    Output files:  "
                << construct_generic_parallel_file_name(m_settings.get_output_filename(), outputRanks) << std::endl;

    stk::unregister_ostream(diag_stream);
  }
}

} }
