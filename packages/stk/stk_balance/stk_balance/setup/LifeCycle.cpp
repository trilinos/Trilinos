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

#include "stk_balance/mesh/BalanceMesh.hpp"
#include "stk_balance/io/BalanceIO.hpp"
#include "stk_balance/internal/Balancer.hpp"

namespace stk {
namespace balance {

LifeCycle::LifeCycle(MPI_Comm c, int argc, const char** argv)
  : m_comm(c),
    m_argc(argc),
    m_argv(argv),
    m_exitCode(0),
    m_isProc0(stk::parallel_machine_rank(m_comm) == 0),
    m_validator(m_comm),
    m_parser(m_comm)
{
  set_output_streams();
}

void LifeCycle::run()
{
  try {
    parse();
  }
  catch(const std::exception &e) {
    print_parse_error(e.what());
    m_exitCode = 1;
    return;
  }

  if (serial_no_op()) {
    print_no_op_message();
    return;
  }

  try {
    balance();
  }
  catch(std::exception& e) {
    print_balance_error(e.what());
    m_exitCode = 2;
  }
}

int LifeCycle::exit_code() const
{
  return m_exitCode;
}

void LifeCycle::parse()
{
  m_parser.parse_command_line_options(m_argc, m_argv, m_settings);
  m_validator.require_file_exists(m_settings.get_input_filename());
}

bool LifeCycle::serial_no_op() const
{
  return m_validator.serial_input_equals_output(m_settings.get_input_filename(), m_settings.get_output_filename());
}

void LifeCycle::balance()
{
  print_running_message();
  stk::balance::BalanceIO io(m_comm, m_settings);
  const stk::balance::Balancer balancer(m_settings);

  stk::balance::BalanceMesh& mesh = io.initial_decomp();
  balancer.balance(mesh);
  io.write(mesh);
}

void LifeCycle::set_output_streams()
{
  if (!m_isProc0) {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
  }
  Ioss::Utils::set_output_stream(sierra::Env::outputP0());
}

void LifeCycle::print_parse_error(const char* what) const
{
  if (m_isProc0) {
    std::cerr << what << m_parser.get_quick_error() << std::endl;
  }
}

void LifeCycle::print_balance_error(const char* what) const
{
  if (m_isProc0) {
    std::cerr << what << std::endl;
  }
}

void LifeCycle::print_no_op_message() const
{
  sierra::Env::outputP0() << "Running on 1 MPI rank and input-file (" << m_settings.get_input_filename()
    << ") == output-file (" << m_settings.get_output_filename()
    << "), doing nothing. Specify outputDirectory if you "
    << "wish to copy the input-file to an output-file of the same name." << std::endl;
}

void LifeCycle::print_running_message() const
{
  sierra::Env::outputP0() << "Running stk_balance" << std::endl;
  sierra::Env::outputP0() << "  Input file:  " << m_settings.get_input_filename() << std::endl;
  sierra::Env::outputP0() << "  Output file: " << m_settings.get_output_filename() << std::endl;
}

} }
