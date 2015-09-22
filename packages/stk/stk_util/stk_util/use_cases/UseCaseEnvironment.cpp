// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// 

#include <iostream>
#include <iomanip>

#include <stk_util/diag/Option.hpp>
#include <stk_util/util/Writer.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include <stk_util/util/Bootstrap.hpp>
#include <stk_util/util/IndentStreambuf.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

namespace {

namespace bopt = boost::program_options;

// Build output logging description for binding output streams
std::string
build_log_description(
  const bopt::variables_map &   vm,
  const std::string &           working_directory,
  int                           parallel_rank,
  int                           parallel_size)
{
  std::ostringstream output_description;

  // On processor 0:
  //   [outfile=path] [poutfile=path.n.r] [doutfile=path.n.r] out>{-|cout|cerr|outfile}+pout pout>{null|poutfile} dout>{out|doutfile}

  // On processor 1..n:
  //   [poutfile=path.n.r] [doutfile=path.n.r] out>pout pout>{null|poutfile} dout>{out|doutfile}

  std::string out_path = "-";
  if (vm.count("output-log"))
    out_path = vm["output-log"].as<std::string>();
  if (out_path == "-")
    out_path = "cout";

  std::string out_ostream;

  if (!stk::get_log_ostream(out_path))
    if (out_path.size() && out_path[0] != '/')
      out_path = working_directory + out_path;

  if (parallel_rank == 0) {
    if (!stk::get_log_ostream(out_path)) {
      output_description << "outfile=\"" << out_path << "\"";
      out_ostream = "outfile";
    }
    else
      out_ostream = out_path;
  }
  else
    out_ostream = "null";

  std::string pout_ostream = "null";
  if (vm.count("pout")) {
    std::string pout_path = vm["pout"].as<std::string>();
    if (pout_path == "-") {
      std::ostringstream s;

      if (stk::get_log_ostream(out_path))
        s << working_directory << "sierra.log." << parallel_size << "." << parallel_rank;
      else
        s << out_path << "." << parallel_size << "." << parallel_rank;
      pout_path = s.str();
    }
    else if (pout_path.find("/") == std::string::npos && !stk::get_log_ostream(pout_path)) {
      std::ostringstream s;

      s << working_directory << pout_path << "." << parallel_size << "." << parallel_rank;
      pout_path = s.str();
    }

    if (!stk::get_log_ostream(pout_path)) {
      output_description << " poutfile=\"" << pout_path << "\"";
      pout_ostream = "poutfile";
    }
    else
      pout_ostream = pout_path;
  }

  std::string dout_ostream;
  if (vm.count("dout")) {
    std::string dout_path = vm["dout"].as<std::string>();
    if (!dout_path.empty() && stk::is_registered_ostream(dout_path))
      dout_ostream = dout_path;
    else {
      std::ostringstream s;
      if (dout_path.size() && dout_path[0] != '/')
        s << working_directory << dout_path << "." << parallel_size << "." << parallel_rank;
      else
        s << dout_path << parallel_size << "." << parallel_rank;
      dout_path = s.str();
      output_description << " doutfile=\"" << dout_path << "\"";
      dout_ostream = "doutfile";
    }
  }
  else
    dout_ostream = "out";

  if (parallel_rank == 0)
    output_description << " out>" << out_ostream << "+pout";
  else
    output_description << " out>pout";

  output_description << " pout>" << pout_ostream << " dout>" << dout_ostream;

  return output_description.str();
}

stk::diag::OptionMaskParser dw_option_mask;
stk::diag::OptionMaskParser timer_option_mask;

void
stk_bootstrap()
{
  dw_option_mask.mask("search", use_case::LOG_SEARCH, "log search diagnostics");
  dw_option_mask.mask("transfer", use_case::LOG_TRANSFER, "log transfer diagnostics");
  dw_option_mask.mask("timer", use_case::LOG_TIMER, "log timer diagnostics");

  timer_option_mask.mask("mesh", use_case::TIMER_MESH, "mesh operations timers");
  timer_option_mask.mask("meshio", use_case::TIMER_MESH_IO, "mesh I/O timers");
  timer_option_mask.mask("transfer", use_case::TIMER_TRANSFER, "transfer timers");
  timer_option_mask.mask("search", use_case::TIMER_SEARCH, "search timers");

  boost::program_options::options_description desc("Use case environment options");
  desc.add_options()
    ("help,h", "produce help message")
    ("directory,d", boost::program_options::value<std::string>(), "working directory")
    ("output-log,o", boost::program_options::value<std::string>(), "output log path")
    ("pout", boost::program_options::value<std::string>()->implicit_value("-"), "per-processor log file path")
    ("dout", boost::program_options::value<std::string>()->implicit_value("out"), "diagnostic output stream one of: 'cout', 'cerr', 'out' or a file path")
    ("runtest,r", boost::program_options::value<std::string>(), "runtest pid file");

  stk::get_options_description().add(desc);
}

} // namespace <empty>

namespace use_case {

// Output streams
std::ostream &
out() {
  static std::ostream s_out(std::cout.rdbuf());

  return s_out;
}


std::ostream &
pout() {
  static std::ostream s_pout(std::cout.rdbuf());

  return s_pout;
}


std::ostream &
dout() {
  static std::ostream s_dout(std::cout.rdbuf());

  return s_dout;
}


std::ostream &
tout() {
  static std::ostream s_tout(std::cout.rdbuf());

  return s_tout;
}


std::ostream &
dwout() {
  static stk::indent_streambuf s_dwoutStreambuf(std::cout.rdbuf());
  static std::ostream s_dwout(&s_dwoutStreambuf);

  return s_dwout;
}


// Diagnostic writer
stk::diag::Writer &
dw()
{
  static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);

  return s_diagWriter;
}


// Message reporting
std::ostream &
operator<<(
  std::ostream &	os,
  message_type           type)
{
  switch (type & stk::MSG_TYPE_MASK) {
  case MSG_WARNING:
    os << "Warning";
    break;
  case MSG_FATAL:
    os << "Fatal error";
    break;
  case MSG_INFORMATION:
    os << "Information";
    break;
  case MSG_EXCEPTION:
    os << "Exception";
    break;
  case MSG_PARALLEL_EXCEPTION:
    os << "Parallel exception";
    break;
  }
  return os;
}


void
report_handler(
  const char *		message,
  int                   type)
{
  if (type & stk::MSG_DEFERRED)
    pout() << "Deferred " << static_cast<message_type>(type) << ": " << message << std::endl;

  else
    out() << static_cast<message_type>(type) << ": " << message << std::endl;
}


// Timers
stk::diag::TimerSet &
timerSet()
{
  static stk::diag::TimerSet s_timerSet(TIMER_ALL);

  return s_timerSet;
}


stk::diag::Timer &timer() {
  static stk::diag::Timer s_timer = stk::diag::createRootTimer("Use Cases", timerSet());

  return s_timer;
}


UseCaseEnvironment::UseCaseEnvironment(
  int *         argc,
  char ***      argv)
  : m_comm(stk::parallel_machine_init(argc, argv)),
    m_need_to_finalize(true)
{
  initialize(argc, argv);
}

UseCaseEnvironment::UseCaseEnvironment(
  int *         argc,
  char ***      argv,
  stk::ParallelMachine comm)
  : m_comm(comm),
    m_need_to_finalize(false)
{
  initialize(argc, argv);
}

void UseCaseEnvironment::initialize(int* argc, char*** argv)
{
  stk::register_log_ostream(std::cout, "cout");
  stk::register_log_ostream(std::cerr, "cerr");

  stk::register_ostream(out(), "out");
  stk::register_ostream(pout(), "pout");
  stk::register_ostream(dout(), "dout");
  stk::register_ostream(tout(), "tout");

  static_cast<stk::indent_streambuf *>(dwout().rdbuf())->redirect(dout().rdbuf());

  stk::set_report_handler(report_handler);

  stk::Bootstrap::bootstrap();
  stk_bootstrap();

  for (int i = 0; i < *argc; ++i) {
    const std::string s((*argv)[i]);
    if (s == "-h" || s == "-help" || s == "--help") {
      std::cout << "Usage: " << (*argv)[0] << " [options...]" << std::endl;
      std::cout << stk::get_options_description() << std::endl;
      return; // So application can handle app-specific options.
    }
  }

  // Broadcast argc and argv to all processors.
  int parallel_rank = stk::parallel_machine_rank(m_comm);
  int parallel_size = stk::parallel_machine_size(m_comm);

  stk::BroadcastArg b_arg(m_comm, *argc, *argv);

  // Parse broadcast arguments
  bopt::variables_map &vm = stk::get_variables_map();
  try {
    bopt::store(bopt::command_line_parser(b_arg.m_argc, b_arg.m_argv).options(stk::get_options_description()).allow_unregistered().run(), vm);
    bopt::notify(vm);
  }
  catch (std::exception &x) {
    stk::RuntimeDoomedSymmetric() << x.what();
  }

  // Set working directory
  m_workingDirectory = "./";
  if (vm.count("directory"))
    m_workingDirectory = vm["directory"].as<std::string>();
  if (m_workingDirectory.length() && m_workingDirectory[m_workingDirectory.length() - 1] != '/')
    m_workingDirectory += "/";

  std::string output_description = build_log_description(vm, m_workingDirectory, parallel_rank, parallel_size);

  stk::bind_output_streams(output_description);

  dout() << "Output log binding: " << output_description << std::endl;

  // Start use case root timer
  timer().start();
}

UseCaseEnvironment::~UseCaseEnvironment()
{
  stk::report_deferred_messages(m_comm);

// Stop use case root timer
  timer().stop();

  stk::diag::printTimersTable(out(), timer(), stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME, false, m_comm);

  stk::diag::deleteRootTimer(timer());

  static_cast<stk::indent_streambuf *>(dwout().rdbuf())->redirect(std::cout.rdbuf());

  stk::unregister_ostream(tout());
  stk::unregister_ostream(dout());
  stk::unregister_ostream(pout());
  stk::unregister_ostream(out());

  stk::unregister_log_ostream(std::cerr);
  stk::unregister_log_ostream(std::cout);

  if (m_need_to_finalize) {
    stk::parallel_machine_finalize();
  }
}

bool print_status(stk::ParallelMachine comm, bool success)
{
  int error_flag = success ? 0 : 1;
  stk::all_reduce( comm , stk::ReduceMax<1>( & error_flag ) );
  bool all_success = !error_flag;

  int rank = stk::parallel_machine_rank(comm);
  if (rank == 0) {
    std::cout << ( all_success ? "STKUNIT_ALL_PASS" : "Use case failed.") << std::endl;
  }

  return all_success;
}

} // namespace use_case
