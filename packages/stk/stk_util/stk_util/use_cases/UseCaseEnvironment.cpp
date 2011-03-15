/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>

#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include <stk_util/util/Bootstrap.hpp>
#include <stk_util/util/IndentStreambuf.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

namespace {

namespace bopt = boost::program_options;

// Parse command line bit masks and produce -h documentation. (Probably moved to Util at some point)
typedef unsigned long OptionMask;

struct OptionMaskName
{
  OptionMaskName()
    : m_name(""),
      m_mask(0),
      m_description("")
  {}

  OptionMaskName(const std::string &name, const OptionMask &mask, const std::string &description = "No description available")
    : m_name(name),
      m_mask(mask),
      m_description(description)
  {}

  virtual ~OptionMaskName()
  {}

  std::string		m_name;
  OptionMask		m_mask;
  std::string		m_description;
};


class OptionMaskNameMap: public std::map<std::string, OptionMaskName>
{
public:
  void mask(const std::string &name, const OptionMask mask, const std::string &description) {
    iterator it = find(name);
    if (it == end())
      insert(std::make_pair(name, OptionMaskName(name, mask, description)));
    else {
      (*it).second.m_mask = mask;
      (*it).second.m_description = description;
    }
  }
};

class OptionMaskParser
{
public:
  typedef OptionMask Mask;		///< Mask for this option

public:
  /**
   * Creates a new <b>OptionMaskParser</b> instance.
   *
   */
  OptionMaskParser(const std::string &description)
    : m_optionMaskNameMap(),
      m_description(description),
      m_optionMask(0),
      m_status(true)
  {}

  virtual ~OptionMaskParser()
  {}

  Mask parse(const char *mask) const;

  virtual void parseArg(const std::string &name) const;

  std::string describe() const {
    std::ostringstream strout;
    strout << m_description << std::endl;
    for (OptionMaskNameMap::const_iterator it = m_optionMaskNameMap.begin(); it != m_optionMaskNameMap.end(); ++it)
      strout << "  " << (*it).first << std::setw(14 - (*it).first.size()) << " " << (*it).second.m_description << std::endl;
    return strout.str();
  }

  void mask(const std::string &name, const Mask mask, const std::string &description) {
    m_optionMaskNameMap.mask(name, mask, description);
  }

protected:
  OptionMaskNameMap		m_optionMaskNameMap;	///< Mask name vector
  std::string                   m_description;          ///< Help description
  mutable OptionMask		m_optionMask;		///< Most recently parsed mask
  mutable bool			m_status;		///< Result of most recent parse
};


OptionMaskParser::Mask
OptionMaskParser::parse(
  const char *          mask) const
{
  if (mask) {
    const std::string mask_string(mask);

    m_status = true;

    std::string::const_iterator it0 = mask_string.begin();
    std::string::const_iterator it1;
    std::string::const_iterator it2;
    std::string::const_iterator it3;
    do {
      // Trim preceeding spaces
      while (it0 != mask_string.end() && *it0 == ' ')
        it0++;

      if (it0 == mask_string.end())
        break;

      for (it1 = it0; it1 != mask_string.end(); ++it1) {
        if (*it1 == '(' || *it1 == ':' || *it1 == ',')
          break;
      }

      // Trim trailing spaces
      it2 = it1;
      while (it2 != it0 && *(it2 - 1) == ' ')
        --it2;

      std::string name(it0, it2);

      // Get argument list
      if (*it1 == '(') {
        it2 = it1 + 1;

        // Trim preceeding spaces
        while (it2 != mask_string.end() && *it2 == ' ')
          ++it2;

        int paren_count = 0;

        for (; it1 != mask_string.end(); ++it1) {
          if (*it1 == '(')
            ++paren_count;
          else if (*it1 == ')') {
            --paren_count;
            if (paren_count == 0)
              break;
          }
        }
        it3 = it1;

        // Trim trailing spaces
        while (it3 != it2 && *(it3 - 1) == ' ')
          --it3;

        // Find next argument start
        for (; it1 != mask_string.end(); ++it1)
          if (*it1 == ':' || *it1 == ',')
            break;
      }
      else
        it2 = it3 = it1;

      const std::string arg(it2, it3);

      parseArg(name);

      it0 = it1 + 1;
    } while (it1 != mask_string.end());
  }

  return m_optionMask;
}


void
OptionMaskParser::parseArg(
  const std::string &	name) const
{
  OptionMaskNameMap::const_iterator mask_entry = m_optionMaskNameMap.find(name);

  if (mask_entry != m_optionMaskNameMap.end()) m_optionMask |= (*mask_entry).second.m_mask;
  else {
    Mask	mask_hex = 0;
    std::istringstream mask_hex_stream(name.c_str());
    if (mask_hex_stream >> std::resetiosflags(std::ios::basefield) >> mask_hex)
      m_optionMask |= mask_hex;
    else
      m_status = false;
  }
}

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

OptionMaskParser dw_option_mask("use case diagnostic writer");
OptionMaskParser timer_option_mask("use case timers");

void
bootstrap()
{
  /// \todo REFACTOR  Put these program_options in a function
  ///                 that can be called without the bootstrapping.
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
    ("dw", boost::program_options::value<std::string>(), dw_option_mask.describe().c_str())
    ("timer", boost::program_options::value<std::string>(), timer_option_mask.describe().c_str())
    ("runtest,r", boost::program_options::value<std::string>(), "runtest pid file");

  stk::get_options_description().add(desc);
}

stk::Bootstrap x(bootstrap);

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
    pout() << "Deferred " << (message_type) type << ": " << message << std::endl;

  else
    out() << (message_type) type << ": " << message << std::endl;
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

  // Parse diagnostic messages to display
  if (vm.count("dw"))
    dw().setPrintMask(dw_option_mask.parse(vm["dw"].as<std::string>().c_str()));

  // Parse timer metrics and classes to display
  stk::diag::setEnabledTimerMetricsMask(stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME);
  if (vm.count("timer"))
    timerSet().setEnabledTimerMask(timer_option_mask.parse(vm["timer"].as<std::string>().c_str()));

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
    std::cout << ( all_success ? "STK_USECASE_PASS" : "Use case failed.") << std::endl;
  }

  return all_success;
}

} // namespace use_case
