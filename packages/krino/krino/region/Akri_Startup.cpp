// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Startup.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_RegisterProduct.hpp>

#include <mpi.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string>

#include <Kokkos_Core.hpp>
#include <stk_util/diag/WriterRegistry.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/ParseCommandLineArgs.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/environment/RuntimeMessage.hpp>
#include <stk_util/environment/Trace.hpp>
#include <stk_util/registry/ProductRegistry.hpp>

namespace krino{

namespace {

std::string
get_program_path(const char *program)
{
  // If we already have the full path, just return it
  if (program[0] == '/')
    return program;

  char full_path[PATH_MAX];
  if (strchr(program, '/') != nullptr) {
    realpath(program, full_path);
    return full_path;
  }

  char *PATH = getenv("PATH");
  while (PATH && *PATH) {
    // Get the character past the end of the next directory in PATH, i.e.
    // either the '/' or the '\0'
    char *end = strchr(PATH, ':');
    if (!end) {
      end = PATH+strlen(PATH);
    }

    // Set current = directory + '/' + program
    strncpy(full_path, PATH, end-PATH);
    full_path[end-PATH] = '/';
    strcpy(&full_path[end-PATH+1], program);

    // Check whether possible exists
    if (access(full_path, X_OK) == 0)
      return full_path;

    // Advance to the next directory
    PATH = *end ? end+1 : end;
  }

  // Not found; this shouldn't happen, but maybe the executable got deleted
  // after it was invoked before we got here -- or we have some crazy
  // parallel machine where the executable is inaccessible on the compute
  // nodes despite it somehow having been loaded.  No big deal, just return
  // the non-absolute path.
  return program;
}

} // namespace <unnamed>

Startup::Startup(int argc, char ** argv)
  : my_flag_exit_early(false)
{
  if ( MPI_SUCCESS != MPI_Init( &argc , &argv ) ) {
    throw std::runtime_error("MPI_Init failed");
  }

  stk::EnvData::instance().m_parallelComm = MPI_COMM_WORLD;
  MPI_Comm_size(stk::EnvData::parallel_comm(), &stk::EnvData::instance().m_parallelSize);
  MPI_Comm_rank(stk::EnvData::parallel_comm(), &stk::EnvData::instance().m_parallelRank);

  Kokkos::initialize(argc, argv);

  stk::register_message_type(stk::MSG_WARNING, 10000000, "Warning");
  stk::register_message_type(stk::MSG_DOOMED, 10000000, "Parser error");
  stk::register_message_type(stk::MSG_EXCEPTION, 1000000, "Exception");
  sierra::Diag::registerWriter("krinolog", krinolog, theDiagWriterParser());

  // Register this so that -version will know this region is used.
  register_product();

  stk::set_report_handler(report_handler);

  setup_commandline_options();

  // parse command line options
  const stk::OptionsSpecification &desc = stk::get_options_specification();
  stk::ParsedOptions &parsedOptions = stk::get_parsed_options();
  stk::parse_command_line_args(argc, const_cast<const char**>(argv), desc, parsedOptions);

  for (auto && writer : sierra::Diag::getWriterRegistry())
  {
    if (parsedOptions.count(writer.first))
    {
      writer.second.first->setPrintMask(writer.second.second->parse(parsedOptions[writer.first].as<std::string>().c_str()));
    }
  }

  if ( parsedOptions.count("help") ) {
    if (0 == stk::EnvData::parallel_rank())
    {
      std::cerr << desc << std::endl;
    }
    my_flag_exit_early = true;
  }

  if (parsedOptions.count("version")) {
    if (0 == stk::EnvData::parallel_rank())
    {
      std::cerr << "Version: Krino1.0" << std::endl;
    }
    my_flag_exit_early = true;
  }

  int local_flag_exit_early = my_flag_exit_early;
  MPI_Allreduce(&local_flag_exit_early, &my_flag_exit_early, 1, MPI_INT, MPI_MAX, stk::EnvData::parallel_comm());
  if (my_flag_exit_early) return;

  const std::string inputFileName = sierra::Env::get_param("input-deck");
  if (stk::EnvData::instance().m_inputFileRequired)
  {
    if(inputFileName == "")
    {
      throw std::runtime_error("No input file specified.  An input file must be specified with the '-i' option");
    }

    stk::EnvData::setInputFileName(inputFileName);
  }

  // setup logfile as inputFileName.log
  std::string logFileName = parsedOptions["output-log"].as<std::string>();
  if (logFileName == "")
  {
    int dotPos = inputFileName.rfind(".");
    if ( -1 == dotPos ) {
      // lacking extension
      logFileName = inputFileName + ".log";
    }
    else {
      // with extension; swap with .log
      logFileName = inputFileName.substr(0, dotPos) + ".log";
    }
  }

  // set up output streams
  const std::string output_description =
      (0 == stk::EnvData::parallel_rank()) ?
      "outfile=" + logFileName + " out>outfile+pout dout>out" :
      "out>pout dout>out";
  std::string parallel_output_description = " pout>null";
  if (parsedOptions.count("pout"))
  {
    std::ostringstream s;
    s << " poutfile=" << logFileName << "." << stk::EnvData::parallel_size() << "." << stk::EnvData::parallel_rank() << " pout>poutfile";
    parallel_output_description = s.str();
  }
  stk::bind_output_streams(output_description+parallel_output_description); // necessary for krinolog to work, otherwise you may get segfault
  stk::EnvData::instance().m_outputP0 = &sierra::out();

  // output run info
  const std::string program_path = get_program_path(argv[0]);
  const char * build_date_time = __DATE__ " " __TIME__;
  stk::EnvData::instance().m_executablePath = program_path;

  stk::EnvData::instance().m_productName = krino::get_product_name();
  stk::ProductRegistry::AttributeMap &product_attributes = stk::ProductRegistry::instance().getProductAttributeMap(sierra::Env::product_name());
  product_attributes[stk::ProductRegistry::BUILD_TIME] = build_date_time;
  product_attributes[stk::ProductRegistry::EXECUTABLE] = program_path;

  sierra::Env::outputP0()
     << sierra::Env::subsection_separator()
     << std::endl
     << "Executable Name   = " << sierra::Env::product_name() << std::endl
     << "Executable Version= " << stk::ProductRegistry::version() << std::endl
     << "Executable Date   = " << stk::ProductRegistry::instance().getProductAttribute(sierra::Env::product_name(), stk::ProductRegistry::BUILD_TIME) << std::endl
     << "Executable File   = " << sierra::Env::executable_file() << std::endl
     << "Run Start Date    = " << sierra::Env::startup_date() << std::endl << std::endl
     << "Working Directory = " << sierra::Env::working_directory() << std::endl
     << "Parsed Input File = " << sierra::Env::getInputFileName() << std::endl
     << sierra::Env::subsection_separator()
     << std::endl << std::endl ;
}

void Startup::setup_commandline_options()
{
  stk::OptionsSpecification desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("input-deck,i", "Analysis input file", stk::DefaultValue<std::string>(""))
    ("output-log,o", "Output log file path, one of : 'cout', 'cerr', or a file path", stk::DefaultValue<std::string>(""))
    ("aprepro,a", "Process (on) or don't process (off) input with aprepro. Default=on.", stk::ImplicitValue<std::string>("on"))
    ("define,D", "Define symbols for use in aprepro processing of input file", stk::ValueType<std::string>())
    ("pout", "use separate output log file for each MPI process");

  for (auto && writer : sierra::Diag::getWriterRegistry())
  {
    std::ostringstream str;
    str << "Diagnostic writer " << writer.first << std::endl;
    writer.second.second->describe(str);
    desc.add_options()(writer.first.c_str(), str.str().c_str(), stk::ValueType<std::string>());
  }

  stk::get_options_specification().add(desc);
}

Startup::~Startup()
{
  krino::MasterElementDeterminer::clear_master_elements();  // needed at least with MasterElementIntrepid to destruct Views on Basis objs before Kokkos::finalize()
  Kokkos::finalize();
  MPI_Finalize();
}

void Startup::handle_exception(const char * what, const bool is_parsing)
{
  krinolog << stk::diag::dendl;
  krinolog << "Exception: " << what << stk::diag::dendl;
  krinolog << stk::diag::Traceback::printTraceback(stk::diag::Traceback::snapshot()) << stk::diag::dendl;
  if (is_parsing)
  {
    if (0 == stk::EnvData::parallel_rank())
    {
      std::cerr << "*** Parser Error ***\n*** check log file for more information ***" << std::endl;
    }
    std::cout << std::flush;
    std::cerr << std::flush;

    ::sleep(1);                                 // Give the other processors a chance at
                                                // catching up, seems to help hanging problems.

    MPI_Finalize();
    std::exit(1);
  }

  if (0 == stk::EnvData::parallel_rank())
  {
    std::cerr << "*** ABORT on P0  ***\n*** check log file for more information ***" << std::endl;
  }
  else
  {
    // give root proc chance to die first
    ::sleep(2);
    std::cerr << "*** ABORT on P" <<stk::EnvData::parallel_rank() << "  ***\n"
        "*** may need to re-run with --pout option and check processor-specific log file for more information ***" << std::endl;
  }
  std::cout << std::flush;
  std::cerr << std::flush;

  const int errorCode{-1};
  MPI_Abort(stk::EnvData::parallel_comm(), errorCode);
}

void Startup::report_handler(const char * message, int type)
{
  /* %TRACE[SPEC]% */ Tracespec trace__("krino::Startup::report_handler( const char * message, int type)"); /* %TRACE% */

  if ((type & stk::MSG_TYPE_MASK) == stk::MSG_EXCEPTION)
  {
    sierra::pout() << stk::get_message_name(type) << " On Processor " << stk::EnvData::parallel_rank()  << ":" << std::endl
        << message << std::endl << std::endl;
    if (type & stk::MSG_SYMMETRIC)
    {
      if (stk::EnvData::parallel_rank() == 0)
      {
        std::cerr << stk::get_message_name(type) << " on all processors:" << std::endl
            << message << std::endl << std::endl;
      }
    }
    else
    {
      std::cerr << stk::get_message_name(type) << " on processor " << stk::EnvData::parallel_rank() << ":"
          << std::endl << message << std::endl << std::endl;
    }
  }

  if (type & stk::MSG_DEFERRED)
  {
    sierra::pout() << std::endl << "Deferred " << stk::get_message_name(type) << " On Processor " << stk::EnvData::parallel_rank()  << ":" << std::endl
      << message << std::endl << std::endl;
  }
  else
  {
    if (!(type & stk::MSG_SYMMETRIC))
    {
      sierra::pout() << "Ad Hoc " << stk::get_message_name(type) << " On Processor " << stk::EnvData::parallel_rank()  << ":" << std::endl
        << message << std::endl << std::endl;
    }
  }

  sierra::out() << stk::get_message_name(type) << ": " << message << std::flush;
}

} // namespace krino
