/**   ------------------------------------------------------------
 *    Copyright 2001 - 2010 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <pwd.h>
#include <unistd.h>

#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <locale>
#include <map>

#include <stk_util/util/Null_Streambuf.hpp>
#include <stk_util/parallel/mpi_filebuf.hpp>

#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterRegistry.hpp>
#include <stk_util/diag/Env.hpp>
#include <stk_util/diag/Platform.hpp>
#include <stk_util/diag/Signal.hpp>
#include <stk_util/parallel/Exception.hpp>
#include <stk_util/parallel/ExceptionReport.hpp>
#include <stk_util/parallel/MPI.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ProductRegistry.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <stk_util/diag/UserPlugin.hpp>
#include <stk_util/parallel/mpih.hpp>
#include <stk_util/diag/PreParse.hpp>

#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/RuntimeMessage.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/Bootstrap.hpp>
#include <stk_util/util/IndentStreambuf.hpp>

namespace sierra {
namespace Env {

namespace {

void bootstrap()
{
  // Add my command line options to the option descriptions.
  boost::program_options::options_description desc("Runtime environment", 120);
  desc.add_options()
    ("help,h", "Display command line options")
    ("directory,d", boost::program_options::value<std::string>()->default_value("./"), "Set working directory")
    ("output-log,o", boost::program_options::value<std::string>()->default_value(""), "Output log file path, one of : 'cout', 'cerr', or a file path")
    ("logfile,l", boost::program_options::value<std::string>()->default_value(""), "Output log file path, one of : 'cout', 'cerr', or a file path")
    ("pout", boost::program_options::value<std::string>()->implicit_value("-"), "Per-processor log file path")
    ("dout", boost::program_options::value<std::string>()->implicit_value("out"), "Diagnostic output stream one of: 'cout', 'cerr', 'out' or a file path")
//    ("timer", boost::program_options::value<std::string>(), "Wall and CPU time options") // , &Diag::Timer::theTimerParser())
    ("version", "Display version information")
    ("jamsub", boost::program_options::value<std::string>(), "Display user subroutine build command")
    ("runtest", boost::program_options::value<std::string>()->implicit_value("pid"), "Record process host and pid to this file")
    ("developer-mode", "Activate developer specific features")
    ("architecture", boost::program_options::value<std::string>(), "Specifies the architecture running the sierra application");

  stk::get_options_description().add(desc);
}

stk::Bootstrap x(&bootstrap);

struct EnvData
{
  typedef std::map<ExecType, ExecInfo>    ExecMap;

  static EnvData &instance() {
    static EnvData s_env;

    return s_env;
  }

  EnvData()
    : m_productName("not specified"),
      m_vm(stk::get_variables_map()),
      m_nullBuf(),
      m_outputNull(&m_nullBuf),
      m_outputP0(&std::cout),
      m_output(),
      m_startTime((double) ::time(NULL)),
      m_executablePath(),
      m_shutdownRequested(false),
      m_inputFileRequired(true),
      m_checkSubCycle(false),
      m_worldComm(MPI_COMM_NULL),
      m_parallelComm(MPI_COMM_NULL),
      m_parallelSize(-1),
      m_parallelRank(-1),
      m_emptyString(),
      m_onString(PARAM_ON),
      m_inputFile("")
  {
    m_execMap[EXEC_TYPE_LAG].m_master      = -1;
    m_execMap[EXEC_TYPE_LAG].m_groupComm   = MPI_COMM_NULL;
    m_execMap[EXEC_TYPE_FLUID].m_master    = -1;
    m_execMap[EXEC_TYPE_FLUID].m_groupComm = MPI_COMM_NULL;
    stk::register_log_ostream(std::cout, "cout");
    stk::register_log_ostream(std::cerr, "cerr");
    
    stk::register_ostream(sierra::out(), "out");
    stk::register_ostream(sierra::pout(), "pout");
    stk::register_ostream(sierra::dout(), "dout");
    stk::register_ostream(sierra::tout(), "tout");
    
    static_cast<stk::indent_streambuf *>(sierra::dwout().rdbuf())->redirect(sierra::dout().rdbuf());
  }

  ~EnvData()
  {
    static_cast<stk::indent_streambuf *>(sierra::dwout().rdbuf())->redirect(std::cout.rdbuf());
  
    stk::unregister_ostream(tout());
    stk::unregister_ostream(dout());
    stk::unregister_ostream(pout());
    stk::unregister_ostream(out());

    stk::unregister_log_ostream(std::cerr);
    stk::unregister_log_ostream(std::cout);
  }
  
  std::string           m_productName;

  boost::program_options::variables_map & m_vm;
  
  null_streambuf	m_nullBuf;
  std::ostream		m_outputNull;
  std::ostream *	m_outputP0;
  std::ostringstream	m_output;

  double		m_startTime;
  std::string		m_executablePath;

  bool			m_shutdownRequested;
  bool                  m_inputFileRequired;
  bool                  m_checkSubCycle;

  MPI_Comm		m_worldComm;
  
  MPI_Comm		m_parallelComm;
  int			m_parallelSize;
  int			m_parallelRank;

  ExecMap               m_execMap;
  
  const std::string	m_emptyString;
  const std::string	m_onString;

  std::string           m_inputFile;
};

} // namespace <unnamed>

const std::string &
product_name()
{
  return EnvData::instance().m_productName;
}


const std::string &
executable_file()
{
  return EnvData::instance().m_executablePath;
}


const std::string &
executable_date()
{
  static std::string executable_date;

  if (executable_date.empty()) 
    executable_date = ProductRegistry::instance().getProductAttribute(EnvData::instance().m_productName, ProductRegistry::BUILD_TIME);

  return executable_date;
}


const std::string &
startup_date()
{
  static std::string startup_date;

  if (startup_date.empty())
    startup_date = format_time(EnvData::instance().m_startTime).c_str();

  return startup_date;
}


double
start_time()
{
  return EnvData::instance().m_startTime;
}


bool
developer_mode()
{
  return !get_param("developer-mode").empty();
}


void setInputFileName(std::string name) {
  EnvData::instance().m_inputFile = name;
}

std::string getInputFileName() {
  return EnvData::instance().m_inputFile;
}

void set_input_file_required(bool value)
{
    EnvData::instance().m_inputFileRequired = value;
}

void set_check_subcycle(bool value)
{
    EnvData::instance().m_checkSubCycle = value;
}


const std::string &
architecture()
{
  return get_param("architecture");
}


const std::string
working_directory() {
  char cwd[PATH_MAX];
  std::string directory = get_param("directory");
  if (directory[0] != '/' && getcwd(cwd, PATH_MAX) != NULL) {
    directory = cwd;
    directory += '/';
  }
  return directory;
}


std::ostream &
output()
{
  return EnvData::instance().m_output;
}


std::ostream &
outputP0()
{
  return *EnvData::instance().m_outputP0;
}


std::ostream &
outputNull() {
  return EnvData::instance().m_outputNull;
}


const char *
section_separator()
{
  static const char *s_sectionSeparator = "+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----";

  return s_sectionSeparator;
}


const char *
subsection_separator()
{
  static const char *s_subsectionSeparator = "---------------------------------------------------";

  return s_subsectionSeparator;
}


std::string
section_title(
  const std::string &	title)
{
  static size_t s_sectionSeparatorLength = std::strlen(section_separator());

  std::ostringstream strout;

  strout << std::left << std::setw(s_sectionSeparatorLength - 20) << title << std::right << std::setw(20) << format_time(Env::wall_now());
  return strout.str();
}


int parallel_size() {
  return EnvData::instance().m_parallelSize;
}

int parallel_rank() {
  return EnvData::instance().m_parallelRank;
}

MPI_Comm
parallel_comm()
{
  return EnvData::instance().m_parallelComm;
}

MPI_Comm
parallel_world_comm()
{
  return EnvData::instance().m_worldComm;
}

int parallel_lag_master() {
  return EnvData::instance().m_execMap[EXEC_TYPE_LAG].m_master;
}

int parallel_fluid_master() {
  return EnvData::instance().m_execMap[EXEC_TYPE_FLUID].m_master;
}

int peer_group() {
  return EnvData::instance().m_execMap[EXEC_TYPE_PEER].m_master;
}

std::string
get_program_path(const char *program)
{
  // If we already have the full path, just return it
  if (program[0] == '/')
    return program;

  char full_path[PATH_MAX];
  if (strchr(program, '/') != NULL) {
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

void parse_options(MPI_Comm comm, int *argc, char ***argv);
void startup_multi_exec(MPI_Comm world_comm, ExecType my_executable_type, const std::vector<int> *peer_sizes);



bool StartupSierra(int *			  argc,
  char ***		  argv,
  const char *		  product_name,
  const char *		  build_time,
  ExecType                mpi_key,
  const std::vector<int> *peer_sizes) {
  bool returnValue = false;
 
  stk::Bootstrap::bootstrap();
  
  EnvData &env_data = EnvData::instance();

  env_data.m_executablePath = get_program_path(*argv[0]);
  env_data.m_productName = product_name;

  ProductRegistry::instance().setProductName(product_name);

  ProductRegistry::AttributeMap &product_attributes = ProductRegistry::instance().getProductAttributeMap(product_name);
  product_attributes[ProductRegistry::BUILD_TIME] = build_time;
  product_attributes[ProductRegistry::EXECUTABLE] = env_data.m_executablePath;

  // Add Utility runtime library to the product registry
  sierra::register_product();

  // Add mpih to the product registry
  sierra::mpih::register_product();

  // Add operating system information to the product registry.
  ProductRegistry::AttributeMap &attr_map = ProductRegistry::instance().addProduct(osname().c_str());
  attr_map[ProductRegistry::VERSION]      = osversion().c_str();

  // Process the broadcast command line arguments
  namespace opt = boost::program_options;
    
  opt::variables_map &vm = stk::get_variables_map();
  opt::options_description &od = stk::get_options_description();
  {
    boost::program_options::options_description desc("Diagnostic writers", 120);
    
    for (Diag::WriterRegistry::iterator it = Diag::getWriterRegistry().begin(); it != Diag::getWriterRegistry().end(); ++it) {
      std::ostringstream str;
      str << "Diagnostic writer " << (*it).first << std::endl;
      (*it).second.second->describe(str);  
      desc.add_options()((*it).first.c_str(), boost::program_options::value<std::string>(), str.str().c_str());
    }
    
    std::ostringstream str;
    str << "Wall and CPU time options" << std::endl;
    Diag::theTimerParser().describe(str);  
    desc.add_options()("timer", boost::program_options::value<std::string>(), str.str().c_str());
    
    od.add(desc);
  }

  for (int i = 0; i < *argc; ++i) {
    const std::string s((*argv)[i]);
    if (s == "-h" || s == "-help" || s == "--help") {
      std::cout << std::endl
                << "Sierra Usage: sierra " << lower(product_name) << " [sierra-options...] -O \"[" << lower(product_name) << "-options...]\"" << std::endl << std::endl
//                << "Usage: (MPI run) " << env_data.m_executablePath << " [options...]" << std::endl
                << "For example:" << std::endl
                << "" << std::endl
                << "  sierra " << lower(product_name) << " -i input_deck.i -o sierra.log" << std::endl
                << "    This creates the normal output file sierra.log" << std::endl
                << "" << std::endl
                << "  sierra " << lower(product_name) << " -i input_deck.i -o sierra.log -O \"--pout=pp.log\"" << std::endl
                << "    The per-processor output is written to pp.log.n.r for each rank, r, of n processors." << std::endl
                << "" << std::endl
                << "  sierra " << lower(product_name) << " -i input_deck.i -o sierra.log -O \"--fmwkout=field,parameters\"" << std::endl
                << "    Enable the framework field and parameter diagnostics" << std::endl
                << "" << std::endl
                << "  sierra " << lower(product_name) << " -i input_deck.i -o sierra.log -O \"--timer=all\"" << std::endl
                << "    Enable the all timers" << std::endl
                << std::endl
                << "  For additional information see:" << std::endl
                << "      http://sierra-dev.sandia.gov/stk/group__stk__util__output__log__detail.html#stk_util_output_log_howto_use_in_sierra_app" << std::endl << std::endl
                << product_name << " options are:" << std::endl
                << stk::get_options_description() << std::endl;
      std::exit(0);
    }
  }

  for (int i = 0; i < *argc; ++i) {
    const std::string s((*argv)[i]);
    if (s == "-jamsub" || s == "--jamsub") {
      const char *t = (*argv)[i + 1];
      const char **symbol = sierra::Plugin::Registry::getsym<const char **>(t);
      if (symbol) {
        std::cout << *symbol << std::endl;
        std::exit(0);
      }
      else
        std::exit(1);
    }
  }

  try {
    startup_preparallel_platform();

    // Communicator has not been set, initialize MPI if not already initialized
    int mpi_init_val = 0 ;
    if ( MPI_SUCCESS != MPI_Initialized( &mpi_init_val ) ) {
      throw RuntimeError() << "MPI_Initialized failed";
    }

    // Default startup communicator
    MPI_Comm startup_mpi_comm = MPI_COMM_WORLD;

    // If we are initializing the comm, see if there are differing
    // executables running.  If there are, find our partition and the
    // leads of the other partitions.
    if ( mpi_init_val == 0 ) {
      if ( MPI_SUCCESS != MPI_Init( argc , argv ) ) {
	throw RuntimeError() << "MPI_Init failed";
      }

      returnValue = true ;

      if (mpi_key != EXEC_TYPE_WORLD) startup_multi_exec(startup_mpi_comm, mpi_key, peer_sizes);
    }
    
    // Ready to reset the environment from NULL, we are the Lagrangian application at this point.
    MPI_Comm new_comm = mpi_key != EXEC_TYPE_WORLD ? env_data.m_execMap[mpi_key].m_groupComm : MPI_COMM_WORLD;
    reset(new_comm);
  }
  catch (const std::exception &x) {
    std::cerr << "SIERRA execution failed during mpi initialization with the following exception:" << std::endl
	      << x.what() << std::endl;
    MPI_Abort(env_data.m_parallelComm , MPI_ERR_OTHER);
  }
  catch (...) {
    std::cerr << "SIERRA execution failed during mpi initialization  with unknown exception:" << std::endl;

    MPI_Abort(env_data.m_parallelComm, MPI_ERR_OTHER);
  }

  parse_options(env_data.m_parallelComm, argc, argv);
  
  {
    std::ostringstream output_description;

    // On processor 0:
    //   [outfile=path] [poutfile=path.n.r] [doutfile=path.n.r] out>{-|cout|cerr|outfile}+pout pout>{null|poutfile} dout>{out|doutfile}

    // On processor 1..n:
    //   [poutfile=path.n.r] [doutfile=path.n.r] out>pout pout>{null|poutfile} dout>{out|doutfile}
    
    std::string out_path1 = vm["output-log"].as<std::string>();
    std::string out_path2 = vm["logfile"].as<std::string>();


    
    std::string originalFileName = Env::get_param("input-deck");
    std::string modifiedFileName = originalFileName;

    if(originalFileName == "") {
      //
      //  If no input file specified, error out (unless just running the --version or --help option)
      //         
      if ( get_param("version").empty() && get_param("help").empty() ) {
        if (env_data.m_inputFileRequired) {
          throw RuntimeError() << "No input file specified.  An input file must be specified with the '-i' option";
        } else {
          std::cerr << "WARNING: No input file specified.  An input file should be specified with the '-i' option!" << std::endl;
        }
      }
    } else if ( env_data.m_checkSubCycle ) {
      // Alter input-deck if subcycle present
      bool debugSubCycleSplit = false;
      std::string subCycleRegexp("^\\s*subcycle\\s+blocks\\s*=");
      bool subCycleSet = CaseInSensitiveRegexInFile(subCycleRegexp, originalFileName, debugSubCycleSplit);
      std::string coarseRegionRegexp("^\\s*begin\\s+presto\\s+region\\s+\\w+_AutoCoarseRegion\\>");
      bool coarseRegionMade = CaseInSensitiveRegexInFile( coarseRegionRegexp, originalFileName, debugSubCycleSplit);
      std::string fineRegionRegexp("^\\s*begin\\s+presto\\s+region\\s+\\w+_AutoFineRegion\\>");
      bool fineRegionMade = CaseInSensitiveRegexInFile( fineRegionRegexp, originalFileName, debugSubCycleSplit);
      if ( subCycleSet ) {
        if ( !coarseRegionMade && !fineRegionMade ) {
          modifiedFileName = CreateSubCycleInputFile( originalFileName );
        } else {
          if(Env::parallel_rank() == 0) {
	    std::cout<<"Input File: " << originalFileName << " Appears to have already been converted for subcycling.  ";
	    std::cout<<"Skipping input conversion " << std::endl;
	  }
        }
      }
    }

    setInputFileName(modifiedFileName);


    std::string trueOut;
    if(out_path2 != "") {
      trueOut = out_path2;
    } else if(out_path1 != "") {
      //
      //  Old syntax compatibility, access the old output-file executable option if the logfile is not defined
      //
      trueOut = out_path1;
    } else {
     //
      //  If log file name is unspecified, default it to (Base Input File Name).log
      //  Use the following logic:
      //   If the input file has an extension, replace the last ".extension" with ".log"
      //   If the input file has no extension, append ".log" to the input file name
      //   If the input file contains the word '.aprepro', assume aprepro was used to convert and strip out the aprepro
      //   If the input file contains any directory movement (like ../) strip them out so log file is written to current director
      //

      int dotPos = originalFileName.rfind(".");

      if(dotPos == -1) {  //No extension
        trueOut = originalFileName + ".log";
      } else {  //Extension found
        trueOut = originalFileName.substr(0, dotPos) + ".log";
      }
      //
      //  If the output path contains a ".aprepro" tag get rid of it
      //
      int apreproPos = trueOut.rfind(".aprepro"); 
      if(apreproPos != -1) {
        trueOut.erase(apreproPos, 8);
      }
      //
      //  If the output path contains a "aaa/input.i" pull off the initial directory redirects so that the log file is written int the current directory
      //
      int lastSlashPos = trueOut.rfind("/");

      if(lastSlashPos != -1) {
        trueOut.erase(0,lastSlashPos+1);
      }


    }
 
    std::string out_path = trueOut;

    if (out_path == "-")
      out_path = "cout";
    
    std::string out_ostream;

    if (!stk::get_log_ostream(out_path))
      if (out_path.size() && out_path[0] != '/')
        out_path = working_directory() + out_path;

    if (parallel_rank() == 0) {
      if (!stk::get_log_ostream(out_path)) {
        output_description << "outfile=\"" << out_path << "\"";
        out_ostream = "outfile";
      }
      else {
        out_ostream = out_path;
      }
    }
    else
      out_ostream = "null";

    std::string pout_ostream = "null";
    if (vm.count("pout")) {
      std::string pout_path = vm["pout"].as<std::string>();
      if (pout_path == "-") {
        std::ostringstream s;

        if (stk::get_log_ostream(out_path))
          s << working_directory() << "sierra.log." << parallel_size() << "." << parallel_rank();
        else
          s << out_path << "." << parallel_size() << "." << parallel_rank();
        pout_path = s.str();
      }
      else if (pout_path.find("/") == std::string::npos && !stk::get_log_ostream(pout_path)) {
        std::ostringstream s;

        s << working_directory() << pout_path << "." << parallel_size() << "." << parallel_rank();
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
          s << working_directory() << dout_path << "." << parallel_size() << "." << parallel_rank();
        else
          s << dout_path << parallel_size() << "." << parallel_rank();;
        dout_path = s.str();
        output_description << " doutfile=\"" << dout_path << "\"";
        dout_ostream = "doutfile";
      }
    }
    else
      dout_ostream = "out";

    if (parallel_rank() == 0)
      output_description << " out>" << out_ostream << "+pout";
    else
      output_description << " out>pout";

    output_description << " pout>" << pout_ostream << " dout>" << dout_ostream;


    stk::bind_output_streams(output_description.str());
  }
  
  env_data.m_outputP0 = &sierra::out();
  
#ifdef SIERRA_EXPORT_CONTROL_EAR99
  // If you are using an EAR99 export controlled version of Sierra,
  // any attempt to modify or bypass this section of code is a 
  // violation of U.S. Export Control Regulations and subject to
  // criminal prosecution.
  if (parallel_size() > SIERRA_EXPORT_CONTROL_EAR99) {
    if (parallel_rank() == 0) {
      std::cerr << "ERROR: You are running an EAR99 export controlled version of\n";
      std::cerr << "       Sierra. For this export control level, a maximum of\n";
      std::cerr << "       "<<SIERRA_EXPORT_CONTROL_EAR99<<" processors is permitted\n";
    }
    MPI_Abort(env_data.m_parallelComm, MPI_ERR_OTHER);
  }
#endif  

  try {
    // Create pid file if runtest command line option specified
    if ( !get_param("runtest").empty() ) {

      mpi_filebuf mpi_buf;

      mpi_buf.open(env_data.m_parallelComm, 0, std::ios::out, get_param("runtest").c_str());

      if ( ! mpi_buf.is_open() )
	throw RuntimeError() << "failed to open pid file " << get_param("runtest");

      std::ostream s( &mpi_buf );
      s << parallel_rank() << ":" << hostname() << domainname() << ":" << pid() << ":" << pgrp() << std::endl;
    }

    // Enable the timers
    if (!get_param("timer").empty()) {
      Diag::TimerParser	parser;

      Diag::sierraTimerSet().setEnabledTimerMask(parser.parse(get_param("timer").c_str()));
    }

    // Enable parallel exception handling, waited until now because it needs the Env output streams
    register_stl_parallel_exceptions();
  }
  catch (const std::exception &x) {
    std::cerr << "SIERRA execution failed during diagnostic and timer initialization with the following exception:" << std::endl
	      << x.what() << std::endl;
    abort();
  }
  catch (...) {
    std::cerr << "SIERRA execution failed during diagnostic and timer initialization with unknown exception:" << std::endl;
    abort();
  }

// Setup the hangup, segmentation violation, illegal instruction, bus error and
//    terminate signal handlers.
  if (get_param("nosignal").empty())
    activate_signals();


  return returnValue;
}




void
Startup::startup(
  int *			  argc,
  char ***		  argv,
  const char *		  product_name,
  const char *		  build_time,
  ExecType                mpi_key,
  const std::vector<int> *peer_sizes) {
  m_mpiInitFlag = StartupSierra(argc, argv, product_name, build_time, mpi_key, peer_sizes);
}

          
Startup::Startup(
  int *                 argc,
  char ***              argv,
  const char *          product_name,
  const char *          build_date_time,
  ExecType              mpi_key,
  const std::vector<int> *peer_sizes)
  : m_mpiInitFlag(false)
{
  startup(argc, argv, product_name, build_date_time, mpi_key, peer_sizes);
}


void ShutDownSierra(bool mpiInitFlag) {
  if (get_param("nosignal").empty())
    deactivate_signals();

  mpih::Delete_Handles();

  EnvData &env_data = EnvData::instance();
  mpih::Keyval_delete(env_data.m_parallelComm);

  reset(MPI_COMM_NULL);

  if (mpiInitFlag)
    MPI_Finalize();
}



Startup::~Startup() {
  ShutDownSierra(m_mpiInitFlag);
}


void parse_options(MPI_Comm  comm,
		   int *     argc,
		   char ***  argv)
{
  try {
    char ** argv2 = new char *[*argc];
    for (int i = 0; i < *argc; ++i) {
      if (std::strlen((*argv)[i]) > 2 && (*argv)[i][0] == '-' && (*argv)[i][1] != '-') {
        argv2[i] = new char[std::strlen((*argv)[i]) + 2];
        argv2[i][0] = '-';
        std::strcpy(&argv2[i][1], (*argv)[i]);
      }
      else {	
        argv2[i] = new char[std::strlen((*argv)[i]) + 1]; 
	std::strcpy(argv2[i], (*argv)[i]);
      }      
    }
  
    // Broadcast argc and argv to all processors.
    stk::BroadcastArg b_arg(comm, *argc, argv2);

    for (int i = 0; i < *argc; ++i)
      delete[] argv2[i];
    delete[] argv2;

    namespace opt = boost::program_options;
    opt::variables_map &vm = stk::get_variables_map();
    opt::options_description &od = stk::get_options_description();
    opt::store(opt::parse_command_line(b_arg.m_argc, b_arg.m_argv, od, opt::command_line_style::unix_style), vm);
    opt::notify(vm);

    for (Diag::WriterRegistry::iterator it = Diag::getWriterRegistry().begin(); it != Diag::getWriterRegistry().end(); ++it)
      if (vm.count((*it).first.c_str()))
        (*it).second.second->parse(vm[(*it).first.c_str()].as<std::string>().c_str());
    

    // Must have a working directory
    const std::string &working_dir = get_param("directory");
    if ( working_dir.empty() || working_dir == PARAM_ON )
      throw RuntimeError() << "working directory must be specified";
    if (working_dir[working_dir.length() - 1] != '/')
      const_cast<std::string &>(working_dir) += '/';
    
  }
  catch (const std::exception &x) {
    std::cerr << "SIERRA execution failed during command line processing with the following exception:" << std::endl
	      << x.what() << std::endl;
    MPI_Abort(comm, MPI_ERR_OTHER);
  }
  catch (...) {
    std::cerr << "SIERRA execution failed during command line processing with unknown exception:" << std::endl;

    MPI_Abort(comm, MPI_ERR_OTHER);
  }
} 

void
startup_multi_exec(MPI_Comm                world_comm,
		   ExecType                my_executable_type,
                   const std::vector<int> *peer_sizes)  // can be NULL.
{
  EnvData &env_data = EnvData::instance();

  // MPI interface construction
  int world_size = -1 ;
  int world_rank = -1 ;
      
  if ( MPI_Comm_size(world_comm, &world_size) != MPI_SUCCESS)
    throw RuntimeError() << "MPI_Comm_size failed";

  if ( MPI_Comm_rank(world_comm, &world_rank) != MPI_SUCCESS || -1 == world_rank )
    throw RuntimeError() << "MPI_Comm_rank failed";

  if (my_executable_type == EXEC_TYPE_FLUID || my_executable_type == EXEC_TYPE_LAG) {
    // This is specific for gemini.  Gemini performs three broadcasts, one for the
    // EXEC_TYPE_FLUID and one for the EXEC_TYPE_LAG.  Also note that the ranks of processors must
    // be ordered such that all gemini processors come first.  Gemini mandates that it;s master is
    // processor 0 and use ranks through its size.
    int lag_master = 0;
    int lag_rank_size = -1;
    int fluid_master = 0;
      
    if (world_rank == 0) {
      typedef std::map<ExecType, std::vector<int> > ExecTypeRanks;

      ExecTypeRanks exec_type_ranks;

      exec_type_ranks[my_executable_type].push_back(0);

      for (int i = 1; i < world_size; ++i) {
        MPI_Status status;
        int proc_stat[2];         // rank, ExecType
        if (MPI_Recv(proc_stat, 2, MPI_INTEGER, i, MPI_ANY_TAG, world_comm, &status) != MPI_SUCCESS)
          throw RuntimeError() << "MPI_Recv failed";

        exec_type_ranks[(ExecType) proc_stat[1]].push_back(proc_stat[0]);
      }        

      std::vector<int> &fluid_ranks = exec_type_ranks[EXEC_TYPE_FLUID];
      if (fluid_ranks.size())
        fluid_master = fluid_ranks.front();

      if (MPI_Bcast(&fluid_master, 1, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Bcast failed";

      std::vector<int> &lag_ranks = exec_type_ranks[EXEC_TYPE_LAG];      
      if (lag_ranks.size())
        lag_master = lag_ranks.front();

      if (MPI_Bcast(&lag_master, 1, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Bcast failed";

      lag_rank_size = lag_ranks.size();
      if (MPI_Bcast(&lag_rank_size, 1, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Bcast failed";
    }
    else {
      int proc_stat[2];
      proc_stat[0] = world_rank;
      proc_stat[1] = my_executable_type;

      if (MPI_Send(proc_stat, 2, MPI_INTEGER, 0, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Send failed";

      if (MPI_Bcast(&fluid_master, 1, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Bcast failed";

      if (MPI_Bcast(&lag_master, 1, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Bcast failed";

      if (MPI_Bcast(&lag_rank_size, 1, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Bcast failed";
    }

    MPI_Comm lag_comm   = world_comm;
    MPI_Comm fluid_comm = MPI_COMM_NULL;
    const int fluid_rank_size = world_size - lag_rank_size;
    if (fluid_rank_size) {

      MPI_Group world_group;
      MPI_Group lag_group;
      MPI_Group fluid_group;

      if (MPI_Comm_group(world_comm, &world_group) != MPI_SUCCESS) 
        throw RuntimeError() << "MPI_Comm_group failed";

      std::vector<int> lag_ranks;
      for (int i = 0; i < lag_rank_size; ++i)
        lag_ranks.push_back(lag_master + i);

      if (MPI_Group_incl(world_group, lag_ranks.size(), &lag_ranks[0], &lag_group) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Group_incl failed";
      if (MPI_Comm_create(world_comm, lag_group, &lag_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Comm_create failed";

      std::vector<int> fluid_ranks;
      for (int i = 0; i < fluid_rank_size; ++i)
        fluid_ranks.push_back(fluid_master + i);

      if (MPI_Group_incl(world_group, fluid_ranks.size(), &fluid_ranks[0], &fluid_group) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Group_incl failed";
      if (MPI_Comm_create(world_comm, fluid_group, &fluid_comm) != MPI_SUCCESS)
        throw RuntimeError() << "MPI_Comm_create failed";
    }

    env_data.m_worldComm                            = world_comm;
    env_data.m_execMap[EXEC_TYPE_LAG].m_master      = lag_master;
    env_data.m_execMap[EXEC_TYPE_LAG].m_groupComm   = lag_comm;
    env_data.m_execMap[EXEC_TYPE_FLUID].m_master    = fluid_master;
    env_data.m_execMap[EXEC_TYPE_FLUID].m_groupComm = fluid_comm;
  }
  else if (my_executable_type == EXEC_TYPE_PEER) {
    // This executable will run on 2 or more communicators.

    // NOTE: Only 2 communicators is currently supported...

    // If peer_sizes is NULL, then split world_comm into two equal
    // size communicators (peer(1) is larger if world_comm size is
    // odd)
    // If peer_sizes is not NULL, then split world_comm into
    // peer_sizes.size() sub communicators with peer(i) of size
    // peer_sizes(i). 

    // Sync 'peer_sizes' across all processors if non-null
    // For now, we limit the number of peer applications to 2.

    if (peer_sizes != NULL && peer_sizes->size() > 2) {
      throw RuntimeError() << "The total number of peer application processor sizes specfied is "
			   << peer_sizes->size()
			   << ",  but the current limit is 2.";
    }

    // Peer sizes is only set correctly on processor 0 since it was passed in by the
    // main routine prior to MPI_Init being called.  Broadcast the values to all processors.
    int peers[2];
    if (world_rank == 0) {
      if (peer_sizes != NULL) {
	peers[0] = (*peer_sizes)[0];
	peers[1] = (*peer_sizes)[1];
      } else {
	peers[0] = world_size / 2;
	peers[1] = world_size - world_size/2;
      }
    }
    if (MPI_Bcast(peers, 2, MPI_INTEGER, 0, world_comm) != MPI_SUCCESS)
      throw RuntimeError() << "MPI_Broadcast -- peers failed";

    // Check that the number of processes specified is equal to the
    // total number of processes
    int peer_proc_count = peers[0] + peers[1];
    if (peer_proc_count != world_size) {
      throw RuntimeError() << "The total number of peer processors specfied is " << peer_proc_count
			   << " which is not equal to the total number of processors (" << world_size << ").";
    }

    int my_peer_group = MPI_UNDEFINED;
    int sum = 0;
    for (size_t i=0; i < 2; i++) {
      sum += peers[i];
      if (world_rank < sum) {
	my_peer_group = i;
	break;
      }
    }

    MPI_Comm peer_comm;
    if (MPI_Comm_split(world_comm, my_peer_group, world_rank, &peer_comm) != MPI_SUCCESS) {
      throw RuntimeError() << "MPI_Comm_split failed";
    }
    env_data.m_worldComm                           = world_comm;
    env_data.m_execMap[EXEC_TYPE_PEER].m_groupComm = peer_comm;
    env_data.m_execMap[EXEC_TYPE_PEER].m_master    = my_peer_group; // Overloading meaning to peer group.
  }
}

bool
is_comm_valid()
{
  EnvData &env_data = EnvData::instance();
  if (env_data.m_parallelComm == MPI_COMM_NULL) {
    return false;
  } else {
    return true;
  }
}

void
reset(
  MPI_Comm		new_comm)
{
  EnvData &env_data = EnvData::instance();

  // Destroy old comm
  if (env_data.m_parallelComm != MPI_COMM_NULL) {

    if (new_comm != MPI_COMM_NULL) {
      mpih::Sub_Communicator(env_data.m_parallelComm, new_comm);
    }

    env_data.m_parallelComm = MPI_COMM_NULL ;
    env_data.m_parallelSize = -1;
    env_data.m_parallelRank = -1 ;
  }

  setMpiCommunicator(new_comm);
}

void setMpiCommunicator(MPI_Comm communicator)
{
    EnvData &env_data = EnvData::instance();
    if(communicator != MPI_COMM_NULL)
    {
        env_data.m_parallelComm = communicator;

        if(MPI_Comm_size(env_data.m_parallelComm, &env_data.m_parallelSize) != MPI_SUCCESS
           || MPI_Comm_rank(env_data.m_parallelComm, &env_data.m_parallelRank) != MPI_SUCCESS
           || env_data.m_parallelSize == -1
           || env_data.m_parallelRank == -1)
        {
            throw RuntimeError() << "reset given bad MPI communicator";
        }
    }
}

void
output_flush()
{
  EnvData &env_data = EnvData::instance();

  stk::report_deferred_messages(Env::parallel_comm());
  
  stk::all_write_string(Env::parallel_comm(), *env_data.m_outputP0, env_data.m_output.str());
  env_data.m_output.str("");
}


void
request_shutdown(bool shutdown)
{
  EnvData::instance().m_shutdownRequested = shutdown;
}


bool
is_shutdown_requested()
{
  int shutdown_requested_in = EnvData::instance().m_shutdownRequested || Env::HUP_received();
  int shutdown_requested;

  MPI_Allreduce(&shutdown_requested_in, &shutdown_requested, 1, MPI_INT, MPI_SUM, Env::parallel_comm());

  return shutdown_requested != 0;
}


void abort() {
  EnvData &env_data = EnvData::instance();

  // Cannot be sure of parallel synchronization status; therefore, no communications can
  // occur.  Grab and dump all pending output buffers to 'std::cerr'.
  std::cerr << std::endl
            << "*** SIERRA ABORT on P" << EnvData::instance().m_parallelRank << " ***"
            << std::endl
            << "*** check " << get_param("output-log")
            << " file for more information ***"
            << std::endl ;

  if (!env_data.m_output.str().empty()) {
    std::cerr << "Buffer contents of deferred output stream on processor " << parallel_rank()
              << std::endl ;
    std::cerr << env_data.m_output.str();
  }
  
  std::cerr.flush();
  std::cout.flush();

  ::sleep(1);					// Give the other processors a chance at
						// catching up, seems to help hanging problems.
  MPI_Abort(env_data.m_parallelComm, MPI_ERR_OTHER);	// First try to die
  std::exit( EXIT_FAILURE );                    // Second try to die
}


const std::string &
get_param(
  const char * const	option)
{
  if (EnvData::instance().m_vm.count(option)) {
    if (EnvData::instance().m_vm[option].as<std::string>().empty())
      return EnvData::instance().m_onString;
    else
      return EnvData::instance().m_vm[option].as<std::string>();
  }
  else
    return EnvData::instance().m_emptyString;
}


void
set_param(
  const char *          option,
  const std::string &   value) {


  namespace opt = boost::program_options;

  opt::variables_map &vm = stk::get_variables_map();
  opt::options_description &od = stk::get_options_description();

  int argc = 1;
  char *s = std::strcpy(new char[std::strlen(option) + 1], option);
  
  opt::store(opt::parse_command_line(argc, &s, od), vm);
  opt::notify(vm);

  delete [] s;
}

} // namespace Env
} // namespace sierra
