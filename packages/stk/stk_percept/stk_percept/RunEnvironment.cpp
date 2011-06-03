/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>      

#define USE_GETCWD 1
#if USE_GETCWD
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#endif

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/util/Bootstrap.hpp>
#include <stk_util/util/IndentStreambuf.hpp>

#include <stk_percept/RunEnvironment.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/OptionMask.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"

#define STK_PERCEPT_DEBUG_INPUT_ARGS 0

/// copied and edited from stk_util/use_cases/UseCaseEnvironment 

namespace {

  OptionMaskParser dw_option_mask("stk_percept diagnostic writer");
  OptionMaskParser timer_option_mask("stk_percept timers");

  //!stk::Bootstrap x(bootstrap);

} // namespace <empty>

namespace stk {
  namespace percept {

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
      static stk::diag::Timer s_timer = stk::diag::createRootTimer("root timer", timerSet());

      return s_timer;
    }

    std::string RunEnvironment::m_workingDirectory = "";

    RunEnvironment::RunEnvironment(
                                   int *         argc,
                                   char ***      argv, bool debug)
      : m_comm(stk::parallel_machine_init(argc, argv)),
        m_need_to_finalize(true), m_debug(debug), m_processCommandLine_invoked(false)
    {
      internal_initialize(argc, argv);
    }

    RunEnvironment::RunEnvironment(
                                   int *         argc,
                                   char ***      argv,
                                   stk::ParallelMachine comm, bool debug)
      : m_comm(comm),
        m_need_to_finalize(false), m_debug(debug), m_processCommandLine_invoked(false)
    {
      internal_initialize(argc, argv);
    }

    void RunEnvironment::internal_initialize(int* argc, char*** argv)
    {
      // Broadcast argc and argv to all processors.
      int parallel_rank = stk::parallel_machine_rank(m_comm);
      //int parallel_size = stk::parallel_machine_size(m_comm);

      if (m_debug && !parallel_rank)
        {
          for (int i = 0; i < *argc; ++i) {
            const std::string s((*argv)[i]);
            std::cout << "tmp 1 argv["<<i<<"]= " << s << std::endl;
          }
        }

      int pargc = *argc;
      std::string *argv_new = new std::string[pargc];
      int argc_new = 0;
      for (int i = 0; i < pargc; i++)
        {
          std::string pargvi( (*argv)[i] );
          int incr=0;
          if (pargvi == "-d")
            {
              pargvi = "--d="+std::string((*argv)[i+1]);
              incr=1;
            }
          if (pargvi == "-o")
            {
              pargvi = "--o="+std::string((*argv)[i+1]);
              incr=1;
            }
          if (incr)
            {
              i++;
            }
          argv_new[argc_new++] = pargvi;
        }

      *argc = argc_new;

      for (int i = 0; i < *argc; ++i) {
        (*argv)[i] = const_cast<char *>(argv_new[i].c_str());
        if (m_debug && !parallel_rank) std::cout << "modified argv["<<i<<"]= " << argv_new[i] << std::endl;
      }

      output_log_opt = "sierra.output.log";
      dw_opt = "";
      timer_opt = "";
      directory_opt = "";
      help_opt = 0;

      pout_opt = "-";
      dout_opt = "out";
      runtest_opt = "";

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
        if ( s == "-h" || s == "-help" || s == "--help") {
          //std::cout << "Found Help:: Usage: " << (*argv)[0] << " [options...]" << std::endl;
          printHelp();
          std::exit(0);
          return; 
        }
      }

      Util::setRank(parallel_rank);
      stk::BroadcastArg b_arg(m_comm, *argc, *argv);

      bootstrap();

      setSierraOpts(parallel_rank, *argc, *argv);
    }

    void RunEnvironment::processCommandLine(int* argc, char*** argv)
    {
      int parallel_rank = stk::parallel_machine_rank(m_comm);
      int parallel_size = stk::parallel_machine_size(m_comm);

      unsigned failed = 0;

      int success = processCLP(parallel_rank, *argc, *argv);
      failed = success == 0 ? 0u : 1u;

      if (failed)
        {
          if ( !parallel_rank)
            {
              std::cout << "Command Line error: echo of args:" << std::endl;
              for (int ii=0; ii < *argc; ii++)
                {
                  std::cout << "failed = 1, arg[" << ii  << "] = " << (*argv)[ii] << std::endl;
                }
              printHelp();
            }
          MPI_Barrier( m_comm );
          //stk::RuntimeDoomedSymmetric() << "parse_command_line";
          exit(1);
        }

      // Parse diagnostic messages to display
      dw().setPrintMask(dw_option_mask.parse(dw_opt.c_str()));

      // Parse timer metrics and classes to display
      stk::diag::setEnabledTimerMetricsMask(stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME);
      timerSet().setEnabledTimerMask(timer_option_mask.parse(timer_opt.c_str()));

      // Set working directory
      m_workingDirectory = "./";
      m_workingDirectory = directory_opt;

      if (m_workingDirectory.length() && m_workingDirectory[m_workingDirectory.length() - 1] != '/')
        m_workingDirectory += "/";

      std::string output_description = build_log_description( m_workingDirectory, parallel_rank, parallel_size);

      stk::bind_output_streams(output_description);

      dout() << "Output log binding: " << output_description << std::endl;

      // Start stk_percept root timer
      timer().start();
      //std::cout << "RunEnvironment::initialize done, m_workingDirectory= " << m_workingDirectory << std::endl;
      m_processCommandLine_invoked = true;

    }

    RunEnvironment::~RunEnvironment()
    {
      if (!m_processCommandLine_invoked)
        {
          throw std::runtime_error("RunEnvironment:: you must now invoke processCommandLine after constructing a RunEnvironment");
        }
      stk::report_deferred_messages(m_comm);

      // Stop stk_percept root timer
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

    void RunEnvironment::printHelp()
    {
      std::cout << "Usage: stk_percept_mesh_mod  [options...]" << std::endl;
      clp.printHelpMessage("RunEnvironment",std::cout);
    }

    void
    RunEnvironment::bootstrap()
    {
      /// \todo REFACTOR  Put these program options in a function
      ///                 that can be called without the bootstrapping.
      //* dw_option_mask.mask("search", stk::percept::LOG_SEARCH, "log search diagnostics");
      //* dw_option_mask.mask("transfer", stk::percept::LOG_TRANSFER, "log transfer diagnostics");
      //* dw_option_mask.mask("timer", stk::percept::LOG_TIMER, "log timer diagnostics");
      dw_option_mask.mask("all", stk::percept::LOG_ALWAYS, "log all");

      timer_option_mask.mask("mesh", stk::percept::TIMER_MESH, "mesh operations timers");
      //* timer_option_mask.mask("meshio", stk::percept::TIMER_MESH_IO, "mesh I/O timers");
      //* timer_option_mask.mask("transfer", stk::percept::TIMER_TRANSFER, "transfer timers");
      //* timer_option_mask.mask("search", stk::percept::TIMER_SEARCH, "search timers");


      //!stk::get_options_description().add(desc);

    }

    void RunEnvironment::
    setSierraOpts(int procRank, int argc, char* argv[])
    {
      Teuchos::oblackholestream blackhole;
      std::ostream &out = ( procRank == 0 ? std::cout : blackhole );

        if (m_debug)
          out << Teuchos::Teuchos_Version() << std::endl << std::endl;
    
        clp.setDocString("Run environment options" );
        clp.setOption("help",         &help_opt,        "help flag");
        clp.setOption("directory",    &directory_opt,   "working directory");
        clp.setOption("d",            &directory_opt,   "working directory");
        clp.setOption("output-log",   &output_log_opt,  "output log path");
        clp.setOption("o",            &output_log_opt,  "output log path");
        clp.setOption("pout",         &pout_opt,        "per-processor log file path");
        clp.setOption("dout",         &dout_opt,        "diagnostic output stream one of: 'cout', 'cerr', 'out' or a file path");
        clp.setOption("dw",           &dw_opt,          dw_option_mask.describe().c_str());
        clp.setOption("timer",        &timer_opt,       timer_option_mask.describe().c_str());
        clp.setOption("runtest",      &runtest_opt,     "runtest pid file");;
        //("exit", "do nothing and then exit");

    }

    int RunEnvironment::
    processCLP(int procRank, int argc, char* argv[])
    {
      Teuchos::oblackholestream blackhole;
      std::ostream &out = ( procRank == 0 ? std::cout : blackhole );

      bool success = true;
  
      try {

        /* There are also two methods that control the behavior of the
           command line processor.  First, for the command line processor to
           allow an unrecognized a command line option to be ignored (and
           only have a warning printed), use:
        */

        clp.recogniseAllOptions(true);
  
        /* Second, by default, if the parser finds a command line option it
           doesn't recognize or finds the --help option, it will throw an
           std::exception.  If you want prevent a command line processor from
           throwing an std::exception (which is important in this program since
           we don't have an try/catch around this) when it encounters a
           unrecognized option or help is printed, use:
        */
        clp.throwExceptions(true);

        /* We now parse the command line where argc and argv are passed to
           the parse method.  Note that since we have turned off std::exception
           throwing above we had better grab the return argument so that
           we can see what happened and act accordingly.
        */
        Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn= Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ;
        try {
          parseReturn = clp.parse( argc, argv );
        }
        catch (std::exception exc)
          {
            out << "RunEnvironment::processCLP error, exc= " << exc.what() << std::endl;
            return 1;
          }

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {

          return 0;
        }

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION   ) {

          if (m_debug)
            out << "RunEnvironment::processCLP error, unrecognized option" << std::endl;
          return 1; // Error!
        }

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
          if (m_debug)
            out << "RunEnvironment::processCLP success" << std::endl;
        }

        // Here is where you would use these command line arguments but for this example program
        // we will just print the help message with the new values of the command-line arguments.
        if (procRank == 0 && m_debug)
          {
            out << "\nPrinting help message with new values of command-line arguments ...\n\n";

            clp.throwExceptions(false);

            clp.printHelpMessage(argv[0],out);

            clp.throwExceptions(true);
          }

        // Now we will print the option values
        if (procRank == 0 && m_debug) {
          out << "\nPrinting user options after parsing ...\n\n";
          out << " output_log_opt= " <<  output_log_opt << std::endl;
          out << " dw_opt= " <<  dw_opt << std::endl;
          out << " timer_opt= " <<  timer_opt << std::endl;
          out << " directory_opt= " <<  directory_opt << std::endl;
          out << " help_opt= " <<  help_opt << std::endl;
        }

      } // try
      TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
  
      if(success && m_debug)
        out << "\nEnd Result: TEST PASSED" << std::endl;	

      return ( success ? 0 : 1 );
    }

    // Build output logging description for binding output streams
    std::string
    RunEnvironment::build_log_description(
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
      out_path = output_log_opt;
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
      if (true) {
        std::string pout_path = pout_opt;
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
      if (true) {
        std::string dout_path = dout_opt;
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

    // broken - do not use - for reference purposes only
    void RunEnvironment::
    doSierraLoadBalance(stk::ParallelMachine comm, std::string meshFileName)
    {
      unsigned p_size = stk::parallel_machine_size(comm);
      unsigned p_rank = stk::parallel_machine_rank(comm);
      if (!p_rank)
        {
          FILE *fpipe;
          std::string tmp_file = "stk_percept_tmp_loadbal.i";
          std::string command="sierra -i "+tmp_file;
          command += " -j "+toString(p_size);
          command += +" stk_percept_verifier_mesh -O\"--exit\"";
          //const char *command="sierra -i stk_percept_tmp_loadbal.i";
          char line[256];

          std::ofstream fout((std::string("./"+tmp_file)).c_str());
          fout << "Begin Sierra Percept\n"
               << "  Title None\n"
               << "  Begin Finite Element Model 3d_model\n"
               << "     Database Name = " << meshFileName << "\n"
               << "  End\n"
               << "End\n" << std::endl;
          fout.close();

          std::cout << "command= " << command << std::endl;
          if ( !(fpipe = (FILE*)popen(command.c_str(),"r")) )
            {  // If fpipe is NULL
              perror("Problems with pipe");
              exit(1);
            }
	 
          while ( fgets( line, sizeof line, fpipe))
            {
              //printf("%s", line);
              std::cout << line << std::endl;
            }
          pclose(fpipe);
        }
      MPI_Barrier( comm );
    }

#if USE_GETCWD
    static void panic()
    {
      (void)perror("Unexpected library error");
      std::abort();
    }

#define MAX_PATH 256
    //: Returns the name of the current working directory
    std::string RunEnvironment::
    get_working_directory()
    {
      int size = MAX_PATH;

      while(true) {
        char *ptr = new char[size];
        if (getcwd(ptr, size-1) != 0) {
          std::string ret_str(ptr);
          ret_str = ret_str +"/";
          delete [] ptr;
          return ret_str;
        }
        if (errno != ERANGE)
          panic();

        // Need larger char* to store current working directory name.
        delete [] ptr;
        size += MAX_PATH;
      }
    }
#else
    std::string RunEnvironment::get_working_directory()
    {
      return m_workingDirectory;
    }

#endif


    static void runCommand(std::string command)
    {
      char line[256];

      FILE *fpipe;
      if ( !(fpipe = (FILE*)popen(command.c_str(),"r")) )
        {  // If fpipe is NULL
          perror("Problems with pipe");
          exit(1);
        }
	 
      while ( fgets( line, sizeof line, fpipe))
        {
          //printf("%s", line);
          std::cout << line; // << std::endl;
        }
      pclose(fpipe);
    }

    void RunEnvironment::
    setFileNames(std::string& fullmesh, std::string& meshFileName)
    {
      fullmesh = get_working_directory() + meshFileName;
      if (meshFileName[0] == '/')
        {
          // already absolute
          fullmesh = meshFileName;
        }
      else if (meshFileName.length() >= 3 && (meshFileName[0] != '.' && meshFileName[1] != '/'))
        {
          // make it a relative path
          fullmesh = get_working_directory() + meshFileName;
          meshFileName = "./"+meshFileName;
          //fullmesh = get_working_directory() + meshFileName.substr(2,meshFileName.length()-2);
        }
      else if (meshFileName.length() >= 3 && (meshFileName[0] == '.' && meshFileName[1] == '/'))
        {
          fullmesh = get_working_directory() + meshFileName.substr(2,meshFileName.length()-2);
        }
      else
        {
          throw new std::runtime_error("bad format for input file name");
        }

      size_t found = 0;
      // strip off the basename + extension
      found = fullmesh.find_last_of("/");
      if (found != std::string::npos)
        {
          fullmesh = fullmesh.substr(0,found);
        }
      if(fullmesh == ".") 
        fullmesh = "./";
      else
        fullmesh += "/";


    }

    void RunEnvironment::
    doLoadBalance(stk::ParallelMachine comm, std::string meshFileName)
    {

      if (meshFileName.length() == 0) return;

      unsigned p_size = stk::parallel_machine_size(comm);
      unsigned p_rank = stk::parallel_machine_rank(comm);
      if (p_size > 1 && !p_rank)
        {
          std::string fullmesh = meshFileName;

          std::cout << "tmp get_working_directory= " << get_working_directory() << std::endl;

          setFileNames(fullmesh, meshFileName);

          std::cout << "tmp fullmesh= " << fullmesh << std::endl;
          std::cout << "tmp meshFileName= " << meshFileName << std::endl;

          //std::string fullmesh = meshFileName;

          //
          //  file.exo
          //  file.e
          unsigned extension_length = 1;

          size_t found = meshFileName.find_last_of(".");
          if (found == std::string::npos) {
            throw new std::runtime_error("RunEnvironment::doLoadBalance input file name must have an extension");
          }
          extension_length = meshFileName.length() - found - 1;
          std::string base_name = meshFileName.substr(0, meshFileName.length()-(extension_length+1));

          std::string extension = meshFileName.substr(meshFileName.length()-extension_length, meshFileName.length());

          std::cout << "tmp: extension= " << extension << " meshFileName= " << meshFileName << " fullmesh= " << fullmesh << " base_name= " << base_name << std::endl;

          std::string command="loadbal -No_subdirectory -spread -suffix_mesh " + extension+" -suffix_spread "+extension+" -p "+toString(p_size)+ " ";
          command += "-R " +fullmesh+" " + base_name;
          //command += "-R ./ "  + base_name;

          std::cout << "RunEnvironment::doLoadBalance: command= " << command << std::endl;

          runCommand(command);
        }

      MPI_Barrier( comm );
    }
  } // namespace percept
} // namespace stk
