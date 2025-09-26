// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/util/Writer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/util/Bootstrap.hpp>
#include <stk_util/util/IndentStreambuf.hpp>
#include <stk_util/environment/Env.hpp>

#include <percept/RunEnvironment.hpp>
#include <percept/Util.hpp>
#include <percept/OptionMask.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"

#define STK_PERCEPT_DEBUG_INPUT_ARGS 0

/// copied and edited from stk_util/use_cases/UseCaseEnvironment

namespace {

  OptionMaskParser dw_option_mask("Percept diagnostic writer");
  OptionMaskParser timer_option_mask("Percept timers");

  //!stk::Bootstrap x(bootstrap);

} // namespace <empty>

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
      switch (static_cast<int>(type) & static_cast<int>(stk::MSG_TYPE_MASK)) {
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
                                   int  *        argc,
                                   char ***      argv, bool debug)
      : ParallelMachineFinalize(false),
        m_comm( stk::parallel_machine_init(argc, argv)),
        m_need_to_finalize(true), m_debug(debug), m_processCommandLine_invoked(false), m_argv_new(0),m_argc(0),m_argv(0)
        //,m_par_finalize(false)
    {
      //internal_initialize(*argc, *argv);
    }

    RunEnvironment::RunEnvironment(
                                   int   *      /*argc*/,
                                   char ***      /*argv*/,
                                   stk::ParallelMachine comm, bool debug)
      : ParallelMachineFinalize(false),
        m_comm(comm),
        m_need_to_finalize(false), m_debug(debug), m_processCommandLine_invoked(false), m_argv_new(0),m_argc(0),m_argv(0)
        //,m_par_finalize(false)

    {
      //internal_initialize(*argc, *argv);
    }

    void RunEnvironment::internal_initialize(int argc, char** argv)
    {
      // Broadcast argc and argv to all processors.
      int parallel_rank = stk::parallel_machine_rank(m_comm);
      //int parallel_size = stk::parallel_machine_size(m_comm);

      if (m_debug && !parallel_rank)
        {
          for (int i = 0; i < argc; ++i) {
            const std::string s((argv)[i]);
            std::cout << "tmp 1 argv["<<i<<"]= " << s << std::endl;
          }
        }

      int pargc = argc;
      m_argv_new = new std::string[pargc];
      int argc_new = 0;
      for (int i = 0; i < pargc; i++)
        {
          std::string pargvi( (argv)[i] );
          int incr=0;
          if (pargvi == "-d")
            {
              pargvi = "--d="+std::string((argv)[i+1]);
              incr=1;
            }
          if (pargvi == "-o")
            {
              pargvi = "--o="+std::string((argv)[i+1]);
              incr=1;
            }
          if (incr)
            {
              i++;
            }
          m_argv_new[argc_new++] = pargvi;
        }

      m_argc = argc_new;
      m_argv = new char*[m_argc];

      for (int i = 0; i < m_argc; ++i) {
        m_argv[i] = new char[m_argv_new[i].length()+1];
        strcpy(m_argv[i], m_argv_new[i].c_str());
        if (m_debug && !parallel_rank) std::cout << "modified argv["<<i<<"]= " << m_argv_new[i] << std::endl;
      }

      output_log_opt = "percept.output.log";
      dw_opt = "";
      timer_opt = "";
      directory_opt = "";
      help_opt = 0;

      pout_opt = "";
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

      Util::setRank(parallel_rank);
      stk::BroadcastArg b_arg(m_comm, argc, argv);

      bootstrap();

      setSierraOpts(parallel_rank, argc, argv);

    }

    int processCommandLine(Teuchos::CommandLineProcessor& clp_in, int argc, char** argv, int &bad_option)
    {
      int parallel_rank = stk::parallel_machine_rank(sierra::Env::parallel_comm());

      int error = processCLP(clp_in, parallel_rank, argc, argv, bad_option);

      // Start percept root timer
      timer().start();
      //m_processCommandLine_invoked = true;

      return error;
    }

    RunEnvironment::~RunEnvironment()
    {
      if (!m_processCommandLine_invoked)
        {
          std::runtime_error("RunEnvironment:: you must now invoke processCommandLine after constructing a RunEnvironment");
        }
      stk::report_deferred_messages(m_comm);

      // Stop percept root timer
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

      if (m_argv_new) delete[] m_argv_new;
      if (m_argv)
        {
          for (int i=0; i < m_argc; ++i)
            {
              delete[] m_argv[i];
            }
          delete[] m_argv;
        }

      if (m_need_to_finalize) {
        //stk::parallel_machine_finalize();
      }
    }

    void
    RunEnvironment::bootstrap()
    {
      dw_option_mask.mask("all", percept::LOG_ALWAYS, "log all");

      timer_option_mask.mask("mesh", percept::TIMER_MESH, "mesh operations timers");
    }

    void RunEnvironment::
    setSierraOpts(int procRank, int /*argc*/, char* /*argv*/[])
    {
      const bool m_debug2 = false;

      Teuchos::oblackholestream blackhole;
      std::ostream &out = ( procRank == 0 ? std::cout : blackhole );

      if (m_debug2)
        out << Teuchos::Teuchos_Version() << std::endl << std::endl;

      //clp.setDocString("Run environment options" );
      //clp.setOption("help",         &help_opt,        "help flag");
      //clp.setOption("directory",    &directory_opt,   "working directory");
      //clp.setOption("d",            &directory_opt,   "working directory");
      clp.setOption("output-log",   &output_log_opt,  "output log path");
      clp.setOption("o",            &output_log_opt,  "output log path");
      clp.setOption("pout",         &pout_opt,        "per-processor log file path");
      clp.setOption("dout",         &dout_opt,        "diagnostic output stream one of: 'cout', 'cerr', 'out' or a file path");
      //clp.setOption("dw",           &dw_opt,          dw_option_mask.describe().c_str());
      //clp.setOption("timer",        &timer_opt,       timer_option_mask.describe().c_str());
      //clp.setOption("runtest",      &runtest_opt,     "runtest pid file");;
      //("exit", "do nothing and then exit");

    }

    int processCLP(Teuchos::CommandLineProcessor& clp_in, int procRank, int argc, char* argv[], int &bad_option)
    {
      bad_option = 0;
      const bool m_debug = false;

      Teuchos::oblackholestream blackhole;
      std::ostream &out = ( procRank == 0 ? std::cout : blackhole );

      bool success = true;

      try {

        /* There are also two methods that control the behavior of the
           command line processor.  First, for the command line processor to
           allow an unrecognized a command line option to be ignored (and
           only have a warning printed), use:
        */

        clp_in.recogniseAllOptions(true);

        /* Second, by default, if the parser finds a command line option it
           doesn't recognize or finds the --help option, it will throw an
           std::exception.  If you want prevent a command line processor from
           throwing an std::exception (which is important in this program since
           we don't have an try/catch around this) when it encounters a
           unrecognized option or help is printed, use:
        */
        clp_in.throwExceptions(false);

        /* We now parse the command line where argc and argv are passed to
           the parse method.  Note that since we have turned off std::exception
           throwing above we had better grab the return argument so that
           we can see what happened and act accordingly.
        */
        Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn= Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ;
        bool parse_error = true;
        try {
          parseReturn = clp_in.parse( argc, argv );
          parse_error = false;
          //std::cout << "tmp srk parseReturn = " << parseReturn << std::endl;
        }
        catch (const std::exception & exc)
        {
          parse_error = true;
          out << "RunEnvironment::processCLP error, exc= " << exc.what() << std::endl;
        }
        // figure out which argv option caused the problem and let the user know!
        if (parse_error)
        {
          char *argv1 = argv[1];
          for (int i = 1; i < argc; ++i)
          {
            // pass just one option
            argv[1] = argv[i];
            try 
            {
              parseReturn = clp_in.parse( 2, argv ); // just two args at a time
            }
            catch (...)
            {
               out << "Error processing argument " << i << "=\"" << argv[1] << "\"" << std::endl;
               bad_option = i;
            }
            // restore
            argv[i] = argv[1]; 
         }
         argv[1] = argv1;
         return 1;
        }

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {

          //std::cout << "tmp srk parseReturn = PARSE_HELP_PRINTED " << parseReturn << std::endl;

          return 1;
        }

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION   ) {

          // std::cout << "tmp srk parseReturn = PARSE_UNRECOGNIZED_OPTION " << parseReturn << std::endl;

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

            clp_in.throwExceptions(false);

            clp_in.printHelpMessage(argv[0],out);

            clp_in.throwExceptions(true);
          }

        /*
        // Now we will print the option values
        if (procRank == 0 && m_debug) {
          out << "\nPrinting user options after parsing ...\n\n";
          out << " output_log_opt= " <<  output_log_opt << std::endl;
          out << " dw_opt= " <<  dw_opt << std::endl;
          out << " timer_opt= " <<  timer_opt << std::endl;
          out << " directory_opt= " <<  directory_opt << std::endl;
          out << " help_opt= " <<  help_opt << std::endl;
        }
        */

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
      if (pout_opt.length()) {
        std::string pout_path = pout_opt;
        if (pout_path == "-") {
          std::ostringstream s;

          if (stk::get_log_ostream(out_path))
            s << working_directory << "percept.log." << parallel_size << "." << parallel_rank;
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

    static void panic()
    {
      (void)perror("Unexpected library error");
      std::abort();
    }

#define MAX_PATH 256
    //: Returns the name of the current working directory
    std::string get_working_directory()
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

    /**
     *
     * input fullpath=meshFilename = some file name with a path...
     *
     * output: fullpath = /some/file/path/ending/with/slash/
     * output: meshFileName = ./filename.ext
     *
     */

    int setFileNames(std::string& fullpath, std::string& meshFileName, std::string& errString)
    {
      int err=0;
      fullpath = get_working_directory() + meshFileName;
      if (meshFileName[0] == '/')
        {
          // already absolute
          fullpath = meshFileName;
        }
      else if (meshFileName.length() >= 3 && (meshFileName[0] != '.' && meshFileName[1] != '/'))
        {
          // make it a relative path
          fullpath = get_working_directory() + meshFileName;
          meshFileName = "./"+meshFileName;
          //fullpath = get_working_directory() + meshFileName.substr(2,meshFileName.length()-2);
        }
      else if (meshFileName.length() >= 3 && (meshFileName[0] == '.' && meshFileName[1] == '/'))
        {
          fullpath = get_working_directory() + meshFileName.substr(2,meshFileName.length()-2);
        }
      else
        {
          err = 1;
          errString = "RunEnvironment::setFileNames: bad format for input file name, name= "+meshFileName;
          return err;
        }

      size_t found = 0;
      // strip off the basename + extension
      found = fullpath.find_last_of("/");
      if (found != std::string::npos)
        {
          fullpath = fullpath.substr(0,found);
        }
      if(fullpath == ".")
        fullpath = "./";
      else
        fullpath += "/";

      return err;
    }

    void runCommand(std::string command)
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

    void doLoadBalance(stk::ParallelMachine comm, std::string meshFileName)
    {
      if (meshFileName.length() == 0) return;

      unsigned p_size = stk::parallel_machine_size(comm);
      unsigned p_rank = stk::parallel_machine_rank(comm);
      int err=0;
      std::string errString;

      if (p_size > 1 && !p_rank)
        {
          std::string fullpath = meshFileName;

          err = setFileNames(fullpath, meshFileName, errString);

          if (!err)
            {
              std::string command="decomp --nolaunch -p "+toString(p_size)+ " -R " +fullpath+" " + meshFileName;

              std::cout << "RunEnvironment::doLoadBalance: command= " << command << std::endl;

              runCommand(command);
            }
        }

      stk::all_reduce( comm, stk::ReduceSum<1>( &err ) );

      if (err && !p_rank)
        {
          std::cout << "ERROR in RunEnvironment::doLoadBalance: " << errString << std::endl;
        }
      if (err)
        throw std::runtime_error("ERROR in RunEnvironment::doLoadBalance: " + errString );

#if defined(STK_HAS_MPI)
      MPI_Barrier( comm );
#endif
    }

    void printHelp(Teuchos::CommandLineProcessor& clp_in)
    {
      clp_in.printHelpMessage("RunEnvironment",std::cout);
    }

  const std::string &OptionDescription::get_help( HelpType help_type ) const
  {
    switch (help_type)
    {
      case LIST:
      return keyword;
      break;
      case ONELINE:
      default:
      return help_oneline;
      break;
      case FULL: 
      return help_extra;
      break;
    }
  }

  int get_help(ParserSystem &ps, const std::string &keyword, std::string &help_string, HelpType help_type )
  {
    auto h = ps.advanced_option_help.find(keyword);
    if (h == ps.advanced_option_help.end())
    {
      h      = ps.basic_option_help.find(keyword);
      if (h == ps.basic_option_help.end())
      {
        // keyword not found
        help_string = std::string();
        return 1;
      }
    } 
    // get the string for the keyword
    help_string = h->second.get_help( help_type ); 
    return 0;
  }
 
  } // namespace percept
