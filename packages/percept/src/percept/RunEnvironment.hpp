// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_RunEnvironment_hpp
#define percept_RunEnvironment_hpp

/// copied and edited from stk_util/use_cases/UseCaseEnvironment


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/util/Writer_fwd.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <iosfwd>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

  namespace percept {


    enum LogMask {
      LOG_ALWAYS          = stk::LOG_ALWAYS,
      LOG_TRACE           = stk::LOG_TRACE,
      LOG_TRACE_STATS     = stk::LOG_TRACE_STATS,
      LOG_TRACE_SUB_CALLS = stk::LOG_TRACE_SUB_CALLS,
      LOG_MEMBERS         = stk::LOG_MEMBERS,

      LOG_APPLICATION     = 0x0000010  // use this as the base for additional masks
      //LOG_SEARCH          = 0x0000010,
      //LOG_TRANSFER        = 0x0000020,
      //LOG_TIMER           = 0x0000040

    };

    /**
     * @brief Class <b>message_type</b> ...
     *
     */
    enum message_type {
      MSG_WARNING = stk::MSG_WARNING,
      MSG_FATAL   = stk::MSG_DOOMED,
      MSG_INFORMATION,
      MSG_EXCEPTION,
      MSG_PARALLEL_EXCEPTION
    };


    /**
     * @brief Class <b>type</b> ...
     *
     */
    enum message_throttle_type {
      MSG_APPLICATION = stk::MSG_APPLICATION,
      MSG_TIME_STEP
    };

    enum TimerSetMask {
      TIMER_MESH     = 0x00000001,		///< Enable mesh timers
      //      TIMER_MESH_IO  = 0x00000002,		///< Enable mesh I/O timers
      //      TIMER_SEARCH   = 0x00000004,		///< Enable search timers
      //    TIMER_TRANSFER = 0x00000008,		///< Enable transfer timers
      TIMER_ALL      = 0xFFFFFFFF,		///< Force timer to be active

      TIMER_FORCE    = 0x00000000		///< Force timer to be active
    };

    std::ostream &out();                ///< Normal output stream
    std::ostream &dout();               ///< Diagnostic output stream
    std::ostream &pout();               ///< Per-processor output stream (See RuntimeDeferredx)
    std::ostream &tout();               ///< Regression test textual output stream

    std::ostream &dwout();              ///< Diagnostic writer stream

    // dw() definition
    stk::diag::Writer &dw();
#define DWENDL stk::diag::dendl

    stk::diag::TimerSet &timerSet();

    stk::diag::Timer &timer();

    void my_report_handler(const char *message, int type);

    // this little class is simply here to force an ordering in the ~RunEnvironment() dtor, so that
    //   parallel_machine_finalize gets called after m_comm is destroyed; but, it still only works
    //   if this is invoked after returning from an enclosing block of RunEnvironment.  why?
    class ParallelMachineFinalize
    {
      bool m_need_to_finalize;
    public:
      ParallelMachineFinalize(bool need_to_finalize=false) : m_need_to_finalize(need_to_finalize) {}
      ~ParallelMachineFinalize()
      {
        if (m_need_to_finalize)
          {
            stk::parallel_machine_finalize();
          }
      }
    };

    class RunEnvironment : public ParallelMachineFinalize
    {

    public:
      // Will initialize a comm
      RunEnvironment(int *argc, char ***argv, bool debug=false);

      // Assumes already-initialized comm
      RunEnvironment(int *argc, char ***argv, stk::ParallelMachine comm, bool debug=false);

      ~RunEnvironment();

      std::string
      build_log_description(const std::string &           working_directory,
                            int                           parallel_rank,
                            int                           parallel_size);

      int get_argc() { return m_argc; }
      char **get_argv() { return m_argv; }



      // command line options
      Teuchos::CommandLineProcessor clp;

      std::string output_log_opt;
      std::string dw_opt;
      std::string timer_opt;
      std::string directory_opt;

      std::string pout_opt;
      std::string dout_opt;
      std::string runtest_opt;

      int help_opt;

      // data
      const stk::ParallelMachine    m_comm;
      static std::string            m_workingDirectory;

    private:

      bool                          m_need_to_finalize;
      bool                          m_debug;
      bool                          m_processCommandLine_invoked;
      std::string                  *m_argv_new;
      int                           m_argc;
      char                        **m_argv;

      //ParallelMachineFinalize       m_par_finalize;

      // shared constructor implementation; do not call directly
      void internal_initialize(int argc, char **argv);
      void bootstrap();

      void setSierraOpts(int procRank, int argc, char* argv[]);
    };

    int processCommandLine(Teuchos::CommandLineProcessor& clp_in, int argc, char **argv, int &bad_option);
    inline
    int processCommandLine(Teuchos::CommandLineProcessor& clp_in, int argc, char **argv)
    {int bad_option(0); return processCommandLine(clp_in, argc, argv, bad_option);}

    int processCLP(Teuchos::CommandLineProcessor& clp_in, int procRank, int argc, char* argv[], int &bad_option);

    std::string get_working_directory();

    int setFileNames(std::string& fullmesh, std::string& meshFileName, std::string& errString);

    void runCommand(std::string command);

    void doLoadBalance(stk::ParallelMachine comm, std::string meshFileName);

    void printHelp(Teuchos::CommandLineProcessor& clp_in);

    // interface to Teuchos_CommandLine setOption
    const int num_operations=5;
    enum OptionType {HELP=0, IO, OPERATION, PARAMETER, INFO = num_operations-1};
    const std::string operation_name[num_operations] = {"HELP", "IO", "OPERATION", "PARAMETER", "INFO"};
    const int num_categories=3;
    enum Category {SIMPLE=0, COMMON, ALL=num_categories-1};
    const std::string category_name[num_categories] = {"SIMPLE", "COMMON", "ALL"};

    // help strings for each keyword
    enum HelpType {LIST, ONELINE, FULL};
    class OptionDescription
    {
    public:
      std::string keyword, help_oneline, help_extra;
      OptionType op = HELP;
      OptionDescription() {}
      OptionDescription(const OptionDescription& that) 
      { keyword = that.keyword; help_oneline = that.help_oneline; help_extra = that.help_extra; op = that.op; }
      OptionDescription(const std::string &kw, OptionType opt, 
        const std::string &help_ol, const std::string &help_ext ) :
        keyword(kw), help_oneline(help_ol), help_extra(help_ext), op(opt) {}
      const std::string &get_help(HelpType help_type = ONELINE) const;
    };  

    class ParserSystem
    {
    public:
      Teuchos::CommandLineProcessor clp;
      std::map<std::string,OptionDescription> basic_option_help; 
      std::map<std::string,OptionDescription> advanced_option_help;
    };
  
    // returns 0 and non-empty strings if the keyword was found
    int get_help(ParserSystem &ps, const std::string &keyword, std::string &help_string, HelpType help_type );

    // variable_string, convert a variable to strings of its type name and value
    // overloaded, rather than template, since the functions to converted to a string are type-specific
  inline void variable_string(const std::string *var, std::string &variable_type, std::string &value_string)
  { variable_type = "string"; value_string = *var; }
  inline void variable_string(const int         *var, std::string &variable_type, std::string &value_string)
  { variable_type = "int";    value_string = std::to_string(*var); }
  inline void variable_string(const double      *var, std::string &variable_type, std::string &value_string)
  { variable_type = "double"; value_string = std::to_string(*var); }

  // right pad with spaces up to desired length
  inline void pad( std::string &s, int desired_length )
  {
    int padl = desired_length - (int) s.length();
    if (padl>0)
      s += std::string( padl, ' ' );
  }

  template <class VarType>
  void set_option(
    ParserSystem &ps,
    const char     option_name[],
    const VarType       &variable,
    const OptionType    &op,
    const std::string   &doc_base,
    const std::string   &doc_oneline = "",
    const std::string   &doc_extra   = "",
    bool show_default = false,
    bool is_advanced  = true
    )
    {

    std::string variable_type_string, value_string;
    variable_string(variable, variable_type_string, value_string);

    // right-pad keyword so they are all the same length
    const int max_keyword_length = 30; // based on longest keyword
    std::string keyword(option_name), keyword_padded( std::string("--") + option_name);
    pad(keyword_padded, max_keyword_length);

    const int max_value_length = 8;  // based on longest "double", "string", "int", ...
    pad( variable_type_string, max_value_length );

    // find longest "operation" name, e.g. "PARAMETER"
    int max_op_length = 0;
    for (int i = 0; i < num_operations; ++i)
      max_op_length = std::max( max_op_length, (int) operation_name[i].length() );
    std::string ops( operation_name[op] );
    pad(ops, max_op_length+2);

    std::string doc_prefix = keyword_padded + variable_type_string + ops;
    int indent_length = doc_prefix.length();
    doc_prefix += doc_base;
    const std::string ds = (show_default ? std::string("\n") + std::string( indent_length, ' ') + std::string("Default=") + value_string : std::string());
    std::string help_oneline = doc_prefix + doc_oneline + ds + "\n";
    std::string help_extra = doc_prefix + ds + doc_extra + "\n";

    if (is_advanced)
      ps.advanced_option_help[keyword] = OptionDescription( keyword, op, help_oneline, help_extra);
    else
      ps.basic_option_help[keyword]    = OptionDescription( keyword, op, help_oneline, help_extra);

    ps.clp.setOption(option_name, variable, (doc_base + doc_oneline).c_str());
  }


  } // namespace percept

#endif // percept_RunEnvironment_hpp
