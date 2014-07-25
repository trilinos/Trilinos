/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_percept_RunEnvironment_hpp
#define stk_percept_RunEnvironment_hpp

/// copied and edited from stk_util/use_cases/UseCaseEnvironment 


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/diag/Writer_fwd.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <iosfwd>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

namespace stk_classic {
  namespace percept {


    enum LogMask {
      LOG_ALWAYS          = stk_classic::LOG_ALWAYS,
      LOG_TRACE           = stk_classic::LOG_TRACE,
      LOG_TRACE_STATS     = stk_classic::LOG_TRACE_STATS,
      LOG_TRACE_SUB_CALLS = stk_classic::LOG_TRACE_SUB_CALLS,
      LOG_MEMBERS         = stk_classic::LOG_MEMBERS,

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
      MSG_WARNING = stk_classic::MSG_WARNING,
      MSG_FATAL   = stk_classic::MSG_DOOMED,
      MSG_INFORMATION,
      MSG_EXCEPTION,
      MSG_PARALLEL_EXCEPTION
    };


    /**
     * @brief Class <b>type</b> ...
     *
     */
    enum message_throttle_type {
      MSG_APPLICATION = stk_classic::MSG_APPLICATION,
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
    stk_classic::diag::Writer &dw();
#define DWENDL stk_classic::diag::dendl

    stk_classic::diag::TimerSet &timerSet();

    stk_classic::diag::Timer &timer();

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
            stk_classic::parallel_machine_finalize();
          }
      }
    };

    class RunEnvironment : public ParallelMachineFinalize
    {

    public:
      // Will initialize a comm
      RunEnvironment(int *argc, char ***argv, bool debug=false);

      // Assumes already-initialized comm
      RunEnvironment(int *argc, char ***argv, stk_classic::ParallelMachine comm, bool debug=false);

      int processCommandLine(int argc, char **argv);
      int processCommandLine() { return processCommandLine(m_argc, m_argv); }

      ~RunEnvironment();

      void printHelp();

      static void 
      doSierraLoadBalance(stk_classic::ParallelMachine comm, std::string meshFileName);

      static void
      doLoadBalance(stk_classic::ParallelMachine comm, std::string meshFileName);

      std::string
      build_log_description(const std::string &           working_directory,
                            int                           parallel_rank,
                            int                           parallel_size);

      int get_argc() { return m_argc; }
      char **get_argv() { return m_argv; }

      static std::string
      get_working_directory();


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
      const stk_classic::ParallelMachine    m_comm;
      static std::string            m_workingDirectory;

    private:

      static int setFileNames(std::string& fullmesh, std::string& meshFileName, std::string& errString);

      bool                          m_need_to_finalize;
      bool                          m_debug;
      bool                          m_processCommandLine_invoked;
      std::string                  *m_argv_new;
      int                           m_argc;
      char                        **m_argv;

      //ParallelMachineFinalize       m_par_finalize;

      int  processCLP(int procRank, int argc, char* argv[]);
      // shared constructor implementation; do not call directly
      void internal_initialize(int argc, char **argv);
      void bootstrap();

      void setSierraOpts(int procRank, int argc, char* argv[]);
    };

  } // namespace percept
} // namespace stk_classic

#endif // stk_percept_RunEnvironment_hpp
