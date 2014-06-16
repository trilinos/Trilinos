/*--------------------------------------------------------------------*/
/*    Copyright 2003 - 2010 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_Env_h
#define STK_UTIL_DIAG_Env_h

#include <ios> // #include <iosfwd> // ios_base is not defined in the forward file
#include <string>

#include <stk_util/stk_config.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif
#include <stk_util/parallel/Parallel.hpp>

#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/util/Bootstrap.hpp>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/Option.hpp>
#include <stk_util/diag/Writer_fwd.hpp>

namespace sierra {

/// @brief Namespace <b>Env</b> contains the runtime environment bootstrap for the
/// MPI, logging, command line argument parsing, runtime information and signal handling.
///
/// <H3>MPI Initialization</H3>
///
/// <H3>Output Log File</H3>
///
/// <H3>Command Line Options</H3>
///
/// <H3>Runtime Information</H3>
///
/// <H3>Signals and Long Jump</H3>
///
namespace Env {

static const std::string PARAM_ON = "on";  ///< Option value when command line option specified without a parameter

///
/// @addtogroup EnvDetail
/// @{
///

/**
 * @brief Enumeration ExecutableType defines the known types of coordinated executables that operate
 * with a sierra application.  Unfortunately, this scheme for coordination is currently defined by
 * Gemini whose implementation forces a limit of two executables, namely it and a fluid code.  The
 * startup_multi_exec() function handles the creation of groups which are contiguous processor
 * groups, each with lead processor being the least ranked processor in the group.
 *
 * Modification of the startup_multi_exec() function would need to be made to enable more than the
 * two executable types.
 */
enum ExecType {
  EXEC_TYPE_WORLD = 0,            ///< Generic application using entire communicator (MPI_COMM_WORLD)
  EXEC_TYPE_FLUID = 1,            ///< Gemini Euler application
  EXEC_TYPE_LAG   = 2,            ///< Sierra Lagrangian application
  EXEC_TYPE_PEER  = 3             ///< Split communicator application; non-Gemini
};

struct ExecInfo
{
  MPI_Comm              m_groupComm;
  int                   m_master;
};



  /**
   * @brief Initialize MPI related operations for sierra, outputs banner, etc.  
   *        returns 1 if MPI was initialized, 0 otherwise
   *
   * @param argc		an <b>int</b> pointer to the main argc.
   *
   * @param argv		a <b>char</b> pointer to the main argv.
   *
   * @param product_name	a <b>char</b> const pointer to the name of the
   *				product.
   *
   * @param build_date_time	a <b>char</b> const pointer to
   *				__DATE__ " " __TIME__
   *
   * @param mpi_key	        an optional <b>ExecType</b> enumeration
   *                            specifying the type of executable.
   *                            Default is a single executable using
   *                            MPI_COMM_WORLD. Other options result
   *                            in a split communicator.
   *
   * @param peer_sizes          an optional <b>std::vector<int></b>
   *                            const pointer containing the number of
   *                            processors that each peer communicator
   *                            (\sa parrallel_peer_comm()) will
   *                            support. The number of peers is
   *                            determined by the size of the vector.
   *                            Only used if mpi_key = EXEC_TYPE_PEER
   */
 bool StartupSierra(int *argc, char ***argv, const char *product_name, const char *build_date_time,
	           ExecType mpi_key = EXEC_TYPE_WORLD, const std::vector<int> *peer_sizes = NULL);

 //
 //  Cleanup any MPI stuff that was created from the startup sierra call, pass in the MPI initialization flag
 //  that StartupSierra generated
 //
 void ShutDownSierra(bool mpiInitFlag);

  
/**
 * @ingroup EnvMPIDetail EnvOutputDetail EnvCommandLineDetail
 * @brief Class <b>Startup</b> is a sentry class for starting the application.  It
 * ensures that the command line arguments, platform and MPI are ready to go at the start
 * of the application.
 *
 */
class Startup
{
public:

  /**
   * @brief Creates a new <b>Startup</b> instance.
   *
   * @param argc		an <b>int</b> pointer to the main argc.
   *
   * @param argv		a <b>char</b> pointer to the main argv.
   *
   * @param product_name	a <b>char</b> const pointer to the name of the
   *				product.
   *
   * @param build_date_time	a <b>char</b> const pointer to
   *				__DATE__ " " __TIME__
   *
   * @param mpi_key	        an optional <b>ExecType</b> enumeration
   *                            specifying the type of executable.
   *                            Default is a single executable using
   *                            MPI_COMM_WORLD. Other options result
   *                            in a split communicator.
   *
   * @param peer_sizes          an optional <b>std::vector<int></b>
   *                            const pointer containing the number of
   *                            processors that each peer communicator
   *                            (\sa parrallel_peer_comm()) will
   *                            support. The number of peers is
   *                            determined by the size of the vector.
   *                            Only used if mpi_key = EXEC_TYPE_PEER
   */
  Startup(int *argc, char ***argv, const char *product_name, const char *build_date_time,
	  ExecType mpi_key = EXEC_TYPE_WORLD, const std::vector<int> *peer_sizes = NULL);

  /**
   * @brief Destroys a <b>Startup</b> instance.  IT closes all logging output streams and
   * will finalize MPI only if the Startup::Startup initialized the MPI.
   *
   */
  ~Startup();

private:
  /**
   * @brief Member function <b>startup</b>
   *  @li initializes MPI if not already initialized,
   *  @li queries the environment and command line for environment options,
   *  @li broadcasts the environment options from MPI_COMM_WORLD processor #0 to all other
   *      processors,
   *  @li opens the default parallel output stream.
   *
   * @param argc		an <b>int</b> pointer to the main argc.
   *
   * @param argv		a <b>char</b> pointer to the main argv.
   *
   * @param product_name	a <b>char</b> const pointer to the name of the
   *				product.
   *
   * @param build_time		a <b>char</b> const pointer to __DATE__ " " __TIME__
   *
   * @param exec_options	an <b>OptionMap</b> const reference to additional
   *				options available for the command line/environment.
   *
   */
  void startup(int *argc, char ***argv, const char *product_name, const char *build_time,
	       ExecType mpi_key, const std::vector<int> *peer_sizes);
  
private:
  bool m_mpiInitFlag;			///< True if Startup initialized MPI
};


/**
 * @ingroup EnvMPIDetail
 * @brief Function <b>reset</b> determines new parallel_size and parallel_rank.
 * Flushes, closes, and reopens log files.
 *
 * This is a very dangerous method, if the previous communicator had been queried and
 * given to an external library then Sierra and that external library will be out-of-sync.
 * It is strongly recommended that parser-interface singletons do not query and save the
 * communicator.
 *
 * This capability is provided to support dakota.
 *
 * @param mpi_comm		a <b>MPI_Comm</b> value for the new mpi
 *				communicator.
 *
 * @param work_dir		a <b>char</b> const pointer to the new working
 *				directory.
 *
 */
void reset(MPI_Comm mpi_comm); // , const char * const work_dir = NULL);
void setMpiCommunicator(MPI_Comm communicator);
bool is_comm_valid();


/**
 * @ingroup EnvCommandLineDetail
 * @brief Function <b>query_env_param</b> searches the command line options for the
 * specified option.  If not found, a const reference to an empty string is returned.
 *
 * @param option		a <b>char</b> const pointer to the option to
 *				search for.
 *
 * @return			a <b>std::string</b> const reference to the options value
 *				or an empty string if not found.
 */
const std::string &get_param(const char * const option);


/**
 * @brief Member function <b>set_param</b> assigns the value to the parameter option.
 *
 * @param option		a <b>char</b> const pointer to the option to
 *				search for.
 *
 * @param value			a <b>char</b> const pointer to the value to assign to the
 *				param.
 *
 */
void set_param(const char *option, const std::string &value);

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>product_name</b> returns the product's name.  This name is
 * used to query the product registry for information concerning this product.
 *
 * @return			a <b>std::string</b> const reference to the product's name.
 */
const std::string &product_name();


/**
 * @ingroup EnvRuntimeInformationDetail 
 * @brief Function <b>developer_mode</b> returns true if the
 * --developer option was specified on the application command line.
 *
 * @return			a <b>std::string</b> const reference to the applications
 *				operating platform string.
 */
bool developer_mode();

/**
 * @brief Function <b>set_input_file_required<b> sets whether lack of an input
 * file specification will automatically cause failure.  The default behavior
 * corresponds to true.
 */
void set_input_file_required(bool value);

void setInputFileName(std::string value);
std::string getInputFileName();


/**
 * @brief Function <b>set_check_subcycle<b> sets whether to check input
 * file for subcycling.  The default behavior corresponds to false.
 */
void set_check_subcycle(bool value);

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>architecture</b> returns the platform executing this product.
 * This is obtained during startup by searching for a file which contains this
 * information.
 *
 * @return			a <b>std::string</b> const reference to the applications
 *				operating platform string.
 */
const std::string &architecture();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>executable_file</b> returns the path of this executable file.
 * information.
 *
 * @return			a <b>std::string</b> const reference to the executable
 *				file path string.
 */
const std::string &executable_file();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>executable_date</b> returns the build date of the executable
 * file as a string in the form Mmm dd yyyy hh:mm::ss.
 *
 * @return			a <b>std::string</b> const reference to the executable
 *				file's build date and time.
 */
const std::string &executable_date();

/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>startup_date</b> returns the startup date of this application
 * execution.
 *
 * @return			a <b>std::string</b> const reference to the application
 *				execution start date and time.
 */
const std::string &startup_date();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Function <b>start_time</b> returns the start time of this application
 * execution.
 *
 * @return			a <b>double</b> value of the application execution
 *				start time.
 */
double start_time();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Member function <b>wall_now</b> returns the epoch as a double precision
 * value in seconds to "millisecond" accuracy.
 *
 * @return a <b>double</b> ...
 */
double wall_now();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Member function <b>cpu_now</b> returns the accumlated cpu time for the
 * process as a double precision value in seconds to "millisecond" accuracy.
 *
 * @return a <b>double</b> ...
 */
double cpu_now();


/**
 * @ingroup EnvRuntimeInformationDetail
 * @brief Member function <b>vm_now</b> returns the virtual memory in use for the
 * process as a double precision value in bytes.
 *
 * @return			a <b>double</b> value of the number of bytes of virtual
 *				memory currently in use.
 */
double vm_now();


/**
 * @ingroup EnvRuntimeInformationDetail, EnvOutputDetail
 * @brief Function <b>working_directory</b> returns the current working directory of
 * this application execution.
 *
 * @return			a <b>std::string</b> to the application
 *				execution working directory.
 */
const std::string working_directory();


/**
 * @ingroup EnvOutputDetail
 * @brief Function <b>output</b> returns the processor output log stream.  This
 * stream is connected via an mpi_filebuf to processor 0.  Upon
 * <b>output_flush()</b> the output from all processors is collected on processor 0
 * in a sequential by process and is logged to output file in a non-jumbled manner.
 *
 * @return			a <b>std::ostream</b> reference to the application
 *				log stream.
 */
std::ostream &output();


/**
 * @ingroup EnvOutputDetail
 * @brief Function <b>outputP0</b> returns the processor output log stream on
 * processor 0 and the null log stream on all other processors.
 *
 * @return			a <b>std::ostream</b> reference to the application
 *				log stream on processor 0 or null stream on the
 *				processors.
 */
std::ostream &outputP0();


/**
 * @ingroup EnvOutputDetail
 * @brief Function <b>outputNull</b> returns the null output stream.  All data is
 * simply discarded by the buffer associated with this stream.
 *
 * @return			a <b>std::ostream</b> reference to the null output
 *				stream.
 */
std::ostream &outputNull();


/**
 * @ingroup EnvOutputDetail
 * @brief Function <b>output_open</b> opens an output file on processor zero for
 * synchronous data output from all processors via the <b>mpi_filebuf</b> class.
 * The output is synchronized via the <b>output_flush()</b> function and maintain in
 * the output stream list so that it is flushed and closed on application rundown.
 *
 * Must be executed concurrently on all processor in the current
 * <b>Env::parallel_comm()</b> group.
 *
 * @param filename		a <b>char</b> const pointer to the path of file to
 *				open.
 *
 * @return			a <b>std::ostream</b> reference to the newly opened
 *				output stream.
 */
//std::ostream &output_open(const char * const filename);


/**
 * @brief Function <b>section_separator</b> returns a c-style string to be used as a
 * output section separator.
 *
 * @return			a <b>char</b> const pointer to the section separator
 *				string.
 */
const char *section_separator();


/**
 * @brief Function <b>subsection_separator</b> returns a c-style string to be used as a
 * output subsection separator.
 *
 * @return			a <b>char</b> const pointer to the subsection separator
 */
const char *subsection_separator();


/**
 * @brief Function <b>section_title</b> returns a section title.  The title has date and
 * time concatenated and right justified to the length of the <b>section_separator</b>.
 * The date and time is 20 characters wide, so adjust your titles accordingly.
 *
 * @param title			a <b>std::string</b> const reference to the title string.
 *
 * @return			a <b>std::string</b> value with the date and time right
 *				justified to <b>section_separator</b>'s length.
 */
std::string section_title(const std::string &title);


/**
 * @ingroup EnvOutputDetail
 * @brief Function <b>output_flush</b> flushes all output on all currently open
 * synchronous outptu files which were opened via <b>output_open</b>.
 *
 * Must be executed concurrently on all processor in the current
 * <b>Env::parallel_comm()</b> group.
 *
 */
void output_flush();


/**
 * @ingroup EnvOutputDetail
 * @brief Function <b>output_flush</b> synchronously flushes <b>stream</b>
 * which was created using an <b>mpi_filebuf</b>.
 *
 * Must be executed concurrently on all processor in the current
 * <b>Env::parallel_comm()</b> group.
 *
 * @param stream		a <b>std::ostream</b> reference to the synchronous
 *				output stream to flush.
 */
void output_flush(std::ostream &stream);

void request_shutdown(bool shutdown = true);

bool is_shutdown_requested();

/**
 * @ingroup EnvRuntimeInformationDetail EnvMPIDetail EnvOutputDetail
 * @brief Function <b>abort</b> aborts the execution of the sierra application.
 *
 */
void abort();


/**
 * @brief Function <b>parallel_comm</b> returns the current MPI communicator used by
 * the sierra environment.
 *
 * @return			a <b>MPI_Comm</b> value of the current MPI
 *				communicator.
 */
MPI_Comm parallel_comm();

/**
 * @brief Function <b>parallel_world_comm</b> returns the MPI_COMM_WORLD communicator used by
 * the sierra environment in a MPMD parallel application.
 *
 * @return			a <b>MPI_Comm</b> value of the current MPI_COMM_WORLD
 *				communicator.
 */
MPI_Comm parallel_world_comm();


/**
 * @brief Function <b>peer_group</b> returns the peer group rank for an application of type
 *        EXEC_TYPE_PEER.
 *
 * @return			a <b>int</b> value of the peer group for the peer application.
 */
int peer_group();


/**
 * @brief Function <b>parallel_lag_master</b> returns the global rank of the Gemini Euler application.
 *
 * @return			a <b>int</b> value of the global rank of the Gemini Euler
 *				application.
 */
int parallel_fluid_master();


/**
 * @brief Function <b>parallel_lag_master</b> returns the global rank of the Sierra lagrangian application.
 *
 * @return			a <b>int</b> value of the global rank of the Sierra lagrangian
 *				application.
 */
int parallel_lag_master();


/**
 * @ingroup EnvMPIDetail
 * @brief function <b>parallel_size</b> returns the number of processors
 * in the current mpi communicator.
 *
 * @return			an <b>int</b> value of the number of processors in
 *				the current mpi communicator.
 */
int parallel_size();


/**
 * @ingroup EnvMPIDetail
 * @brief function <b>parallel_rank</b> returns the rank of this processor
 * in the current mpi communicator.
 *
 * @return			an <b>int</b> value of the rank for this processors
 *				in the current mpi communicator.
 */
int parallel_rank();

// /**
//  * @brief Function <b>set_current_diag_stream</b> set the diagnostic writer current out
//  * put stream.  The special names - and cout attach the diag stream to the std::cout
//  * stream.  The special name cerr attaches the stream to the std::cerr stream, output
//  * attaches to the stream Env::output() and outputp0 to the stream outputP0().  Otherwise,
//  * the file is opened with <b>mode</b> ios mode flags and attached to the diagnostic
//  * stream.
//  *
//  * @param path			a <b>char</b> const pointer to the file path to
//  *				open.
//  *
//  * @param mode			a <b>std::ios::openmode</b> value of the ios mode to
//  *				open the file.
//  *
//  */
// void set_current_diag_stream(const char *path, std::ios_base::openmode mode = std::ios_base::out);

// /**
//  * @brief Function <b>getCurrentStream</b> returns the current diagnostic stream.
//  *
//  * @return			a <b>std::ostream</b> reference to the current diagnostic
//  *				stream.
//  */
// // Diag::Stream &get_current_diag_stream();
// std::ostream &get_current_diag_stream();

///
/// @}
///

} // namespace Env
} // namespace sierra

#endif // STK_UTIL_DIAG_Env_h
