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

#ifndef STK_UTIL_DIAG_Env_h
#define STK_UTIL_DIAG_Env_h

#include <stk_util/stk_config.h>
#if defined ( STK_HAS_MPI )
#  include <mpi.h>                        // for MPI_Comm
#endif
#include <ios>                          // for ostream
#include <string>                       // for string



namespace sierra {

/**
 * @brief Function <b>format_time</b> encodes the time using the format specified.
 * The format is described in <b>stdftime</b>.
 *
 * @param t		a <b>time_t</b> value of the time to format.
 *
 * @param format	a <b>char</b> const pointer to the format.
 *
 * @return		a <b>String</b> value of the encoded time.
 */
std::string format_time(double t, const char *format = "%b %e %Y %H:%M:%S");

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


//
//  Provide single switch point for gemini interface options
//
  enum GeminiSCIVersion {
    GEMINI_SCI_UNKNOWN = 0,
    GEMINI_SCI_1 = 1,
    GEMINI_SCI_2 = 2
  };

  GeminiSCIVersion GetGeminiVersion(GeminiSCIVersion ver=GEMINI_SCI_UNKNOWN);



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
 * @brief Function <b>set_zapotec<b> sets whether this code is
 * zaptoec.  The default behavior corresponds to false.
 * Send all function-related hate mail to Arne Gullerud.
 */
void set_zapotec(bool value);

/**
 * @brief Function <b>is_zapotec<b> returns whether this code is
 * zaptoec. Send all function-related hate mail to Arne Gullerud.
 */
bool is_zapotec();

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
