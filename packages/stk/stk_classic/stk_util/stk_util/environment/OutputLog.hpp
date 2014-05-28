/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_OUTPUTLOG_HPP
#define STK_UTIL_ENVIRONMENT_OUTPUTLOG_HPP

#include <iosfwd>
#include <string>

// #include <stk_util/util/TeeStreambuf.hpp>

namespace stk {

template<class Ch, class Tr>
class basic_tee_streambuf;

/// Tee stream buffer for char
typedef stk::basic_tee_streambuf<char, std::char_traits<char> > tee_streambuf;

///
/// @addtogroup output_log_detail
/// @{
///

/**
 * @file
 *
 * The logging is implemented using two maps.  One maps names to log files and streams and one maps
 * names to output streams.  Log files and streams are the ultimate destinations for output.  The
 * output streams are the data providers.  The output streams may also serve as destinations
 * streams, passing stream data to their destinations.
 *
 * To implement this strategy, the output streams, once registered, create a tee streambuf and
 * attach this as the rdbuf of the registered stream.  The original streambuf becomes the
 * destination of the tee streambuf.  The tee streambuf sinks data from a stream and sources it to
 * all of its destination streambufs just as the tee utility tees input from stdin to a file set and
 * stdout.
 *
 * The tee streambuf class provides functions to add and remove destination streams.
 *
 * When a output stream is unregistered, the tee streambuf is destroyed and the streams original
 * streambuf is restored.
 *
 * The bind_output_streams() function provides string commands to control the sinks and sources of
 * the log file streams and the output streams.
 *
 * The functions using ostreams could be converted to templates of Ch and std::char_traits<Ch>.
 */

/**
 * @brief Function <b>register_ostream</b> registers an output stream with the output stream
 * registry.  The registration process creates an intermediate tee streambuf.
 *
 * @param name			a <b>std::string</b> const reference to the name of the output
 *                              stream. 
 *
 * @param output_stream		a <b>std::ostream</b> reference to the output stream to register.
 *
 */
void register_ostream(std::ostream &os, const std::string &name);

/**
 * @brief Function <b>unregister_ostream</b> unregisters an output stream.
 *
 * @param output_stream		a <b>std::ostream</b> reference to the output stream to unregister.
 *
 */
void unregister_ostream(std::ostream &os);

/**
 * @brief Function <b>bind_output_streams</b> parses the output_description and opens and registers
 * the log streams and binds the registered output streams to the registered log streams.
 *
 * The output description is defined as a white space separated string of command phrases.  There
 * are two types of command phrases, log file creation and output stream direction.
 *
 * The log file creation phrase consists of a name, and equal sign (=) and a file path.  The file
 * path is opened and the log file stream is registered with the log streams.
 *
 * The output stream phrase consists of a name, a greater than sign (>) and a output stream
 * selector.  The selector is a list of log streams or output streams which may be prefixed with a
 * plus sign (+) or a minus sign (-).  No prefix removes all current log and output streams
 * from the named output stream before adding the new log or output stream. A plus prefix adds the
 * new log or output stream to the named output stream and a minus prefix removes the log or output
 * stream from the named output stream.
 *
 * @param output_description	a <b>std::string</b> const reference to the output desciption.
 *
 */
void bind_output_streams(const std::string &output_description);

/**
 * @brief Function <b>register_log_ostream</b> takes an existing std::ostream and makes it available
 * for output redirection.
 *
 * @param output_stream		a <b>std::ostream</b> reference to the output stream to register.
 *
 * @param name			a <b>std::string</b> const reference to the name of this log
 *                              stream. 
 *
 */
void register_log_ostream(std::ostream &os, const std::string &name);

/**
 * @brief Function <b>register_log_ostream</b> takes an existing std::ostream and makes it available
 * for output redirection.
 *
 * @param output_stream		a <b>std::ostream</b> reference to the output stream to register.
 *
 */
void unregister_log_ostream(std::ostream &os);

/**
 * @brief Function <b>create_log_file</b> opens a log file at the specified path and adds it to the
 * registry of log files with the specified name.  This name is be used at to locate the path of the
 * file or it's output stream from the registry using the get_log_file_path() and
 * get_log_file_ostream() functions.
 *
 * @param name			a <b>std::string</b> const reference to the name to give the log
 *                              file.
 *
 * @param path			a <b>std::string</b> const reference to the path of the log file to
 *                              create. 
 *
 */
void create_log_file(const std::string &name, const std::string &path);

/**
 * @brief Function <b>close_log_file</b> close the log file with the specified name and
 * removes it from the registry of log files.
 *
 * @param name			a <b>std::string</b> const reference to the name of the log file to
 *                              close. 
 *
 */
void close_log_file(const std::string &name);

/**
 * @brief Function <b>is_registered_ostream</b> returns true if an output stream of the
 * specified name is registered.
 *
 * @param name			a <b>std::string</b> const reference to the output stream to test
 *                              existence in the output stream registry.
 *
 * @return			a <b>bool</b> value of true if the specified output stream is in the
 *                              registry. 
 */
bool is_registered_ostream(const std::string &name);

/**
 * @brief Function <b>get_log_path</b> returns the file path of the log file with the specified name
 * from the log file registry.  If the specified name does not exist in the registry, an empty
 * string is returned.
 *
 * @param name			a <b>std::string</b> const reference to the name of the log file to
 *                              return its file path.
 *
 * @return			a <b>std::string</b> const reference to the file path of the log
 *                              file.
 */
const std::string &get_log_path(const std::string &name);

/**
 * @brief Function <b>get_log_file_ostream</b> return the output stream of the log file with the
 * specified name from the log registry.  If the specified name does not exist in the registry, a
 * null pointer is returned.
 *
 * @param name			a <b>std::string</b> const reference to the name of the log file to
 *                              return its output stream.
 *
 * @return			a <b>std::ostream</b> pointer to the output stream of the log file.
 */
std::ostream *get_log_ostream(const std::string &name);

/**
 * @brief Function <b>get_ostream_streambuf</b> locates the output stream registered with
 * the specified name.  If the specified output stream does not exist in the registry, a null
 * pointer is returned.
 *
 * @param name			a <b>std::string</b> const reference of the name of the outputstream
 *                              to return its tee streambuf.
 *
 * @return			a <b>std::ostream</b> pointer to the output stream.
 */
std::ostream *get_ostream_ostream(const std::string &name);

/**
 * @brief Function <b>get_ostream_tee_streambuf</b> locates the tee streambuf registered
 * with the specified name.  If the specified output stream does not exist in the registry, a null
 * pointer is returned.
 *
 * @param name			a <b>std::string</b> const reference of the name of the outputstream
 *                              to return its tee streambuf.
 *
 * @return			a <b>tee_streambuf</b> pointer to the tee streambuf.
 */
std::ostream *get_ostream_tee_ostream(const std::string &name);

///
/// @}
///

} // namespace stk

namespace sierra {

std::ostream &out();                ///< Normal output stream
std::ostream &dout();               ///< Diagnostic output stream
std::ostream &pout();               ///< Per-processor output stream (See RuntimeDeferredx)
std::ostream &tout();               ///< Regression test textual output stream

std::ostream &dwout();              ///< Diagnostic writer stream

} // namespace sierra


#endif // STK_UTIL_ENVIRONMENT_OUTPUTLOG_HPP
