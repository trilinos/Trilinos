/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_RUNTIMEMESSAGE_HPP
#define STK_UTIL_ENVIRONMENT_RUNTIMEMESSAGE_HPP

#include <iosfwd>
#include <vector>
#include <cstddef>

#include <stk_util/parallel/Parallel.hpp>

namespace stk {

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @file
 */

/**
 * @brief Typedef <b>MessageId</b> defines a message identifier.
 *
 * Message identifiers must be consist from reference to reference, unique for each instance, and
 * yet consist within each instance across multiple processors.  To meet these criteria, the message
 * identifier is implemented as a static memory location.  It must be declared by the application
 * developer as static.  This results in the linker selecting an address for each instance, which
 * never changes from reference to reference, is unique for each instance and the same regardless of
 * executable (assuming the executable is mapped into the same memory location for each process).
 * In order to remove the pointer-ness of the static memory location, the address is cast to a
 * pointer difference type which is an integral type by subtracting the zero pointer from it.
 *
 */
typedef ptrdiff_t MessageId;
// typedef std::vector<void *>::difference_type MessageId;

/**
 * @brief Enumeration <b>MessageType</b> declares the global message types.
 *
 * Currently warning and doomed (error) message types are defined.  Additional type may be added
 * after MSG_DOOMED.  The MSG_SYMMETRIC bit indicates that the message was generated identically
 * across all processors.
 *
 */
enum MessageType {
  MSG_WARNING           = 0,                    ///< Message is a warning
  MSG_DOOMED            = 1,                    ///< Message is a fatal error
  MSG_EXCEPTION         = 2,                    ///< Message is an exception
  MSG_INFORMATION       = 3,                    ///< Message is informational

  MSG_TYPE_MASK         = 0x0FFFFFFF,           ///< Mask of levels
  MSG_SYMMETRIC         = 0x80000000,           ///< Message is symmetrical
  MSG_DEFERRED          = 0x40000000,           ///< Message is deferred
  MSG_UNUSED0           = 0x20000000,
  MSG_UNUSED1           = 0x10000000
};

/**
 * @brief Enumeration <b>ThrottleGroup</b> lists defined throttling groups.
 *
 * When messages are throttled, the throttle count may be reset at varior points during an
 * application run.  Some throttles defined for the application, while other may be reset at each
 * time step or other interval.  This allows warnings to be repeated at each time step rather than
 * cut off.
 *
 */
enum ThrottleGroup {
  MSG_APPLICATION       = 0,
  MSG_TIME_STEP         = 1
};

/**
 * @brief Class <b>Throttle</b> describes the cutoff limits for a message throttle.
 *
 */
struct Throttle
{
  /**
   * Creates a new <b>Throttle</b> instance.
   *
   * @param cutoff		a <b>size_t</b> value to display before the message is no longer
   *                            displayed. 
   *
   * @param group		an <b>int</b> value to identify the throttle group that this message
   *                            belongs to. 
   *
   */
  Throttle(size_t cutoff, int group)
    : m_cutoff(cutoff),
      m_group(group),
      m_count(0)
  {}
  
  size_t        m_cutoff;                       ///< Maximum number to display
  int           m_group;                        ///< Throttle group of message
  size_t        m_count;                        ///< Number which have been displayed
};

/**
 * @brief Class <b>MessageCode</b> declares a message identifier and throttle characteristics for a
 * message.  THESE MUST BE DECLARED STATIC.
 *
 * All messages have an associated message code.  This message code is used to identify a message
 * for throttling and aggregation.
 *
 * Message identifiers must be consist from reference to reference, unique for each instance, and
 * yet consist within each instance across multiple processors.  To meet these criteria, the message
 * identifier is implemented as a static memory location.  It must be declared by the application
 * developer as static.  This results in the linker selecting an address for each instance, which
 * never changes from reference to reference, is unique for each instance and the same regardless of
 * executable (assuming the executable is mapped into the same memory location for each process).
 * In order to remove the pointer-ness of the static memory location, the address is cast to a
 * pointer difference type which is an integral type by subtracting the zero pointer from it.
 *
 */
struct MessageCode
{
  /**
   * Creates a new <b>MessageCode</b> instance.
   *
   * @param throttle_cutoff	a <b>size_t</b> value to display before the message is no longer
   *                            displayed. 
   *
   * @param throttle_group	an <b>int</b> value to identify the throttle group that this message
   *                            belongs to. 
   *
   */
  MessageCode(size_t throttle_cutoff = 5, int throttle_group = MSG_APPLICATION)
    : m_id(&m_id - (MessageId *) 0),
      m_throttle(throttle_cutoff, throttle_group)
  {}

  /**
   * Creates a new <b>MessageCode</b> instance.  Be particularly careful when usnig this
   * constructor.  The message_id value must be the same on all processors for deferred message
   * reporting to work properly.  
   *
   * @param message_id		a <b>MessageId</b> value of the message id.  This value must be the
   *                            same for each message across all processors.
   *
   * @param throttle_cutoff	a <b>size_t</b> value to display before the message is no longer
   *                            displayed. 
   *
   * @param throttle_group	an <b>int</b> value to identify the throttle group that this message
   *                            belongs to. 
   *
   */
  MessageCode(MessageId message_id, size_t throttle_cutoff, int throttle_group)
    : m_id(message_id),
      m_throttle(throttle_cutoff, throttle_group)
  {}

  static MessageCode    s_defaultMessageCode;   ///< Default message code
  
  MessageId             m_id;                   ///< Message identifier
  Throttle              m_throttle;             ///< Throttle characteristics
};

/**
 * @brief Member function <b>get_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
unsigned get_message_count(unsigned message_type);

/**
 * @brief Member function <b>reset_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 */
void reset_message_count(unsigned message_type);

/**
 * @brief Member function <b>set_max_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @param max_count		an <b>unsigned int</b> ...
 *
 */
void set_max_message_count(unsigned message_type, unsigned max_count);

/**
 * @brief Member function <b>get_max_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
unsigned get_max_message_count(unsigned message_type);

/**
 * @brief Member function <b>get_message_name</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @return			a <b>std::string</b> ...
 */
const std::string &get_message_name(unsigned message_type);

/**
 * @brief Member function <b>set_message_name</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @param max_count		an <b>unsigned int</b> ...
 *
 * @param name			a <b>std::string</b> const ...
 *
 */
void register_message_type(unsigned message_type, unsigned max_count, const char *name);

/**
 * @brief Function <b>reset_message_group</b> sets the count to zero of all messages in the
 * specified throttle group.
 *
 * @param throttle_group	an <b>int</b> value of the throttle group to reset.
 *
 */
void reset_throttle_group(int throttle_group);

/**
 * @brief Member function <b>report_message</b> ...
 *
 * @param message		an <b>char</b> const pointer ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @param message_code		a <b>MessageCode</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
void report_message(const char *message, unsigned message_type, const MessageCode &message_code);

/**
 * @brief Function <b>add_deferred_message</b> adds a message to the deferred message queue.
 *
 * @param message_type		an <b>int</b> value of the message type, usually WARNING or DOOMED
 *
 * @param message_id		a <b>MessageId</b> value of the message identifier
 *
 * @param throttle_cutoff	a <b>size_t</b> value to display before the message is no longer
 *                              displayed. 
 *
 * @param throttle_group	an <b>int</b> value to identify the throttle group that this message
 *                              belongs to. 
 *
 * @param header		a <b>char</b> const pointer to the message header string.
 *
 * @param aggegrate		a <b>char</b> const pointer to the message aggregation string.
 *
 */
void add_deferred_message(int message_type, MessageId message_id, size_t throttle_cutoff, int throttle_group, const char *header, const char *aggegrate);

/**
 * @brief Function <b>report_deferred_messages</b> aggregates and reports the message on the root
 * processor. 
 *
 * @param comm			a <b>ParallelMachine</b> communicator.
 *
 */
void report_deferred_messages(ParallelMachine comm);

/**
 * @brief Function <b>aggregate_messages</b> writes a message message to the output string by
 * joining the messages from each processor, in order.  Each message is separated by the specified
 * separation string.
 *
 * @param comm			a <b>ParallelMachine</b> communicator.
 *
 * @param os			a <b>std::ostream</b> reference to the output stream to receive the
 *                              aggregated message.
 *
 * @param separator		a <b>char</b> const pointer to the separation string.
 *
 */
void aggregate_messages(ParallelMachine comm, std::ostringstream &os, const char *separator = ", ");

/**
 * @brief Function <b>operator<<</b> writes the message type name to the output stream.  If the
 * symmetric bit is set, "parallel" is prefixed to the name.
 *
 * @param os		a <b>std::ostream</b> reference to the output stream.
 *
 * @param message_type	a <b>MessageType</b> const reference to write the name of.
 *
 * @return		a <b>std::ostream</b> reference to os
 */
std::ostream &operator<<(std::ostream &os, const MessageType &message_type);

///
/// @}
///

} // namespace stk

#endif // STK_UTIL_ENVIRONMENT_RUNTIMEMESSAGE_HPP
