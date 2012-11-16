/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_Signal_hpp
#define STK_UTIL_DIAG_Signal_hpp

#include <string>

#include <setjmp.h>

/**
 * @file
 *
 * Signal handling is enabled via the <b>activate_signals()</b> function.  The
 * following signal behaviors are actived:
 *
 *	HUP		Sets the shutdown_request flag.
 *	TERM		Long jumps with Termiante signal received message.
 *	SEGV		Long jumps with Segmentation violation error message.
 *	BUS		Long jumps with Bus error message.
 *	ILL		Long jumps with Illegal instruction error message.
 *
 *
 */

namespace sierra {
namespace Env {

///
/// @addtogroup EnvDetail
/// @{
///

/**
 * @brief Function <b>activate_signals</b> enables the signal handlers.
 *
 */
void activate_signals();

/**
 * @brief Function <b>deactivate_signals</b> disables the signal handlers.
 *
 */
void deactivate_signals();

/**
 * @brief Function <b>get_sigjmpbuf</b> enables signal handling and returns a
 * pointer to the jump buffer for <b>::sigsetjmp</b> and
 * <b>::siglongjmp()</b>.
 *
 * if (::sigsetjmp(*sierra::Env::get_signalException(), 1))
 *   throw sierra::RuntimeError(sierra::Env::get_signal_message());
 *
 * @return			a <b>sigjmp_buf</b> pointer to the jmp buffer.
 */
sigjmp_buf *get_sigjmpbuf();

// /**
//  * @brief Function <b>disable_siglongjmp</b> disables the long jump buffer.  When
//  * signals are received, they return to the caller without long jumping to the set jump point.
//  *
//  */
// void disable_siglongjmp();

/**
 * @brief Function <b>get_signal_message</b> returns the message associated with the
 * most recent signal.
 *
 * @return			a <b>std::string</b> const reference to the most
 *				recent signal message.
 */
const std::string &get_signal_message();

/**
 * @brief Function <b>request_shutdown</b> sets the shutdown requested flag so that
 * future calls to <b>shutdown_requested()</b> return true;
 *
 */
bool HUP_received();

/**
 * @brief Function <b>shutdown_requested</b> returns true if an application shutdown
 * has requested via the <b>request_shutdown</b> has been called.
 *
 * @return			a <b>bool</b> value of true if application has been
 *				requested to shutdown.
 */
bool shutdown_requested();

///
/// @}
///

} // namespace Env
} // namespace sierra

#endif // STK_UTIL_DIAG_Signal_hpp
