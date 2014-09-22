/*--------------------------------------------------------------------*/
/*    Copyright (c) 2013, Sandia Corporation.
/*    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*    the U.S. Governement retains certain rights in this software.
/*    
/*    Redistribution and use in source and binary forms, with or without
/*    modification, are permitted provided that the following conditions are
/*    met:
/*    
/*        * Redistributions of source code must retain the above copyright
/*          notice, this list of conditions and the following disclaimer.
/*    
/*        * Redistributions in binary form must reproduce the above
/*          copyright notice, this list of conditions and the following
/*          disclaimer in the documentation and/or other materials provided
/*          with the distribution.
/*    
/*        * Neither the name of Sandia Corporation nor the names of its
/*          contributors may be used to endorse or promote products derived
/*          from this software without specific prior written permission.
/*    
/*    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*    
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
