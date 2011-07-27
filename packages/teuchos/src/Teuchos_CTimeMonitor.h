/* @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#ifndef TEUCHOS_CTIMEMONITOR_H
#define TEUCHOS_CTIMEMONITOR_H

/*! \file Teuchos_CTimeMonitor.hpp
    \brief Timer functions for C that starts and stops timers for C code.
*/

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Start a timer with a given name and ID.
 *
 * \param  timerName
 *           [in] Globally unique null-terminated string name of the timer.
 *           This is only significant on the first call.
 *
 * \param  timerID
 *           [in] On first call, <tt>timerID</tt> should be less than
 *           <tt>0</tt> On future calls, it should be what was returned by the
 *           first call.
 *
 * \returns on first call <tt>returnVal</tt> gives the ID of a newly created
 * timer of the given globally unique name <tt>timerName</tt>.  On future
 * calls, <tt>returnVal==timerID</tt>.
 *
 * \note You can not start the same timer more than once.  You must stop a
 * timer with <tt>Teuchos_stopTimer()</tt> before you can call this function
 * to start it again.
 */
int Teuchos_startTimer( char timerName[], int timerID );

/** \brief Stop a timer that was started with <tt>Teuchos_startTimer()</tt>.
 *
 * \param  timerID
 *           [in] Must be the ID returned from a prior call to
 *           <tt>Teuchos_startTimer()</tt>.
 *
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>timerID >= 0</tt> and it must have been
 *   created by a prior call to <tt>Teuchos_startTimer()</tt>.
 * </ul>
 *
 * \note It is okay to stop a timer more than once (i.e. stop a timer that is
 * not running).  But, the timer must actually exist.
 */
void Teuchos_stopTimer( int timerID );

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TEUCHOS_CTIMEMONITOR_H */
