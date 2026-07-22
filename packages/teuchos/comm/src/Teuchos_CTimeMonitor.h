// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
