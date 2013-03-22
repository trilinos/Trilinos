/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
/*-------------------------------------------------------------------------*/
/**  @file trios_signal.c
 *
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 560 $
 *   $Date: 2006-02-28 14:02:02 -0700 (Tue, 28 Feb 2006) $
 *
 */

#include "Trios_config.h"

#include <stdlib.h>
#include <signal.h>

#include "Trios_logger.h"
#include "Trios_signal.h"
#include "Trios_threads.h"

#include <utility>
#include <list>

using namespace std;


/* local variables */
static int volatile _exit_now = 0;

/* --------------------- Private methods ------------------- */

static void sighandler(int sig)
{
    log_warn(LOG_UNDEFINED, "Caught signal %d, setting exit_now flag", sig);
    trios_abort();
    //exit(sig);
}


/**
 * @brief Install signal handlers.
 */
int trios_install_sighandler()
{
    struct sigaction new_action, old_action;

    new_action.sa_handler = sighandler;
    sigemptyset (&new_action.sa_mask);
    new_action.sa_flags = 0;

    sigaction (SIGINT, NULL, &old_action);
    if (old_action.sa_handler != SIG_IGN) {
        sigaction (SIGINT, &new_action, NULL);
    }
    sigaction (SIGHUP, NULL, &old_action);
    if (old_action.sa_handler != SIG_IGN) {
        sigaction (SIGHUP, &new_action, NULL);
    }
    sigaction (SIGTERM, NULL, &old_action);
    if (old_action.sa_handler != SIG_IGN) {
        sigaction (SIGTERM, &new_action, NULL);
    }
    sigaction (SIGABRT, NULL, &old_action);
    if (old_action.sa_handler != SIG_IGN) {
        sigaction (SIGABRT, &new_action, NULL);
    }
    return 0;
}


/**
 * @brief Return the exit_now variable.  If set,
 * it is time to exit the service.
 */
int trios_exit_now() {
    return _exit_now;
}


/**
 * @brief Cleanly abort the running service.
 *
 * The abort function kills a running service by sending a
 * SIGINT signal to the running process.  If the service
 * has already started,  the signal is caught by the
 * signal handler.
 */
void trios_abort()
{
    log_level debug_level = LOG_UNDEFINED;

    /* kill(0,SIGINT) */
    log_debug(debug_level, "Received abort()... setting exit_now flag");
    _exit_now = 1;
}
