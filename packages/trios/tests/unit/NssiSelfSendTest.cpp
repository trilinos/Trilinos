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
/*
 * SelfSendTest.c
 *
 *  Created on: March 20, 2014
 *      Author: thkorde
 */

#include "Trios_nssi_client.h"
#include "Trios_nssi_server.h"
#include "Trios_nssi_fprint_types.h"

#include "Trios_logger.h"
#include "Trios_timer.h"

#include <unistd.h>
#include <pthread.h>

#include <iostream>

log_level selfsend_debug_level = LOG_UNDEFINED;

std::string  my_url(NSSI_URL_LEN, '\0');

pthread_barrier_t barrier;

bool success=true;

static void *do_wait(void *args)
{
    nssi_service svc;

    // service init
    nssi_service_init(NSSI_DEFAULT_TRANSPORT, NSSI_SHORT_REQUEST_SIZE, &svc);
    nssi_get_url(NSSI_DEFAULT_TRANSPORT, &my_url[0], NSSI_URL_LEN);

    // client is waiting for us to initialize
    pthread_barrier_wait(&barrier);

    svc.max_reqs = -1;
    int rc = nssi_service_start(&svc);
    if (rc != NSSI_OK) {
        log_info(selfsend_debug_level, "exited selfsend_svc: %s",
                nssi_err_str(rc));
    }

    /* finalize the service */
    nssi_service_fini(&svc);

    return(NULL);
}

static pthread_t wait_thread;
static void launch_wait_thread()
{
    /* Start up polling thread */
    fprintf (stdout, "Start wait thread.\n");
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int rc = pthread_create(&wait_thread, &attr, do_wait, NULL);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot create polling thread, err code = %d\n", rc);
        return;
    }
    pthread_attr_destroy(&attr);
}
static void join_wait_thread()
{
    void *status;
    int rc = pthread_join(wait_thread, &status);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot join wait thread, err code = %d\n", rc);
    } else {
        fprintf (stdout, "Joined wait thread.\n");
    }
}


int main(int argc, char *argv[])
{
    int rc;
    nssi_service svc;

    logger_init(LOG_ERROR, NULL);

    pthread_barrier_init(&barrier, NULL, 2);

    // common init
    nssi_rpc_init(NSSI_DEFAULT_TRANSPORT, NSSI_DEFAULT_ENCODE, NULL);

    launch_wait_thread();

    // wait for the service to initialize
    pthread_barrier_wait(&barrier);

    // client init
    nssi_init(NSSI_DEFAULT_TRANSPORT);

    for (int i=0; i < 3; i++) {
        log_debug(selfsend_debug_level, "Try to connect to server: attempt #%d", i);
        rc=nssi_get_service(NSSI_DEFAULT_TRANSPORT, my_url.c_str(), 5000, &svc);
        if (rc == NSSI_OK)
            break;
        else if (rc != NSSI_ETIMEDOUT) {
            log_error(selfsend_debug_level, "could not get svc description: %s",
                    nssi_err_str(rc));
            break;
        }
    }

    sleep(1);

    // shutdown the service
    rc = nssi_kill(&svc, 0, 5000);
    if (rc != NSSI_OK) {
        log_error(selfsend_debug_level, "Error in nssi_kill");
        goto err;
    }

    // finalize the client
    nssi_fini(NSSI_DEFAULT_TRANSPORT);

    join_wait_thread();

    // finalize the NSSI library
    rc = nssi_rpc_fini(NSSI_DEFAULT_TRANSPORT);
    if (rc != NSSI_OK) {
        log_error(selfsend_debug_level, "Error in nssi_rpc_fini");
        goto err;
    }

    logger_fini();

    std::cout << "\nEnd Result: TEST PASSED" << std::endl;

    return 0;

err:
    std::cout << "\nEnd Result: TEST FAILED" << std::endl;

    return 1;
}
