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

#include <unistd.h>
#include <stdint.h>
#include "Trios_timer.h"

int main(int argc, char *argv[])
{
    int sleeptime = 3;


    uint64_t start_ns=0, end_ns=0;
    uint64_t start_us, end_us;
    uint64_t start_ms, end_ms;
    double start_sec, end_sec;

    printf("Test trios_timer (This is not a resolution test)\n");

    /* output the timer implementation we are using */
    printf("Timer implementation = %s\n\n", trios_timer_getimpl());

    trios_get_time_ns();
    sleep(1);


    start_ns = trios_get_time_ns();
    start_us = trios_get_time_us();
    start_ms = trios_get_time_ms();
    start_sec = trios_get_time();

    /* sleep three seconds */
    sleep(sleeptime);

    end_ns = trios_get_time_ns();
    end_us = trios_get_time_us();
    end_ms = trios_get_time_ms();
    end_sec = trios_get_time();

    printf("slept for %d seconds:\n", sleeptime);
    printf("\tns = %lu\n", (unsigned long)(end_ns-start_ns));
    printf("\tus = %lu\n", (unsigned long)(end_us-start_us));
    printf("\tms = %lu\n", (unsigned long)(end_ms-start_ms));
    printf("\tsec = %f\n", end_sec-start_sec);


    fprintf(stdout, "\nCalling trios_timer_test() ... ");

    if (trios_timer_test() != 0) {
        fprintf(stdout, "failed\n");
        return -1;
    }
    else {
        fprintf(stdout, "passed\n");
    }

    return 0;
}
