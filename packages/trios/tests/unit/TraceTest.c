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
#include "Trios_trace.h"

static int gid;

static void testFunction( int n )
{
  	int	j;
	int pid = 0;

  	trace_event(gid, 3, pid, "entered testFunction");

  	for (j = 0; j < n; j++) {
    	   trace_event(gid, 2, pid, "inside for(j) loop");
  	}

  	trace_event(gid, 4, pid, "exiting testFunction");
}

int main(int argc, char *argv[])
{
	int i;
	int end = 3;
	int event_id = 0;
	int subevent_id = 1;
	int thread_id = 0;
	char fname[256];


	/* assign a file name for the trace file */
	sprintf(fname, "%s.sddf", argv[0]);
	trace_init(fname, TRACE_SDDF_BINARY);


	/* register a trace group */
	trace_register_group("test", &gid);


	/* turn on tracing for the test group */
	trace_enable_gid(gid);

	/* write a generic trace event */
	trace_event(gid, event_id, thread_id, "start");

	trace_start_interval(gid, thread_id);

	for (i = 1; i <= end; i++) {
		trace_start_interval(gid, thread_id);
		trace_inc_count(gid, 2, thread_id, "loop counter");

		trace_event(gid, 1, thread_id, "top of for(i) loop");
		testFunction(i);
		/*sleep(end-i);*/
		sleep(1);
		trace_end_interval(gid, subevent_id, thread_id, "internal loop interval");
	}

	trace_end_interval(gid, event_id, thread_id, "outside loop interval");

	trace_put_all_counts(gid, thread_id, "put all counts");

	trace_fini();
	return 0;
}

