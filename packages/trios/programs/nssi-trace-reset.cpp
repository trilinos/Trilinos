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

#include "nssi_config.h"

#include <time.h>

#include <stdlib.h>
#include <string.h>
#include <string>
using namespace std;

#include "nssi_types.h"
#include "nssi_client.h"
#include "nssi_logger.h"
#include "nssi_CommandLineProcessor.hpp"
using namespace Nessie;

log_level debug_level = LOG_UNDEFINED;

enum TraceFtype {BINARY=0, ASCII=1};

int
main (int argc, char *argv[])
{
    int rc;
    nssi_service svc;

    // Command-line options
    string addr;
    int retries = 0;
    int sig = 0;
    int timeout = 1000;
    log_level debug_level = LOG_ERROR;
    std::string logfile("");
    std::string trace_enable("");
    std::string trace_path("");
    enum TraceFtype trace_ftype = BINARY;

    std::string server_url;
    char my_url[NSSI_URL_LEN];

    try {
        Nessie::CommandLineProcessor parser;

        // init parser
        parser.setDocString("Ping a remote Nessie server");

        parser.setOption("verbose", (int *)(&debug_level), "Debug level.");
        parser.setOption("logfile", &logfile, "Path to logfile");
        parser.setOption("server-url", &server_url, "Portals network ID");
        parser.setOption("timeout", &timeout, "Timout for contacting services (ms)");
        parser.setOption("retries", &retries, "Number of times to retry before exiting");
        parser.setOption("sig", &sig, "Signal to use for the kill command");
        parser.setOption("trace-enable", &trace_enable,
                "Enable tracing on specified groups (comma-separated list)");
        parser.setOption("trace-path", &trace_path, "Path to trace file");

        const int num_ftypes = 2;
        const char *ftype_opt_names[] = {"binary","ascii"};
        const enum TraceFtype ftype_opt_values[] = {BINARY, ASCII};
        parser.setOption("trace-ftype", &trace_ftype, num_ftypes,
                ftype_opt_values, ftype_opt_names, "File type of trace file");

        parser.recogniseAllOptions();
        parser.throwExceptions();

        Nessie::CommandLineProcessor::EParseCommandLineReturn
        parseReturn= parser.parse( argc, argv );
        if( parseReturn == Nessie::CommandLineProcessor::PARSE_HELP_PRINTED ) {
            return 0;
        }
        if( parseReturn != Nessie::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
            return 1; // Error!
        }
    }
    catch (...) {
        exit(-1);
    }


    /* initialize the logger */
    logger_init(debug_level, logfile.c_str());

    nssi_rpc_init(NSSI_DEFAULT_TRANSPORT, NSSI_DEFAULT_ENCODE, NULL);

    nssi_get_url(NSSI_DEFAULT_TRANSPORT, my_url, NSSI_URL_LEN);


    rc=nssi_get_service(NSSI_DEFAULT_TRANSPORT, server_url.c_str(), timeout, &svc);
    if (rc != NSSI_OK) {
        log_error(debug_level, "could not get svc description: %s",
                nssi_err_str(rc));
        return rc;
    }


    /* reset the tracing on the remote service */
    rc = nssi_trace_reset(&svc,
            trace_path.c_str(),
            trace_ftype,
            trace_enable.c_str(),
            timeout);
    if (rc != NSSI_OK) {
        log_error(debug_level, "could not reset tracing: %s",
                nssi_err_str(rc));
        return rc;
    }

    nssi_rpc_fini(NSSI_DEFAULT_TRANSPORT);

    return 0;
}
