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

#include "Trios_config.h"
#include "Trios_logger.h"
#include "Trios_nssi_client.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

//#include "nssi_types.h"

#include <string>
using namespace std;

int admin_debug_level = LOG_UNDEFINED;

/* ----------------- COMMAND-LINE OPTIONS --------------- */


static void read_contact_info(const char *fname, char *url, int maxlen)
{
    const char *contact_filename=NULL;
    FILE *cf=NULL;

    if ((fname==NULL) || (fname[0]=='\0')) {
        contact_filename=getenv("NNTI_CONTACT_FILENAME");
    } else {
        contact_filename=fname;
    }
    if (contact_filename==NULL) {
        url[0]='\0';
        return;
    }
    cf=fopen(contact_filename, "r");
    if (cf == NULL) {
        url[0]='\0';
        return;
    }
    fgets(url, maxlen, cf);
    fclose(cf);
}

int kill_svc(nssi_service *svc, int sig, long timeout)
{
    int rc;


    /* try to kill the remote service (TODO: need a timeout) */
    rc = nssi_kill(svc, sig, timeout);
    if (rc != NSSI_OK) {
        log_error(admin_debug_level, "could not kill svc: %s",
                nssi_err_str(rc));
        return rc;
    }
    return rc;
}

int
main (int argc, char *argv[])
{
    int rc;

    // command-line arguments
    int retries = 0;
    int sig = 0;
    int timeout = 1000;
    log_level debug_level = LOG_ERROR;
    string logfile("");

    nssi_service svc;
    char my_url[NSSI_URL_LEN];

    std::string server_url("");
    char        server_str[NSSI_URL_LEN];
    std::string contact_file("");   /* the file where the server's url should be written */


    try {
        Teuchos::CommandLineProcessor parser;

        // init parser
        parser.setDocString("Kill an NSSI server");

        parser.setOption("verbose", (int *)(&debug_level), "Debug level.");
        parser.setOption("logfile", &logfile, "Path to file for debug statements");
        parser.setOption("server-url", &server_url, "URL of NSSI service");
        parser.setOption("contact-file", &contact_file, "Where to read the server's URL");
        parser.setOption("timeout", &timeout, "Timout for contacting services (ms)");
        parser.setOption("retries", &retries, "Number of times to retry before exiting");
        parser.setOption("sig", &sig, "Signal to use for the kill command");

        parser.recogniseAllOptions();
        parser.throwExceptions();

        Teuchos::CommandLineProcessor::EParseCommandLineReturn
        parseReturn= parser.parse( argc, argv );
        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
            return 0;
        }
        if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
            return 1; // Error!
        }
    }
    catch (...) {
        exit(-1);
    }

    /* initialize the logger */
    logger_init(debug_level, logfile.c_str());

    if (server_url.c_str()[0]=='\0') {
        sleep(1);
        log_debug(debug_level, "reading URL from file");
        read_contact_info(contact_file.c_str(), server_str, NSSI_URL_LEN);
    } else {
        log_debug(debug_level, "using URL from command-line");
        strcpy(server_str, server_url.c_str());
    }

    nssi_rpc_init(NSSI_DEFAULT_TRANSPORT, NSSI_DEFAULT_ENCODE, NULL);

    nssi_get_url(NSSI_DEFAULT_TRANSPORT, my_url, NSSI_URL_LEN);

//    sleep(1);
    log_info(debug_level, "\nTrying to get service at %s", server_str);

    rc=nssi_get_service(NSSI_DEFAULT_TRANSPORT, server_str, timeout, &svc);
    if (rc != NSSI_OK) {
        log_error(admin_debug_level, "could not get svc description: %s",
                nssi_err_str(rc));
        return rc;
    }
    rc = kill_svc(&svc, sig, timeout);
    if (rc == NSSI_ETIMEDOUT) {
        fprintf(stderr, "Timed out trying to kill (%s)\n",
                server_url.c_str());
        return rc;
    }
    else if (rc != NSSI_OK) {
        log_error(admin_debug_level, "failed to kill service: %s",
                nssi_err_str(rc));
        return rc;
    }

    nssi_rpc_fini(NSSI_DEFAULT_TRANSPORT);

    return 0;
}
