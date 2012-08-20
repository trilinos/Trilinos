/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/

/*
 * xfer_service_test.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: raoldfi
 */


#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_nssi_xdr.h"
#include "Trios_nssi_debug.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_oblackholestream.hpp"

#include <xfer_service_args.h>
#include "xfer_client.h"
#include "xfer_debug.h"
#include "xfer_util.h"


#include <mpi.h>
#include <assert.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

// Prototypes for client and server codes
int xfer_server_main(nssi_rpc_transport transport, MPI_Comm server_comm);
int xfer_client_main (struct xfer_args &args, nssi_service &xfer_svc, MPI_Comm client_comm);


/* -------------- private methods -------------------*/


int print_args(
        std::ostream &out,
        const struct xfer_args &args,
        const char *prefix)
{

    if (args.client_flag && args.server_flag)
        out << prefix << " ------------  ARGUMENTS (client and server) ----------- " << std::endl;
    else if (args.client_flag && !args.server_flag)
        out << prefix << " ------------  ARGUMENTS (client) ----------- " << std::endl;
    else if (!args.client_flag && args.server_flag)
        out << prefix << " ------------  ARGUMENTS (server) ----------- " << std::endl;

    out << prefix << " \tserver-url       = " << args.server_url.c_str() << std::endl;

    if (args.client_flag) {
        out << prefix << " \ttransport        = " << args.transport_name << std::endl;
        out << prefix << " \tio-method        = " << args.io_method_name << std::endl;
        out << prefix << " \tnum-trials       = " << args.num_trials << std::endl;
        out << prefix << " \tnum-reqs         = " << args.num_reqs << std::endl;
        out << prefix << " \tlen              = " << args.len << std::endl;
        out << prefix << " \tvalidate         = " << ((args.validate_flag)?"true":"false") << std::endl;
        out << prefix << " \tresult-file      = " <<
                (args.result_file.empty()?"<stdout>":args.result_file) << std::endl;
        out << prefix << " \tresult-file-mode = " << args.result_file_mode << std::endl;
    }
    out << prefix << " \tdebug            = " << args.debug_level << std::endl;
    out << prefix << " \tlogfile          = " << args.logfile << std::endl;
    out << prefix << " ----------------------------------- " << std::endl;

    return 0;
}





int main(int argc, char *argv[])
{
    int np=1, rank=0;
    int splitrank, splitsize;
    int rc = 0;
    nssi_service xfer_svc;

    int server_index=0;
    int rank_in_server=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    MPI_Barrier(MPI_COMM_WORLD);

    Teuchos::oblackholestream blackhole;
    std::ostream &out = ( rank == 0 ? std::cout : blackhole );

    struct xfer_args args;

    const int num_io_methods = 8;
    const int io_method_vals[] = {
            XFER_WRITE_ENCODE_SYNC, XFER_WRITE_ENCODE_ASYNC,
            XFER_WRITE_RDMA_SYNC, XFER_WRITE_RDMA_ASYNC,
            XFER_READ_ENCODE_SYNC, XFER_READ_ENCODE_ASYNC,
            XFER_READ_RDMA_SYNC, XFER_READ_RDMA_ASYNC};
    const char * io_method_names[] = {
            "write-encode-sync", "write-encode-async",
            "write-rdma-sync", "write-rdma-async",
            "read-encode-sync", "read-encode-async",
            "read-rdma-sync", "read-rdma-async"};

    const int num_nssi_transports = 4;
    const int nssi_transport_vals[] = {
            NSSI_RPC_PTL,
            NSSI_RPC_IB,
            NSSI_RPC_GEMINI,
            NSSI_RPC_MPI};
    const char * nssi_transport_names[] = {
            "ptl",
            "ib",
            "gni",
            "mpi"
    };


    // Initialize arguments
    args.transport=NSSI_DEFAULT_TRANSPORT;
    args.len = 1;
    args.delay = 1;
    args.io_method = XFER_WRITE_RDMA_SYNC;
    args.debug_level = LOG_WARN;
    args.num_trials = 1;
    args.num_reqs = 1;
    args.result_file_mode = "a";
    args.result_file = "";
    args.url_file = "";
    args.logfile = "";
    args.client_flag = true;
    args.server_flag = true;
    args.num_servers = 1;
    args.num_threads = 0;
    args.timeout = 500;
    args.num_retries = 5;
    args.validate_flag = true;
    args.block_distribution = true;


    bool success = true;

    /**
     * We make extensive use of the \ref Teuchos::CommandLineProcessor for command-line
     * options to control the behavior of the test code.   To evaluate performance,
     * the "num-trials", "num-reqs", and "len" options control the amount of data transferred
     * between client and server.  The "io-method" selects the type of data transfer.  The
     * server-url specifies the URL of the server.  If running as a server, the server-url
     * provides a recommended URL when initializing the network transport.
     */
    try {

        //out << Teuchos::Teuchos_Version() << std::endl << std::endl;

        // Creating an empty command line processor looks like:
        Teuchos::CommandLineProcessor parser;
        parser.setDocString(
                "This example program demonstrates a simple data-transfer service "
                "built using the NEtwork Scalable Service Interface (Nessie)."
        );

        /* To set and option, it must be given a name and default value.  Additionally,
           each option can be given a help std::string.  Although it is not necessary, a help
           std::string aids a users comprehension of the acceptable command line arguments.
           Some examples of setting command line options are:
         */

        parser.setOption("delay", &args.delay, "time(s) for client to wait for server to start" );
        parser.setOption("timeout", &args.timeout, "time(ms) to wait for server to respond" );
        parser.setOption("server", "no-server", &args.server_flag, "Run the server" );
        parser.setOption("client", "no-client", &args.client_flag, "Run the client");
        parser.setOption("len", &args.len, "The number of structures in an input buffer");
        parser.setOption("debug",(int*)(&args.debug_level), "Debug level");
        parser.setOption("logfile", &args.logfile, "log file");
        parser.setOption("num-trials", &args.num_trials, "Number of trials (experiments)");
        parser.setOption("num-reqs", &args.num_reqs, "Number of reqs/trial");
        parser.setOption("result-file", &args.result_file, "Where to store results");
        parser.setOption("result-file-mode", &args.result_file_mode, "Write mode for the result");
        parser.setOption("server-url-file", &args.url_file, "File that has URL client uses to find server");
        parser.setOption("validate", "no-validate", &args.validate_flag, "Validate the data");
        parser.setOption("num-servers", &args.num_servers, "Number of server processes");
        parser.setOption("num-threads", &args.num_threads, "Number of threads used by each server process");
        parser.setOption("block-distribution", "rr-distribution", &args.block_distribution,
                "Use a block distribution scheme to assign clients to servers");

        // Set an enumeration command line option for the io_method
        parser.setOption("io-method", &args.io_method, num_io_methods, io_method_vals, io_method_names,
                "I/O Methods for the example: \n"
                "\t\t\twrite-encode-sync : Write data through the RPC args, synchronous\n"
                "\t\t\twrite-encode-async: Write data through the RPC args - asynchronous\n"
                "\t\t\twrite-rdma-sync : Write data using RDMA (server pulls) - synchronous\n"
                "\t\t\twrite-rdma-async: Write data using RDMA (server pulls) - asynchronous\n"
                "\t\t\tread-encode-sync : Read data through the RPC result - synchronous\n"
                "\t\t\tread-encode-async: Read data through the RPC result - asynchronous\n"
                "\t\t\tread-rdma-sync : Read data using RDMA (server puts) - synchronous\n"
                "\t\t\tread-rdma-async: Read data using RDMA (server puts) - asynchronous");


        // Set an enumeration command line option for the io_method
        parser.setOption("transport", &args.transport, num_nssi_transports, nssi_transport_vals, nssi_transport_names,
                "NSSI transports (not all are available on every platform): \n"
                "\t\t\tportals : Cray or Schutt\n"
                "\t\t\tinfiniband : libibverbs\n"
                "\t\t\tgemini : Cray\n"
                "\t\t\tmpi : isend/irecv implementation\n"
                );



        /* There are also two methods that control the behavior of the
           command line processor.  First, for the command line processor to
           allow an unrecognized a command line option to be ignored (and
           only have a warning printed), use:
         */
        parser.recogniseAllOptions(true);

        /* Second, by default, if the parser finds a command line option it
           doesn't recognize or finds the --help option, it will throw an
           std::exception.  If you want prevent a command line processor from
           throwing an std::exception (which is important in this program since
           we don't have an try/catch around this) when it encounters a
           unrecognized option or help is printed, use:
         */
        parser.throwExceptions(false);

        /* We now parse the command line where argc and argv are passed to
           the parse method.  Note that since we have turned off std::exception
           throwing above we had better grab the return argument so that
           we can see what happened and act accordingly.
         */
        Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn= parser.parse( argc, argv );

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
            return 0;
        }

        if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
            return 1; // Error!

        }

        // Here is where you would use these command line arguments but for this example program
        // we will just print the help message with the new values of the command-line arguments.
        //if (rank == 0)
        //    out << "\nPrinting help message with new values of command-line arguments ...\n\n";

        //parser.printHelpMessage(argv[0],out);

    }

    TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

    log_debug(args.debug_level, "%d: Finished processing arguments", rank);


    if (!success) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!args.server_flag && args.client_flag) {
        /* initialize logger */
        if (args.logfile.empty()) {
            logger_init(args.debug_level, NULL);
        } else {
            char fn[1024];
            sprintf(fn, "%s.client.%03d.log", args.logfile.c_str(), rank);
            logger_init(args.debug_level, fn);
        }
    } else if (args.server_flag && !args.client_flag) {
        /* initialize logger */
        if (args.logfile.empty()) {
            logger_init(args.debug_level, NULL);
        } else {
            char fn[1024];
            sprintf(fn, "%s.server.%03d.log", args.logfile.c_str(), rank);
            logger_init(args.debug_level, fn);
        }
    } else if (args.server_flag && args.client_flag) {
        /* initialize logger */
        if (args.logfile.empty()) {
            logger_init(args.debug_level, NULL);
        } else {
            char fn[1024];
            sprintf(fn, "%s.%03d.log", args.logfile.c_str(), rank);
            logger_init(args.debug_level, fn);
        }
    }

    log_level debug_level = args.debug_level;

    // Communicator used for both client and server (may split if using client and server)
    MPI_Comm comm;

    log_debug(debug_level, "%d: Starting xfer-service test", rank);

    /**
     * Since this test can be run as a server, client, or both, we need to play some fancy
     * MPI games to get the communicators working correctly.  If we're executing as both
     * a client and a server, we split the communicator so that the client thinks its
     * running by itself.
     */
    int color = 0;  // color=0-->server, color=1-->client
    if (args.client_flag && args.server_flag) {
        if (np < 2) {
            log_error(debug_level, "Must use at least 2 MPI processes for client and server mode");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        // Split the communicators. Put all the servers as the first ranks.
        if (rank < args.num_servers) {
            color = 0;
            log_debug(debug_level, "rank=%d is a server", rank);
        }
        else {
            color = 1;  // all others are clients
            log_debug(debug_level, "rank=%d is a client", rank);
        }

        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    }
    else {
        if (args.client_flag)
            color=1;
        else if (args.server_flag)
            color=0;
        else {
            log_error(debug_level, "Must be either a client or a server");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    }

    MPI_Comm_rank(comm, &splitrank);
    MPI_Comm_size(comm, &splitsize);

    log_debug(debug_level, "%d: Finished splitting communicators", rank);

    /**
     * Initialize the Nessie interface by specifying a transport, encoding scheme, and a
     * recommended URL.  \ref NSSI_DEFAULT_TRANSPORT is usually the best choice, since it
     * is often the case that only one type of transport exists on a particular platform.
     * Currently supported transports are \ref NSSI_RPC_PTL, \ref NSSI_RPC_GNI, and
     * \ref NSSI_RPC_IB.  We only support one type of encoding scheme so NSSI_DEFAULT_ENCODE
     * should always be used for the second argument.   The URL can be specified (as we did for
     * the server, or NULL (as we did for the client).  This is a recommended value.  Use the
     * \ref nssi_get_url function to find the actual value.
     */
    nssi_rpc_init((nssi_rpc_transport)args.transport, NSSI_DEFAULT_ENCODE, NULL);

    // Get the Server URL
    std::string my_url(NSSI_URL_LEN, '\0');
    nssi_get_url((nssi_rpc_transport)args.transport, &my_url[0], NSSI_URL_LEN);

    // If running as both client and server, gather and distribute
    // the server URLs to all the clients.
    if (args.server_flag && args.client_flag) {

        std::string all_urls;

        // This needs to be a vector of chars, not a string
        all_urls.resize(args.num_servers * NSSI_URL_LEN, '\0');

        // Have servers gather their URLs
        if (color == 0) {
            assert(args.num_servers == splitsize);  // these should be equal

            log_debug(debug_level, "%d: Gathering urls: my_url=%s", rank, my_url.c_str());

            // gather all urls to rank 0 of the server comm (also rank 0 of MPI_COMM_WORLD)
            MPI_Gather(&my_url[0], NSSI_URL_LEN, MPI_CHAR,
                    &all_urls[0], NSSI_URL_LEN, MPI_CHAR, 0, comm);
        }

        // broadcast the full set of server urls to all processes
        MPI_Bcast(&all_urls[0], all_urls.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

        log_debug(debug_level, "%d: Bcast urls, urls.size=%d", rank, all_urls.size());

        if (color == 1) {

            // For block distribution scheme use the utility function (in xfer_util.cpp)
            if (args.block_distribution) {
                // Use this utility function to calculate the server_index
                xfer_block_partition(args.num_servers, splitsize, splitrank, &server_index, &rank_in_server);
            }

            // Use a simple round robin distribution scheme
            else {
                server_index   = splitrank % args.num_servers;
                rank_in_server = splitrank / args.num_servers;
            }

            // Copy the server url out of the list of urls
            int offset = server_index * NSSI_URL_LEN;

            args.server_url = all_urls.substr(offset, NSSI_URL_LEN);

            log_debug(debug_level, "client %d assigned to server \"%s\"", splitrank, args.server_url.c_str());
        }


        log_debug(debug_level, "%d: Finished distributing server urls, server_url=%s", rank, args.server_url.c_str());
    }

    // If running as a client only, have to get the list of servers from the urlfile.
    else if (!args.server_flag && args.client_flag){

        sleep(args.delay);  // give server time to get started

        std::vector< std::string > urlbuf;
        xfer_read_server_url_file(args.url_file.c_str(), urlbuf, comm);
        args.num_servers = urlbuf.size();

        // For block distribution scheme use the utility function (in xfer_util.cpp)
        if (args.block_distribution) {
            // Use this utility function to calculate the server_index
            xfer_block_partition(args.num_servers, splitsize, splitrank, &server_index, &rank_in_server);
        }

        // Use a simple round robin distribution scheme
        else {
            server_index   = splitrank % args.num_servers;
            rank_in_server = splitrank / args.num_servers;
        }

        args.server_url = urlbuf[server_index];
        log_debug(debug_level, "client %d assigned to server \"%s\"", splitrank, args.server_url.c_str());
    }

    else if (args.server_flag && !args.client_flag) {
        args.server_url = my_url;

        if (args.url_file.empty()) {
            log_error(debug_level, "Must set --url-file");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        xfer_write_server_url_file(args.url_file.c_str(), my_url.c_str(), comm);
    }

    // Set the debug level for the xfer service.
    xfer_debug_level = args.debug_level;

    // Print the arguments after they've all been set.
    args.io_method_name = std::string(io_method_names[args.io_method]);
    args.transport_name = std::string(nssi_transport_names[args.transport]);

    log_debug(debug_level, "%d: server_url=%s", rank, args.server_url.c_str());

    print_args(out, args, "%");

    log_debug(debug_level, "server_url=%s", args.server_url.c_str());

    //------------------------------------------------------------------------------
    /** If we're running this job with a server, the server always executes on node 0.
     *  In this example, the server is a single process.
     */
    if (color == 0) {
        rc = xfer_server_main((nssi_rpc_transport)args.transport, comm);
        log_debug(debug_level, "Server is finished");
    }

    // ------------------------------------------------------------------------------
     /**  The parallel client will execute this branch.  The root node, node 0, of the client connects
      *   connects with the server, using the \ref nssi_get_service function.  Then the root
      *   broadcasts the service description to the other clients before starting the main
      *   loop of the client code by calling \ref xfer_client_main.
      */
    else {
        int i;
        int client_rank;

        // get rank within the client communicator
        MPI_Comm_rank(comm, &client_rank);

        nssi_init((nssi_rpc_transport)args.transport);

        // Only one process needs to connect to the service
        // TODO: Make get_service a collective call (some transports do not need a connection)
        //if (client_rank == 0) {
        {


            // connect to remote server
            for (i=0; i < args.num_retries; i++) {
                log_debug(debug_level, "Try to connect to server: attempt #%d, url=%s", i, args.server_url.c_str());
                rc=nssi_get_service((nssi_rpc_transport)args.transport, args.server_url.c_str(), args.timeout, &xfer_svc);
                if (rc == NSSI_OK)
                    break;
                else if (rc != NSSI_ETIMEDOUT) {
                    log_error(xfer_debug_level, "could not get svc description: %s",
                            nssi_err_str(rc));
                    break;
                }
            }
        }

        // wait for all the clients to connect
        MPI_Barrier(comm);

        //MPI_Bcast(&rc, 1, MPI_INT, 0, comm);

        if (rc == NSSI_OK) {
            if (client_rank == 0) log_debug(debug_level, "Connected to service on attempt %d\n", i);

            // Broadcast the service description to the other clients
            //log_debug(xfer_debug_level, "Bcasting svc to other clients");
            //MPI_Bcast(&xfer_svc, sizeof(nssi_service), MPI_BYTE, 0, comm);

            log_debug(debug_level, "Starting client main");
            // Start the client code
            xfer_client_main(args, xfer_svc, comm);


            MPI_Barrier(comm);

            // Tell one of the clients to kill the server
            if (rank_in_server == 0) {
                log_debug(debug_level, "%d: Halting xfer service", rank);
                rc = nssi_kill(&xfer_svc, 0, 5000);
            }
        }

        else {
            if (client_rank == 0)
                log_error(debug_level, "Failed to connect to service after %d attempts: ABORTING", i);
            success = false;
            //MPI_Abort(MPI_COMM_WORLD, -1);
        }

        nssi_fini((nssi_rpc_transport)args.transport);

    }

    log_debug(debug_level, "%d: clean up nssi", rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // Clean up nssi_rpc
    rc = nssi_rpc_fini((nssi_rpc_transport)args.transport);
    if (rc != NSSI_OK)
        log_error(debug_level, "Error in nssi_rpc_fini");

    log_debug(debug_level, "%d: MPI_Finalize()", rank);
    MPI_Finalize();

    logger_fini();

    if(success && (rc == NSSI_OK))
      out << "\nEnd Result: TEST PASSED" << std::endl;
    else
        out << "\nEnd Result: TEST FAILED" << std::endl;

    return ((success && (rc==NSSI_OK)) ? 0 : 1 );
}
