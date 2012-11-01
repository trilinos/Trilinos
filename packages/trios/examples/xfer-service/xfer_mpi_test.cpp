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
 * xfer_service_test.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: raoldfi
 */


#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_nssi_xdr.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_oblackholestream.hpp"

#include <xfer_service_args.h>
#include "xfer_client.h"
#include "xfer_debug.h"


#include <mpi.h>
#include <assert.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <fstream>


// Prototypes for client and server codes
int xfer_mpi_server_main(MPI_Comm server_comm);
int xfer_mpi_client_main (struct xfer_args &args, const int server_rank, MPI_Comm client_comm);


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

    out << prefix << " \tserver-url       = " << args.server_url << std::endl;

    if (args.client_flag) {
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

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    MPI_Barrier(MPI_COMM_WORLD);

    Teuchos::oblackholestream blackhole;
    std::ostream &out = ( rank == 0 ? std::cout : blackhole );

    struct xfer_args args;

    const int num_io_methods = 6;
    const int io_method_vals[] = {
            XFER_MPI_SEND, XFER_MPI_ISEND,
            XFER_MPI_RECV, XFER_MPI_IRECV,
            XFER_MPI_PUT, XFER_MPI_GET};
    const char * io_method_names[] = {
            "send", "isend",
            "recv", "irecv",
            "put", "get"};


    // Initialize arguments
    args.len = 1;
    args.delay = 1;
    args.io_method = XFER_MPI_SEND;
    args.debug_level = LOG_WARN;
    args.num_trials = 1;
    args.num_reqs = 1;
    args.result_file_mode = "a";
    args.result_file = "";
    args.url_file = "";
    args.logfile = "";
    args.client_flag = true;
    args.server_flag = true;
    args.timeout = 500;
    args.num_retries = 5;
    args.validate_flag = true;
    args.server_url = "";

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
        parser.setOption("len", &args.len, "The number of structures in an input buffer");
        parser.setOption("debug",(int*)(&args.debug_level), "Debug level");
        parser.setOption("logfile", &args.logfile, "log file");
        parser.setOption("num-trials", &args.num_trials, "Number of trials (experiments)");
        parser.setOption("num-reqs", &args.num_reqs, "Number of reqs/trial");
        parser.setOption("result-file", &args.result_file, "Where to store results");
        parser.setOption("result-file-mode", &args.result_file_mode, "Write mode for the result");
        parser.setOption("validate", "no-validate", &args.validate_flag, "Validate the data");

        // Set an enumeration command line option for the io_method

        parser.setOption("io-method", &args.io_method, num_io_methods, io_method_vals, io_method_names,
                "I/O Methods for the example: \n"
                "\t\t\tsend : Write data using the MPI_Send function\n"
                "\t\t\tisend: Write data using the MPI_Isend function \n"
                "\t\t\trecv : Read data using MPI_Recv function\n"
                "\t\t\tirecv : Read data using the MPI_Irecv function\n"
                "\t\t\tput: Write data using the MPI_Put function\n"
                "\t\t\tget: Read data through the MPI_Get function");




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
    debug_level = xfer_debug_level;

    // Communicator used for both client and server (may split if using client and server)
    MPI_Comm comm;

    log_debug(debug_level, "%d: Starting xfer-mpi test", rank);


    /**
     * We need to play some fancy MPI games to get the communicators working correctly.
     * We're running as both a client and a server in the same app, so we need to
     * split the communicator so that the client can barrier sync without involving the server.
     */
    if (np < 2) {
        log_error(debug_level, "Must use at least 2 MPI processes for client and server mode");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // Split the communicators. Processors with color=0 are servers.
    int server_rank = 0;
    int color = (rank == server_rank) ? 0 : 1; // only one server
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);

    MPI_Comm_rank(comm, &splitrank);
    MPI_Comm_size(comm, &splitsize);

//    std::cout << "rank=" << rank << "/" << np << ", color=" << color <<
//                ", new_rank=" << splitrank << "/" << splitsize << std::endl << std::endl;



    // Set the debug level for the xfer service.
    xfer_debug_level = args.debug_level;

    // Print the arguments after they've all been set.
    args.io_method_name = io_method_names[args.io_method];
    print_args(out, args, "%");


    //------------------------------------------------------------------------------
    /** If we're running this job with a server, the server always executes on node 0.
     *  In this example, the server is a single process.
     */
    if (rank == server_rank) {
        rc = xfer_mpi_server_main(comm);
        log_debug(debug_level, "Server is finished, rc=%d", rc);
    }

    // ------------------------------------------------------------------------------
     /**  The parallel client will execute this branch.  */
    else {
        int i;
        int client_rank;

            log_debug(debug_level, "Starting client main");
            // Start the client code
            rc = xfer_mpi_client_main(args, server_rank, comm);
            log_debug(debug_level, "Client is finished, rc=%d", rc);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    log_debug(debug_level, "%d: MPI_Finalize()", rank);
    MPI_Finalize();


    if (!success) {
        log_warn(debug_level, "%d: success=FALSE", rank);
    }
    if (rc != 0) {
        log_warn(debug_level, "%d: rc=%d", rank, rc);
    }

    if(success && (rc == 0))
      out << "\nEnd Result: TEST PASSED" << std::endl;
    else
        out << "\nEnd Result: TEST FAILED" << std::endl;

    return ((success && (rc==NSSI_OK)) ? 0 : 1 );
}
