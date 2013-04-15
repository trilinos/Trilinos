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

#include "Trios_nssi_fprint_types.h"
#include "Trios_nnti_fprint_types.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include <mpi.h>
#include <libtopomap.hpp>

#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <string>
using namespace std;

int placement_debug_level = LOG_UNDEFINED;

string strategy("ascending");

enum graph_connection_t {
    GRAPH_COMPLETE=0,
    GRAPH_CLIENT_COMPLETE,
    GRAPH_SERVER_COMPLETE,
    GRAPH_CLIENT_SERVER_ONLY
};


/* ----------------- COMMAND-LINE OPTIONS --------------- */

#define DEBUG_PLACEMENT

void construct_graph(
        int *rank_map,
        int *nid_map,
        int num_servers,
        int num_clients,
        int servers_per_node,
        int clients_per_node,
        int server_weight,
        int client_weight,
        int client_server_weight,
        enum graph_connection_t connection,
        int passthru)
{
    int npes, me, i, status;
    int num_neighs;
    std::vector<int> neighbors;
    std::vector<int> weights;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    const int len=128;
    char my_hostname[len];
    int my_nid;

    gethostname(my_hostname, len);
    my_nid = atoi(my_hostname+3);

//#if defined(DEBUG_PLACEMENT)
//    fprintf(stdout, "%i: hostname=%s, nid=%i\n", me, my_hostname, my_nid);
//    sleep(2);
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif

    // Original rank_map
    for (i = 0; i < npes; i++) {
        rank_map[i] = i;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (me < num_servers) {
        if ((connection == GRAPH_COMPLETE) ||
            (connection == GRAPH_SERVER_COMPLETE)) {
            // construct a complete graph for servers-to-server communication
            for (i=0;i<num_servers;i++) {
                if (i != me) {
                    neighbors.push_back(i);
                    weights.push_back(server_weight);
                }
            }
        }
    } else {
        if ((connection == GRAPH_COMPLETE) ||
            (connection == GRAPH_CLIENT_COMPLETE)) {
            // construct a complete graph for client-to-client communication
            for (i=num_servers;i<npes;i++) {
                if (i != me) {
                    neighbors.push_back(i);
                    weights.push_back(client_weight);
                }
            }
        }
        if (num_servers > 0) {
            // clients also communicate with one server
            int client_rank=(me-num_servers);       // my rank within the clients
            int bin_size=(num_clients/num_servers); // clients per server
            int my_server=client_rank/bin_size;
            neighbors.push_back(my_server);
            weights.push_back(client_server_weight);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    num_neighs=neighbors.size();

#if defined(DEBUG_PLACEMENT)
    for (i=0;i<npes;i++) {
        if (i == me) {
            std::vector<int>::iterator iter=neighbors.begin();
            for (;iter<neighbors.end();iter++) {
                fprintf(stdout, "%d: neighbor=%d\n", me, *iter);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Create a MPI distributed graph communicator
    MPI_Comm distgr;
    status = MPIX_Dist_graph_create(
                MPI_COMM_WORLD,
                1,
                &me,
                &num_neighs,
                &neighbors[0],
                &weights[0],
                MPI_INFO_NULL,
                0,
                &distgr);
    if (status != MPI_SUCCESS) {
        fprintf(stdout, "%i: MPI_Dist_graph_create() failed.\n", me);
        exit(-1);
    }

#if defined(DEBUG_PLACEMENT)
{
    // Make sure our rank hasn't changed!
    int dgpe;
    MPI_Comm_rank(distgr, &dgpe);
    if (dgpe != me) {
        fprintf(stdout, "%d: dgpe != me (%i != %i)\n", me, dgpe, me);
        exit(-1);
    }
}
#endif

    // Call libtopomap
    setenv("TPM_ROUTING_FILE", "routes.txt", 1);
    if (passthru) {
        setenv("TPM_STRATEGY", "none", 1);
        setenv("TPM_FIX_GRAPH", "no", 1);
        setenv("TPM_PARMETIS", "no", 1);
        setenv("TPM_FIX_PARMETIS", "no", 1);
        setenv("TPM_ANNEAL", "no", 1);
    } else {
        setenv("TPM_STRATEGY", strategy.c_str(), 1);
    }
    int new_rank;
    status = TPM_Topomap(distgr, "topomap.txt", 0, &new_rank);
    if (status != 0) {
        fprintf(stdout, "%i: TPM_Topomap() failed.\n", me);
        exit(-1);
    }

//#if defined(DEBUG_PLACEMENT)
//    fprintf(stdout, "%i: my new rank = %i\n", me, new_rank);
//    sleep(2);
//    MPI_Barrier(MPI_COMM_WORLD);
//#endif

    MPI_Comm improved_placement_comm;
    MPI_Comm_split(MPI_COMM_WORLD, 1, new_rank, &improved_placement_comm);
    int improved_npes=-1;
    int improved_me=-1;
    MPI_Comm_size(improved_placement_comm, &improved_npes);
    MPI_Comm_rank(improved_placement_comm, &improved_me);
    if (new_rank != improved_me) {
        fprintf(stdout, "%d: new_rank != improved_me (%i != %i)\n", me, new_rank, improved_me);
        exit(-1);
    }

#if defined(DEBUG_PLACEMENT)
    for (i=0;i<npes;i++) {
        if (i == me) {
            fprintf(stdout, "i = %3d ; me = %3d ; new_rank = %3d ; improved_me = %3d\n", i, me, new_rank, improved_me);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Gather the new rank_map
    status = MPI_Allgather(&new_rank, 1, MPI_INT, rank_map, 1, MPI_INT, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stdout, "%i: MPI_Allgather() failed.\n", me);
    }

    // Gather the nid map... nid_map[MPI_COMM_WORLD Rank] = nid
    status = MPI_Allgather(&my_nid, 1, MPI_INT, nid_map, 1, MPI_INT, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stdout, "%i: MPI_Allgather() failed. 2\n", me);
    }

    if (!me) {
        fprintf(stdout, "New Process to Rank/Block map:\n");
        for (i = 0; i < npes; i++) {
            fprintf(stdout, "\trank_map[%d]=%d(nid=%d)\n", i, rank_map[i], nid_map[i]);
        }
        fprintf(stdout, "\n");
    }
}


int
main (int argc, char *argv[])
{
    // command-line arguments
    log_level debug_level = LOG_ERROR;
    string logfile("");

    int npes, me, i;

    int num_servers=1;
    int num_clients=1;

    int servers_per_node=1;
    int clients_per_node=1;

    int client_weight=10;
    int server_weight=10;
    int client_server_weight=5;

    string server_node_file("SNF.txt");
    string client_node_file("CNF.txt");

    const int num_graphs = 4;
    const int graph_vals[] = {
            GRAPH_COMPLETE,
            GRAPH_CLIENT_COMPLETE,
            GRAPH_SERVER_COMPLETE,
            GRAPH_CLIENT_SERVER_ONLY
    };
    const char * graph_names[] = {
            "complete",
            "client-complete",
            "server-complete",
            "client-server-only"
    };
    enum graph_connection_t graph_connection=GRAPH_COMPLETE;

    MPI_Init(&argc, &argv);

    try {
        Teuchos::CommandLineProcessor parser;

        // init parser
        parser.setDocString("Find node placement of server and client ranks");

        parser.setOption("strategy", &strategy, "LibTopoMap strategy (greedy, greedy_route, recursive, rcm, scotch, ascending)");
        parser.setOption("num-servers", (int *)(&num_servers), "Number of servers to place");
        parser.setOption("num-clients", (int *)(&num_clients), "Number of clients to place");
        parser.setOption("servers-per-node", (int *)(&servers_per_node), "Number of server ranks per compute node");
        parser.setOption("clients-per-node", (int *)(&clients_per_node), "Number of client ranks per compute node");
        parser.setOption("server-weight", (int *)(&server_weight), "Edge weight of server-to-server communication");
        parser.setOption("client-weight", (int *)(&client_weight), "Edge weight of client-to-client communication");
        parser.setOption("client-server-weight", (int *)(&client_server_weight), "Edge weight of client-to-server communication");
        parser.setOption("server-node-file", &server_node_file, "Where to write the server placement results");
        parser.setOption("client-node-file", &client_node_file, "Where to write the client placement results");
        parser.setOption("verbose", (int *)(&debug_level), "Debug level");
        parser.setOption("logfile", &logfile, "Path to file for debug statements");
        // Set an enumeration command line option for the connection graph
        parser.setOption("graph-connection", (int*)&graph_connection, num_graphs, graph_vals, graph_names,
                "Graph Connections for the example: \n"
                "\t\t\tcomplete : client-client graph is complete, server-server graph is complete\n"
                "\t\t\tclient-complete: client-client graph is complete, server-server graph is empty\n"
                "\t\t\tserver-complete : client-client graph is empty, server-server graph is complete\n"
                "\t\t\tclient-server-only: client-client graph is empty, server-server graph is empty\n"
                "\t\t\tIn all cases, each client has an edge to one of the servers\n"
                );

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

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    if (me==0) {
        cout << " ----------------  ARGUMENTS --------------- " << std::endl;
        cout << " \tstrategy             = " << strategy << std::endl;
        cout << " \tgraph-connection     = " << graph_names[graph_connection] << std::endl;
        cout << " \tnum-servers          = " << num_servers << std::endl;
        cout << " \tnum-clients          = " << num_clients << std::endl;
        cout << " \tservers-per-node     = " << servers_per_node << std::endl;
        cout << " \tclients-per-node     = " << clients_per_node << std::endl;
        cout << " \tserver-weight        = " << server_weight << std::endl;
        cout << " \tclient-weight        = " << client_weight << std::endl;
        cout << " \tclient-server-weight = " << client_server_weight << std::endl;
        cout << " \tserver-node-file     = " << server_node_file << std::endl;
        cout << " \tclient-node-file     = " << client_node_file << std::endl;
        cout << " \tverbose              = " << debug_level << std::endl;
        cout << " \tlogfile              = " << logfile << std::endl;
        cout << " ------------------------------------------- " << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int *rank_map=(int*)malloc(sizeof(int) * npes);
    int *nid_map=(int*)malloc(sizeof(int) * npes);

    construct_graph(
            rank_map,
            nid_map,
            num_servers,
            num_clients,
            servers_per_node,
            clients_per_node,
            server_weight,
            client_weight,
            client_server_weight,
            graph_connection,
            0);

    if (me == 0) {
        ofstream snf(server_node_file.c_str(), ios_base::out);
        ofstream cnf(client_node_file.c_str(), ios_base::out);

        for (i=0;i<npes;i++) {
            if (rank_map[i] < num_servers)
                snf << nid_map[i] << "\t" << i << "\t" << rank_map[i] << std::endl;
        }
        for (i=0;i<npes;i++) {
            if (rank_map[i] >= num_servers)
                cnf << nid_map[i] << "\t" << i << "\t" << rank_map[i] << std::endl;
        }

        snf.close();
        cnf.close();
    }

    MPI_Finalize();

    return 0;
}
