/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/
#include <map>
#include <set>
#include <unordered_map>
#include <utility>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <type_traits>      // std::is_same

namespace {      // (anonymous)

using std::endl;
using std::max;
using std::min;
using Teuchos::arcp;
using Teuchos::arcpClone;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::null;
using Teuchos::outArg;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;
using Teuchos::VERB_EXTREME;
using Teuchos::VERB_HIGH;
using Teuchos::VERB_LOW;
using Teuchos::VERB_MEDIUM;
using Teuchos::VERB_NONE;
using Tpetra::DynamicProfile;
using Tpetra::ProfileType;
using Tpetra::StaticProfile;
using Tpetra::TestingUtilities::getDefaultComm;
typedef Tpetra::global_size_t GST;

double      errorTolSlack = 1e+1;
std::string filedir;

#define STD_TESTS(graph)                                                                                      \
    {                                                                                                         \
        auto   STCOMM   = graph.getComm();                                                                    \
        auto   STMYGIDS = graph.getRowMap()->getNodeElementList();                                            \
        size_t STMAX    = 0;                                                                                  \
                                                                                                              \
        for(size_t STR = 0; STR < graph.getNodeNumRows(); ++STR)                                              \
        {                                                                                                     \
            TEST_EQUALITY(graph.getNumEntriesInLocalRow(STR), graph.getNumEntriesInGlobalRow(STMYGIDS[STR])); \
            STMAX = std::max(STMAX, graph.getNumEntriesInLocalRow(STR));                                      \
        }                                                                                                     \
        TEST_EQUALITY(graph.getNodeMaxNumRowEntries(), STMAX);                                                \
        GST STGMAX;                                                                                           \
        Teuchos::reduceAll<int, GST>(*STCOMM, Teuchos::REDUCE_MAX, STMAX, Teuchos::outArg(STGMAX));           \
        TEST_EQUALITY(graph.getGlobalMaxNumRowEntries(), STGMAX);                                             \
    }


TEUCHOS_STATIC_SETUP()
{
    Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption("filedir", &filedir, "Directory of expected input files.");
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi",
                  "test-serial",
                  &Tpetra::TestingUtilities::testMpi,
                  "Test MPI (if available) or force test of serial.  In a serial build,"
                  " this option is ignored and a serial comm is always used.");
    clp.setOption("error-tol-slack", &errorTolSlack, "Slack off of machine epsilon used to check test results");
}
// todo: update options so that this test won't run in serial (it's not set up for that currently)




// comm           : a communicator with exactly 2 processes in it.
// edges          : [ (u,v), ... ]
// row_owners     : [ (u, pid), ... ]
// gbl_num_columns: Max # of columns in the matrix-representation of the graph.
//                  This should be >= the highest value of v from all edges (u,v) in edges.
//                  Note: u and v are 0-indexed, so if the highest v is 11, then this should be 12.
template<class LO, class GO, class Node, class comm_t>
Teuchos::RCP<Tpetra::CrsGraph<LO, GO, Node>>
generate_crsgraph(Teuchos::RCP<comm_t>                      comm,
                  const std::vector<std::pair<GO, GO>>&     gbl_edges,
                  const std::vector<std::pair<GO, size_t>>& gbl_row_owners,
                  const size_t                              gbl_num_columns,
                  const bool                                do_fillComplete=true)
{
    using Teuchos::Comm;

    typedef Tpetra::CrsGraph<LO, GO, Node> graph_t;
    typedef Tpetra::Map<LO, GO, Node>      map_t;

    typedef typename std::map<GO, size_t> map_rows_t;      // row_id, num_entries


    typedef typename std::vector<GO>        vec_go_t;
    typedef typename std::map<GO, vec_go_t> map_row_to_cols_t;

    const size_t comm_rank = comm->getRank();

    map_row_to_cols_t gbl_rows;
    for(auto& e: gbl_edges)
    {
        if(gbl_rows.find(e.first) == gbl_rows.end())
        {
            gbl_rows[ e.first ] = vec_go_t();
        }
        gbl_rows[ e.first ].push_back(e.second);
    }

    map_rows_t gbl_row2pid;
    for(auto& p: gbl_row_owners) { gbl_row2pid.insert(p); }

    // Print out some debugging information on what's in the graph
    #if 0
    std::cout << "p=0 | gbl_num_rows: " << gbl_rows.size() << std::endl;
    for(auto& p: gbl_rows)
    {
        std::cout << "p=0 | gbl_row     : " << p.first << " (" << p.second.size() << ") ";
        for(auto& j: p.second)
        {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }
    for(auto& p: gbl_row2pid)
        std::cout << "p=0 | gbl_row2pid : " << p.first << " => " << p.second << std::endl;
    #endif

    GO gbl_num_rows = gbl_rows.size();      // the number of global rows
    LO lcl_num_rows = 0;                    // this will be updated later

    // How many local rows do I have?
    for(auto& r: gbl_row2pid)
    {
        if(comm_rank == r.second)
        {
            lcl_num_rows++;
        }
    }

    #if 0
    std::cout << "p=" << comm_rank << " | " << "lcl_num_rows = " << lcl_num_rows << std::endl;
    #endif

    // Set up global ids
    std::vector<GO> global_ids;
    for(auto& r: gbl_row2pid)
    {
        if(comm_rank == r.second)
        {
            global_ids.push_back(r.first);
        }
    }
    #if 0
    for(size_t i = 0; i < global_ids.size(); i++)
    {
        std::cout << "p=" << comm_rank << " | " << "global_ids[" << i << "] = " << global_ids[i] << std::endl;
    }
    #endif

    // Create the Row Map
    #if 0
    std::cout << "p=" << comm_rank << " | " << "row_map = map_t(" << gbl_num_rows << ", " << "global_ids.data(), "
              << lcl_num_rows << ", 0, comm)" << std::endl;
    #endif
    RCP<const map_t> row_map(new map_t(gbl_num_rows, global_ids.data(), lcl_num_rows, 0, comm));

    Teuchos::ArrayRCP<size_t> num_ent_per_row(lcl_num_rows);
    size_t                    idx = 0;
    for(auto& r: gbl_rows)
    {
        const GO     irow    = r.first;
        const size_t row_pid = gbl_row2pid.find(irow)->second;
        if(row_pid == comm_rank)
        {
            num_ent_per_row[ idx++ ] = r.second.size();
        }
    }

    #if 0
    std::cout << "p=" << comm_rank << " | lcl_num_rows = " << lcl_num_rows << std::endl;
    for(int i=0; i<lcl_num_rows; i++)
    {
        std::cout << "p=" << comm_rank << " | " << "num_ent_per_row[" << i << "] = " << num_ent_per_row[i] << std::endl;
    }
    #endif

    RCP<graph_t> output_graph(new graph_t(row_map, num_ent_per_row, Tpetra::StaticProfile));

    for(auto& r: gbl_rows)
    {
        const GO     irow = r.first;
        const size_t pid  = gbl_row2pid.find(irow)->second;
        if(pid == comm_rank)
        {
            std::vector<GO> gbl_inds;
            for(auto& v: r.second)
            {
                // std::cout << "p=" << comm_rank << " | " << "gbl_inds.push_back(" << v << ")" << std::endl;
                gbl_inds.push_back(v);
            }

            #if 0
            std::cout << "p=" << comm_rank << " | " << "gbl_inds size = " << r.second.size() << std::endl;
            for(size_t i=0; i<gbl_inds.size(); i++)
            {
                std::cout << "p=" << comm_rank << " | " << "gbl_inds[" << i << "] = " << gbl_inds[i] << std::endl;
            }
            #endif

            // std::cout << "p=" << comm_rank << " | " << "insertGlobalIndices(" << irow << ", " << r.second.size() << ",
            // gbl_inds.data())" << std::endl;
            output_graph->insertGlobalIndices(irow, r.second.size(), gbl_inds.data());
        }
    }

    RCP<const map_t> range_map = row_map;

    const GO         index_base = 0;
    RCP<const map_t> domain_map(new map_t(gbl_num_columns, index_base, comm));

    if(do_fillComplete)
    {
        output_graph->fillComplete(domain_map, range_map);
    }

    return output_graph;
}



//
// UNIT TESTS
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, Swap, LO, GO, Node)
{
    using std::endl;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;

    typedef Teuchos::Comm<int>             comm_t;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_t;

    typedef typename std::pair<GO, GO>     pair_edge_t;
    typedef typename std::pair<GO, size_t> pair_owner_t;

    typedef typename std::vector<pair_edge_t>  vec_edges_t;
    typedef typename std::vector<pair_owner_t> vec_owners_t;

    bool verbose = Tpetra::Details::Behavior::verbose();

    auto initialComm = getDefaultComm();
    TEUCHOS_TEST_FOR_EXCEPTION(initialComm->getSize() < 2, std::runtime_error, "This test requires at least two processors.");

    // Set up a communicator that has exactly two processors in it for the actual test.
    const int color = (initialComm->getRank() < 2) ? 0 : 1;
    const int key   = 0;
    auto      comm  = initialComm->split(color, key);

    // If I am involved in this test (i.e., my pid == 0 or 1)
    if(0 == color)
    {
        const int num_procs = comm->getSize();
        assert(num_procs == 2);

        out << ">>> CrsGraph::swap() Unit Test" << std::endl;
        out << ">>> num_procs: " << num_procs << std::endl;

        success = true;

        vec_edges_t  vec_edges  = {pair_edge_t(0, 0), pair_edge_t(0, 11), pair_edge_t(1, 3), pair_edge_t(1, 4),
                                   pair_edge_t(3, 2), pair_edge_t(7, 5),  pair_edge_t(7, 7), pair_edge_t(10, 6)};
        vec_owners_t vec_owners = {pair_owner_t(0, 0), pair_owner_t(1, 0), pair_owner_t(3, 0), pair_owner_t(7, 1),
                                   pair_owner_t(10, 1)};

         if(verbose)
         {
            for(auto& p: vec_edges)  out << "e: " << p.first << ", " << p.second << std::endl;
            for(auto& p: vec_owners) out << "o: " << p.first << ", " << p.second << std::endl;
         }

        out << ">>> create graph_a" << std::endl;
        RCP<graph_t> graph_a = generate_crsgraph<LO, GO, Node, comm_t>(comm, vec_edges, vec_owners, 12);
        //graph_a->describe(out, Teuchos::VERB_DEFAULT);

        out << ">>> create graph_b" << std::endl;
        RCP<graph_t> graph_b = generate_crsgraph<LO, GO, Node, comm_t>(comm, vec_edges, vec_owners, 12);
        //graph_b->describe(out, Teuchos::VERB_DEFAULT);

        vec_edges.clear();
        vec_owners.clear();

        vec_edges  = {pair_edge_t(0, 0), pair_edge_t(0, 11), pair_edge_t(1, 7), pair_edge_t(1, 8),
                      pair_edge_t(3, 1), pair_edge_t(7, 5), pair_edge_t(10, 4) };
        vec_owners = {pair_owner_t(0, 0), pair_owner_t(1, 0), pair_owner_t(3, 1), pair_owner_t(7, 1),
                      pair_owner_t(10, 1)};

        out << ">>> create graph_c" << std::endl;
        RCP<graph_t> graph_c = generate_crsgraph<LO, GO, Node, comm_t>(comm, vec_edges, vec_owners, 12);
        //graph_c->describe(out, Teuchos::VERB_DEFAULT);

        out << ">>> create graph_d" << std::endl;
        RCP<graph_t> graph_d = generate_crsgraph<LO, GO, Node, comm_t>(comm, vec_edges, vec_owners, 12);
        //graph_d->describe(out, Teuchos::VERB_DEFAULT);

        TEST_EQUALITY(graph_a->isIdenticalTo(*graph_b), true);       // graph_a and graph_b should be the same
        TEST_EQUALITY(graph_c->isIdenticalTo(*graph_d), true);       // graph_c and graph_d should be the same
        TEST_EQUALITY(graph_a->isIdenticalTo(*graph_c), false);      // graph_a and graph_c should be different
        TEST_EQUALITY(graph_b->isIdenticalTo(*graph_d), false);      // graph_b and graph_d should be different


        // Swap graph b and c
        out << ">>> swap graph_b and graph_c" << std::endl;
        graph_c->swap(*graph_b);

        TEST_EQUALITY(graph_a->isIdenticalTo(*graph_b), false);    // graph_a and graph_b should be different
        TEST_EQUALITY(graph_c->isIdenticalTo(*graph_d), false);    // graph_c and graph_d should be different
        TEST_EQUALITY(graph_a->isIdenticalTo(*graph_c), true);     // graph_a and graph_c should be the same
        TEST_EQUALITY(graph_b->isIdenticalTo(*graph_d), true);     // graph_b and graph_d should be the same
    }
}


  //
  // INSTANTIATIONS
  //



  // Tests to build and run in both debug and release modes.  We will
  // instantiate them over all enabled local ordinal (LO), global
  // ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP_DEBUG_AND_RELEASE(LO, GO, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, Swap, LO, GO, NODE)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP_DEBUG_AND_RELEASE)

  }      // namespace
