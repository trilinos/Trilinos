// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <vector>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <type_traits>      // std::is_same


namespace Tpetra {

template<class LocalOrdinal, class GlobalOrdinal, class Node>
class crsGraph_Swap_Tester
{
    using Scalar                   = int;                                                   // Not really used for CrsGraph construction but
                                                                                            // we keep this around since the CrsMatrix::swap
                                                                                            // test is constructing its graphs the same way
                                                                                            // so hanging onto this to keep the structs in
                                                                                            // this file happy (this does not get passed to
                                                                                            // anything in Tpetra in this test).
    using comm_type                = Teuchos::RCP<const Teuchos::Comm<int>>;                // The comm type
    using graph_type               = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;   // Tpetra CrsGraph type
    using pair_owner_type          = std::pair<GlobalOrdinal, int>;                         // For row owners, pairs are (rowid, comm rank)
    using vec_owners_type          = std::vector<pair_owner_type>;                          // For vectors of owners
    using tuple_weighted_edge_type = std::tuple<GlobalOrdinal, GlobalOrdinal, Scalar>;      // Weighted edges (u,v,w)
    using vec_weighted_edge_type   = std::vector<tuple_weighted_edge_type>;                 // Vector of weighted edges

  public:
    void execute(Teuchos::FancyOStream& out, bool& success)
    {
        using std::endl;
        using Teuchos::Comm;
        using Teuchos::outArg;
        using Teuchos::RCP;

        bool verbose = Tpetra::Details::Behavior::verbose();

        comm_type initialComm = getDefaultComm();

        TEUCHOS_TEST_FOR_EXCEPTION(initialComm->getSize() < 2, std::runtime_error, "This test requires at least two processors.");

        // Set up a communicator that has exactly two processors in it for the actual test.
        const int color = (initialComm->getRank() < 2) ? 0 : 1;
        const int key   = 0;
        comm_type comm  = initialComm->split(color, key);

        // If I am involved in this test (i.e., my pid == 0 or 1)
        if(0 == color)
        {
            const int num_procs = comm->getSize();
            const int my_rank   = comm->getRank();
            assert(num_procs == 2);

            out << ">>> CrsGraph::swap() Unit Test" << std::endl;
            out << ">>> num_procs: " << num_procs << std::endl;

            success = true;

            //--------------------------------------------------------------------
            //
            // Construct CrsGraphs & CrsMatrices
            //
            //--------------------------------------------------------------------
            out << std::endl
                << ">>> ----------------------------------------------------------" << std::endl
                << ">>> -" << std::endl
                << ">>> - Construct CrsGraphs & CrsMatrices" << std::endl
                << ">>> -" << std::endl
                << ">>> ----------------------------------------------------------" << std::endl
                << std::endl;

            // Create weighted edge pairs for graph/matrix type 1
            vec_weighted_edge_type vec_edges_wgt = {
              std::make_tuple(0, 0, 0),      std::make_tuple(0, 11, 11),      // row 0
              std::make_tuple(1, 3, 103),    std::make_tuple(1, 4, 104),      // row 1
              std::make_tuple(3, 2, 302),                                     // row 3
              std::make_tuple(7, 5, 705),    std::make_tuple(7, 7, 707),      // row 7
              std::make_tuple(10, 6, 1006)                                    // row 10
            };

            vec_owners_type vec_owners = {
              pair_owner_type(0, 0),   pair_owner_type(1, 0),   pair_owner_type(3, 0),    // p0
              pair_owner_type(7, 1),   pair_owner_type(10, 1)                             // p1
            };

            if(verbose)
            {
                for(auto& p: vec_edges_wgt)
                    out << "ew: " << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
                for(auto& p: vec_owners) out << "o : " << p.first << ", " << p.second << std::endl;
            }

            if(my_rank == 0)
                out << ">>> create graph_a" << std::endl;
            RCP<graph_type> graph_a = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_a->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create graph_b" << std::endl;
            RCP<graph_type> graph_b = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_b->describe(out, Teuchos::VERB_DEFAULT);

            vec_edges_wgt.clear();
            vec_owners.clear();

            vec_edges_wgt = {
              std::make_tuple(0, 0, 0),      std::make_tuple(0, 11, 11),      // row 0
              std::make_tuple(1, 7, 107),    std::make_tuple(1, 8, 108),      // row 1
              std::make_tuple(3, 1, 301),                                     // row 3
              std::make_tuple(7, 5, 705),                                     // row 7
              std::make_tuple(10, 4, 1004)                                    // row 10
            };

            vec_owners = {
              pair_owner_type(0, 0),   pair_owner_type(1, 0),                               // p0
              pair_owner_type(3, 1),   pair_owner_type(7, 1),   pair_owner_type(10, 1)      // p1
            };

            if(my_rank == 0)
                out << ">>> create graph_c" << std::endl;
            RCP<graph_type> graph_c = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_c->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create graph_d" << std::endl;
            RCP<graph_type> graph_d = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_d->describe(out, Teuchos::VERB_DEFAULT);

            //--------------------------------------------------------------------
            //
            // Test CrsGraph::swap()
            //
            //--------------------------------------------------------------------
            out << std::endl
                << ">>> ----------------------------------------------------------" << std::endl
                << ">>> -" << std::endl
                << ">>> - Test CrsGraph::swap()" << std::endl
                << ">>> -" << std::endl
                << ">>> ----------------------------------------------------------" << std::endl
                << std::endl;

            // Verify the initial matching state of the graphs and matrices.
            TEST_EQUALITY(graph_a->isIdenticalTo(*graph_b), true);       // graph_a and graph_b should be the same
            TEST_EQUALITY(graph_c->isIdenticalTo(*graph_d), true);       // graph_c and graph_d should be the same
            TEST_EQUALITY(graph_a->isIdenticalTo(*graph_c), false);      // graph_a and graph_c should be different
            TEST_EQUALITY(graph_b->isIdenticalTo(*graph_d), false);      // graph_b and graph_d should be different

            //
            // Swap graph b and c
            //
            out << std::endl << ">>> swap B and C so that A!=B, B!=C, A==C, B==D" << std::endl;
            out << ">>> swap graph_b and graph_c" << std::endl;
            graph_b->swap(*graph_c);

            // Verify that the graph and matrices were swapped correctly.
            TEST_EQUALITY(graph_a->isIdenticalTo(*graph_b), false);      // graph_a and graph_b should be different
            TEST_EQUALITY(graph_c->isIdenticalTo(*graph_d), false);      // graph_c and graph_d should be different
            TEST_EQUALITY(graph_a->isIdenticalTo(*graph_c), true);       // graph_a and graph_c should be the same
            TEST_EQUALITY(graph_b->isIdenticalTo(*graph_d), true);       // graph_b and graph_d should be the same

            //
            // Swap graph b and c again.
            //
            out << std::endl << ">>> swap B and C back to original state." << std::endl;
            out << ">>> swap graph_b and graph_c" << std::endl;
            graph_b->swap(*graph_c);

            // Verify the initial identical-to state of the graphs (they should be back to the original state)
            TEST_EQUALITY(graph_a->isIdenticalTo(*graph_b), true);       // graph_a and graph_b should be the same
            TEST_EQUALITY(graph_c->isIdenticalTo(*graph_d), true);       // graph_c and graph_d should be the same
            TEST_EQUALITY(graph_a->isIdenticalTo(*graph_c), false);      // graph_a and graph_c should be different
            TEST_EQUALITY(graph_b->isIdenticalTo(*graph_d), false);      // graph_b and graph_d should be different
        }      // if(0==color) ...   i.e., if I am in the active communicator group for the test.
    }

  private:
    // comm           : a communicator with exactly 2 processes in it.
    // gbl_wgt_edges  : [ (u,v,w), ... ]  - weighted edges: u->v with weight w.
    // row_owners     : [ (u, pid), ... ]
    // gbl_num_columns: Max # of columns in the matrix-representation of the graph.
    //                  This should be >= the highest value of v from all edges (u,v) in edges.
    //                  Note: u and v are 0-indexed, so if the highest v is 11, then this should be 12.
    Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
          generate_crsgraph(Teuchos::RCP<const Teuchos::Comm<int>>&                              comm,
                            const std::vector<std::tuple<GlobalOrdinal, GlobalOrdinal, Scalar>>& gbl_wgt_edges,
                            const std::vector<std::pair<GlobalOrdinal, int>>&                    gbl_row_owners,
                            const size_t                                                         gbl_num_columns,
                            const bool                                                           do_fillComplete = true)
    {
        using Teuchos::Comm;
        using Teuchos::RCP;

        using map_type             = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;           // Tpetra Map type
        using map_rows_type        = std::map<GlobalOrdinal, int>;                   // map rows to pid's
        using vec_go_type          = std::vector<GlobalOrdinal>;                     // vector of GlobalOrdinals
        using map_row_to_cols_type = std::map<GlobalOrdinal, vec_go_type>;           // Map rows to columns

        const bool verbose   = Tpetra::Details::Behavior::verbose();
        const int  comm_rank = comm->getRank();

        // fill in the row-to-columns map
        map_row_to_cols_type gbl_rows;
        for(auto& e: gbl_wgt_edges)
        {
            if(gbl_rows.find(std::get<0>(e)) == gbl_rows.end())
            {
                gbl_rows[ std::get<0>(e) ] = vec_go_type();
            }
            gbl_rows[ std::get<0>(e) ].push_back(std::get<1>(e));
        }

        // fill in the row owner map
        map_rows_type gbl_row2pid;
        for(auto& p: gbl_row_owners) { gbl_row2pid.insert(p); }

        // Print out some debugging information on what's in the graph
        if(verbose)
        {
            std::cout << "p=0 | gbl_num_rows: " << gbl_rows.size() << std::endl;
            for(auto& p: gbl_rows)
            {
                std::cout << "p=0 | gbl_row     : " << std::get<0>(p) << " (" << std::get<1>(p).size() << ") ";
                for(auto& j: std::get<1>(p)) { std::cout << j << " "; }
                std::cout << std::endl;
            }
            for(auto& p: gbl_row2pid) std::cout << "p=0 | gbl_row2pid : " << p.first << " => " << p.second << std::endl;
        }

        GlobalOrdinal gbl_num_rows = gbl_rows.size();      // the number of global rows
        LocalOrdinal lcl_num_rows = 0;                    // this will be updated later

        // How many local rows do I have?
        for(auto& r: gbl_row2pid)
        {
            if(comm_rank == r.second)
            {
                lcl_num_rows++;
            }
        }

        if(verbose)
            std::cout << "p=" << comm_rank << " | "
                      << "lcl_num_rows = " << lcl_num_rows << std::endl;

        // Set up global ids
        std::vector<GlobalOrdinal> global_ids;
        for(auto& r: gbl_row2pid)
        {
            if(comm_rank == r.second)
            {
                global_ids.push_back(r.first);
            }
        }

        if(verbose)
        {
            for(size_t i = 0; i < global_ids.size(); i++)
            {
                std::cout << "p=" << comm_rank << " | "
                          << "global_ids[" << i << "] = " << global_ids[ i ] << std::endl;
            }
            std::cout << "p=" << comm_rank << " | "
                      << "row_map = map_type(" << gbl_num_rows << ", "
                      << "global_ids.data(), " << lcl_num_rows << ", 0, comm)" << std::endl;
        }

        // Create the Row Map
        RCP<const map_type> row_map(new map_type(gbl_num_rows, global_ids.data(), lcl_num_rows, 0, comm));

        Teuchos::ArrayRCP<size_t> num_ent_per_row(lcl_num_rows);
        size_t                    idx = 0;
        for(auto& r: gbl_rows)
        {
            const GlobalOrdinal  irow    = r.first;
            const int row_pid = gbl_row2pid.find(irow)->second;
            if(comm_rank == row_pid)
            {
                num_ent_per_row[ idx++ ] = r.second.size();
            }
        }

        if(verbose)
        {
            std::cout << "p=" << comm_rank << " | lcl_num_rows = " << lcl_num_rows << std::endl;
            for(int i = 0; i < lcl_num_rows; i++)
            {
                std::cout << "p=" << comm_rank << " | "
                          << "num_ent_per_row[" << i << "] = " << num_ent_per_row[ i ] << std::endl;
            }
        }

        RCP<graph_type> output_graph(new graph_type(row_map, num_ent_per_row ()));

        for(auto& r: gbl_rows)
        {
            const GlobalOrdinal  irow = r.first;
            const int pid  = gbl_row2pid.find(irow)->second;
            if(comm_rank == pid)
            {
                std::vector<GlobalOrdinal> gbl_inds;
                for(auto& v: r.second) { gbl_inds.push_back(v); }
                if(verbose)
                {
                    std::cout << "p=" << comm_rank << " | "
                              << "gbl_inds size = " << r.second.size() << std::endl;
                    for(size_t i = 0; i < gbl_inds.size(); i++)
                    {
                        std::cout << "p=" << comm_rank << " | "
                                  << "gbl_inds[" << i << "] = " << gbl_inds[ i ] << std::endl;
                    }
                    std::cout << "p=" << comm_rank << " | "
                              << "insertGlobalIndices(" << irow << ", " << r.second.size() << ", gbl_inds.data())" << std::endl;
                }
                output_graph->insertGlobalIndices(irow, r.second.size(), gbl_inds.data());
            }
        }

        RCP<const map_type> range_map = row_map;

        const GlobalOrdinal index_base = 0;
        RCP<const map_type> domain_map(new map_type(gbl_num_columns, index_base, comm));

        if(do_fillComplete)
        {
            output_graph->fillComplete(domain_map, range_map);
        }

        return output_graph;
    }

};      // class crsGraph_Swap_Tester

}      // namespace Tpetra





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
using Tpetra::TestingUtilities::getDefaultComm;
typedef Tpetra::global_size_t GST;

double      errorTolSlack = 1e+1;
std::string filedir;

#define STD_TESTS(graph)                                                                                        \
    {                                                                                                           \
        auto   STCOMM   = graph.getComm();                                                                      \
        auto   STMYGIDS = graph.getRowMap()->getLocalElementList();                                              \
        size_t STMAX    = 0;                                                                                    \
                                                                                                                \
        for(size_t STR = 0; STR < graph.getLocalNumRows(); ++STR)                                                \
        {                                                                                                       \
            TEST_EQUALITY(graph.getNumEntriesInLocalRow(STR), graph.getNumEntriesInGlobalRow(STMYGIDS[ STR ])); \
            STMAX = std::max(STMAX, graph.getNumEntriesInLocalRow(STR));                                        \
        }                                                                                                       \
        TEST_EQUALITY(graph.getLocalMaxNumRowEntries(), STMAX);                                                  \
        GST STGMAX;                                                                                             \
        Teuchos::reduceAll<int, GST>(*STCOMM, Teuchos::REDUCE_MAX, STMAX, Teuchos::outArg(STGMAX));             \
        TEST_EQUALITY(graph.getGlobalMaxNumRowEntries(), STGMAX);                                               \
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


//
// UNIT TESTS
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, Swap, LO, GO, Node)
{
    auto crsMatrixTester = Tpetra::crsGraph_Swap_Tester<LO, GO, Node>();
    crsMatrixTester.execute(out, success);
}


//
// INSTANTIATIONS
//



// Tests to build and run in both debug and release modes.  We will
// instantiate them over all enabled local ordinal (LO), global
// ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP_DEBUG_AND_RELEASE(LO, GO, NODE) TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, Swap, LO, GO, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP_DEBUG_AND_RELEASE)

}      // namespace
