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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class crsMatrix_Swap_Tester
{
    using comm_type                = Teuchos::RCP<const Teuchos::Comm<int>>;                          // The comm type
    using graph_type               = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;             // Tpetra CrsGraph type
    using pair_owner_type          = std::pair<GlobalOrdinal, int>;                                   // For row owners, pairs are (rowid, comm rank)
    using vec_owners_type          = std::vector<pair_owner_type>;                                    // For vectors of owners
    using tuple_weighted_edge_type = std::tuple<GlobalOrdinal, GlobalOrdinal, Scalar>;                // Weighted edges (u,v,w)
    using vec_weighted_edge_type   = std::vector<tuple_weighted_edge_type>;                           // Vector of weighted edges
    using matrix_type              = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;    // Tpetra CrsMatrix type

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
              std::make_tuple(0, 0, 0),     std::make_tuple(0, 11, 11),      // row 0
              std::make_tuple(1, 3, 103),   std::make_tuple(1, 4, 104),      // row 1
              std::make_tuple(3, 2, 302),                                    // row 3
              std::make_tuple(7, 5, 705),   std::make_tuple(7, 7, 707),      // row 7
              std::make_tuple(10, 6, 1006)                                   // row 10
            };

            vec_owners_type vec_owners = {
              pair_owner_type(0, 0),   pair_owner_type(1, 0),   pair_owner_type(3, 0),      // p0
              pair_owner_type(7, 1),   pair_owner_type(10, 1)                               // p1
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
                out << ">>> create matrix_a" << std::endl;
            RCP<matrix_type> matrix_a = generate_crsmatrix(comm, graph_a, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                matrix_a->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create graph_b" << std::endl;
            RCP<graph_type> graph_b = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_b->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create matrix_b" << std::endl;
            RCP<matrix_type> matrix_b = generate_crsmatrix(comm, graph_b, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                matrix_b->describe(out, Teuchos::VERB_DEFAULT);

            vec_edges_wgt.clear();
            vec_owners.clear();

            vec_edges_wgt = {
              std::make_tuple(0, 0, 0),    std::make_tuple(0, 11, 11),      // row 0
              std::make_tuple(1, 7, 107),  std::make_tuple(1, 8, 108),      // row 1
              std::make_tuple(3, 1, 301),                                   // row 3
              std::make_tuple(7, 5, 705),                                   // row 7
              std::make_tuple(10, 4, 1004)                                  // row 10
            };

            vec_owners = {
              pair_owner_type(0, 0),   pair_owner_type(1, 0),                             // p0
              pair_owner_type(3, 1),   pair_owner_type(7, 1),   pair_owner_type(10, 1)    // p1
            };

            if(my_rank == 0)
                out << ">>> create graph_c" << std::endl;
            RCP<graph_type> graph_c = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_c->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create matrix_c" << std::endl;
            RCP<matrix_type> matrix_c = generate_crsmatrix(comm, graph_c, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                matrix_c->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create graph_d" << std::endl;
            RCP<graph_type> graph_d = generate_crsgraph(comm, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                graph_d->describe(out, Teuchos::VERB_DEFAULT);

            if(my_rank == 0)
                out << ">>> create matrix_d" << std::endl;
            RCP<matrix_type> matrix_d = generate_crsmatrix(comm, graph_d, vec_edges_wgt, vec_owners, 12);
            if(verbose)
                matrix_d->describe(out, Teuchos::VERB_DEFAULT);

            //--------------------------------------------------------------------
            //
            // Test CrsMatrix::swap()
            //
            //--------------------------------------------------------------------
            out << std::endl
                << ">>> ----------------------------------------------------------" << std::endl
                << ">>> -" << std::endl
                << ">>> - Test CrsMatrix::swap()" << std::endl
                << ">>> -" << std::endl
                << ">>> ----------------------------------------------------------" << std::endl
                << std::endl;

            //
            // Compare the matrices to verify the right initial state:
            // - A == B
            // - C == D
            // - A != C
            // - B != D
            //
            out << std::endl << ">>> Compare matrix_a == matrix_b" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_a, *matrix_b),
                          true);      // matrix_a and matrix_b should be the same

            out << std::endl << ">>> Compare matrix_c == matrix_d" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_c, *matrix_d),
                          true);      // matrix_c and matrix_d should be the same

            out << std::endl << ">>> Compare matrix_a != matrix_c" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_a, *matrix_c),
                          false);      // matrix_a and matrix_c should be different

            out << std::endl << ">>> Compare matrix_b != matrix_d" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_b, *matrix_d),
                          false);      // matrix_b and matrix_d should be different

            //
            // Swap matrix b and c
            //
            out << ">>> swap matrix_b and matrix_c" << std::endl << std::endl;
            matrix_b->swap(*matrix_c);

            //
            // Compare the matrices to verify the right initial state:
            // - A != B
            // - C != D
            // - A == C
            // - B == D
            //
            out << std::endl << ">>> Compare matrix_a != matrix_b" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_a, *matrix_b),
                          false);      // matrix_a and matrix_b should be different

            out << std::endl << ">>> Compare matrix_c != matrix_d" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_c, *matrix_d),
                          false);      // matrix_c and matrix_d should be different

            out << std::endl << ">>> Compare matrix_a == matrix_c" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_a, *matrix_c),
                          true);      // matrix_a and matrix_c should be the same

            out << std::endl << ">>> Compare matrix_b == matrix_d" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_b, *matrix_d),
                          true);      // matrix_b and matrix_d shoudl be the same

            //
            // Swap matrix b and c again
            //
            out << ">>> swap matrix_b and matrix_c" << std::endl << std::endl;
            matrix_c->swap(*matrix_b);

            //
            // Compare the matrices to verify the right initial state:
            // - A == B
            // - C == D
            // - A != C
            // - B != D
            //
            out << std::endl << ">>> Compare matrix_a == matrix_b" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_a, *matrix_b),
                          true);      // matrix_a and matrix_b should be the same

            out << std::endl << ">>> Compare matrix_c == matrix_d" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_c, *matrix_d),
                          true);      // matrix_c and matrix_d should be the same

            out << std::endl << ">>> Compare matrix_a != matrix_c" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_a, *matrix_c),
                          false);      // matrix_a and matrix_c should be different

            out << std::endl << ">>> Compare matrix_b != matrix_d" << std::endl;
            TEST_EQUALITY(compare_crsmatrix(comm, out, *matrix_b, *matrix_d),
                          false);      // matrix_b and matrix_d should be different

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
                      const bool do_fillComplete = true)
    {
        using Teuchos::Comm;
        using Teuchos::RCP;

        using map_type             = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;          // Tpetra Map type
        using map_rows_type        = std::map<GlobalOrdinal, int>;                            // map rows to pid's
        using vec_go_type          = std::vector<GlobalOrdinal>;                              // vector of GlobalOrdinals
        using map_row_to_cols_type = std::map<GlobalOrdinal, vec_go_type>;                    // Map rows to columns

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

        Teuchos::Array<size_t> num_ent_per_row(lcl_num_rows);
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



    // comm           : a communicator with exactly 2 processes in it.
    // gbl_wgt_edges  : [ (u,v,w), ... ]  - weighted edges: u->v with weight w.
    // row_owners     : [ (u, pid), ... ]
    // gbl_num_columns: Max # of columns in the matrix-representation of the graph.
    //                  This should be >= the highest value of v from all edges (u,v) in edges.
    //                  Note: u and v are 0-indexed, so if the highest v is 11, then this should be 12.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    generate_crsmatrix(Teuchos::RCP<const Teuchos::Comm<int>>&                                  comm,
                       const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>& graph,
                       const std::vector<std::tuple<GlobalOrdinal, GlobalOrdinal, Scalar>>&     gbl_wgt_edges,
                       const std::vector<std::pair<GlobalOrdinal, int>>&                        gbl_row_owners,
                       const size_t                                                             gbl_num_columns,
                       const bool                                                               do_fillComplete = true)
    {
        using Teuchos::Comm;
        using Teuchos::RCP;

        using map_type             = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;                    // Tpetra Map type
        using map_rows_type        = std::map<GlobalOrdinal, int>;                                      // map rows to pid's
        using vec_go_type          = std::vector<GlobalOrdinal>;                                        // vector of GlobalOrdinals
        using map_row_to_cols_type = std::map<GlobalOrdinal, vec_go_type>;                              // Map rows to columns

        using vec_scalar_type      = std::vector<Scalar>;                                               // Vector of Scalars
        using map_row_to_vals_type = std::map<GlobalOrdinal, vec_scalar_type>;

        // Get test verbosity from the settings
        const bool verbose = Tpetra::Details::Behavior::verbose();

        // Get my process id
        const auto comm_rank = comm->getRank();

        // fill in the row-to-columns map
        map_row_to_cols_type gbl_rows;      // map row to column ids (edges)
        map_row_to_vals_type gbl_vals;      // map row to values (edge weight)

        for(auto& e: gbl_wgt_edges)
        {
            if(gbl_rows.find(std::get<0>(e)) == gbl_rows.end())      // initialize row entries if it doesn't exist.
            {
                gbl_rows[ std::get<0>(e) ] = vec_go_type();
                gbl_vals[ std::get<0>(e) ] = vec_scalar_type();
            }
            gbl_rows[ std::get<0>(e) ].push_back(std::get<1>(e));      // add the column id (i.e., the edge)
            gbl_vals[ std::get<0>(e) ].push_back(std::get<2>(e));      // add the column value (i.e. the edge weight)
        }

        // fill in the row owner map
        map_rows_type gbl_row2pid;
        for(auto& p: gbl_row_owners) { gbl_row2pid.insert(p); }

        // Print out some debugging information on what's in the graph
        if(verbose && 0 == comm_rank)
        {
            std::cout << "p=" << comm_rank << " | gbl_num_rows: " << gbl_rows.size() << std::endl;
            for(auto& p: gbl_rows)
            {
                const auto  row_id                 = std::get<0>(p);
                const auto  num_row_entries        = std::get<1>(p).size();
                const auto* ptr_vec_col_entries    = &gbl_rows.at(row_id);
                const auto* ptr_vec_scalar_entries = &gbl_vals.at(row_id);

                std::cout << "p=" << comm_rank << " | " << row_id << " (" << num_row_entries << ") :\t ";
                for(size_t i = 0; i < num_row_entries; i++)
                {
                    const auto col_id = ptr_vec_col_entries->at(i);
                    const auto value  = ptr_vec_scalar_entries->at(i);
                    std::cout << col_id << "(" << value << ")\t ";
                }
                std::cout << std::endl;
            }

            for(auto& p: gbl_row2pid)
            { std::cout << "p=" << comm_rank << " | gbl_row2pid : " << p.first << " => " << p.second << std::endl; }
        }

        GlobalOrdinal gbl_num_rows = gbl_vals.size();      // the number of global rows
        LocalOrdinal lcl_num_rows = 0;                     // this will be updated later

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

        Teuchos::Array<size_t> num_ent_per_row(lcl_num_rows);
        size_t                    idx = 0;
        for(auto& r: gbl_vals)
        {
            const GlobalOrdinal  irow = r.first;
            const int row_pid         = gbl_row2pid.find(irow)->second;
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

        // Create the matrix using the input graph.
        RCP<matrix_type> output_matrix(new matrix_type(graph));

        for(auto& p: gbl_rows)
        {
            const auto row_id = std::get<0>(p);
            const auto pid    = gbl_row2pid.find(row_id)->second;
            if(comm_rank == pid)
            {
                const auto  num_row_entries        = std::get<1>(p).size();
                const auto* ptr_vec_col_entries    = &gbl_rows.at(row_id);
                const auto* ptr_vec_scalar_entries = &gbl_vals.at(row_id);

                Teuchos::Array<GlobalOrdinal> cols(gbl_rows.size());
                Teuchos::Array<Scalar>        vals(gbl_vals.size());
                for(size_t i = 0; i < num_row_entries; i++)
                {
                    cols[ i ] = ptr_vec_col_entries->at(i);
                    vals[ i ] = ptr_vec_scalar_entries->at(i);
                }

                output_matrix->sumIntoGlobalValues(row_id, cols, vals);
            }
        }

        // FillComplete the matrix.
        output_matrix->fillComplete();

        return output_matrix;
    }



    void print_crsmatrix(const Teuchos::RCP<const Teuchos::Comm<int>>&                        comm,
                         Teuchos::FancyOStream&                                               out,
                         const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&  m,
                         const std::string                                                    label)
    {
        using Teuchos::Comm;

        const int myRank = comm->getRank();

        // For debugging
        if(true)
        {
            if(0 == myRank)
            {
                out << std::endl;
                out << "=====[" << label << "]===============================================" << std::endl;
                out << "-----[RangeMap]--------------------" << std::endl;
            }
            out << m.getRangeMap().getRawPtr() << std::endl;
            m.getRangeMap()->describe(out, Teuchos::VERB_EXTREME);
            if(0 == myRank)
            {
                out << "=============================================================" << std::endl;
                out << std::endl;
            }
        }
    }



    bool compare_crsmatrix(const Teuchos::RCP<const Teuchos::Comm<int>>&  comm,
                           Teuchos::FancyOStream&                         out,
                           const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& matrix1,
                           const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& matrix2)
    {
        using Teuchos::Comm;

        bool output = true;

        if(!matrix1.isFillComplete() || !matrix2.isFillComplete())
        {
            if(output)
                out << std::endl;
            out << "Compare: FillComplete check failed." << std::endl;
            output = false;
        }
        if(!matrix1.getRangeMap()->isSameAs(*matrix2.getRangeMap()))
        {
            if(output)
                out << std::endl;
            out << "Compare: RangeMap check failed." << std::endl;
            out << matrix1.getRangeMap().getRawPtr() << " : " << matrix2.getRangeMap().getRawPtr() << std::endl;
            output = false;
        }
        if(!matrix1.getRowMap()->isSameAs(*matrix2.getRowMap()))
        {
            if(output)
                out << std::endl;
            out << "Compare: RowMap check failed." << std::endl;
            output = false;
        }
        if(!matrix1.getColMap()->isSameAs(*matrix2.getColMap()))
        {
            if(output)
                out << std::endl;
            out << "Compare: ColMap check failed." << std::endl;
            output = false;
        }
        if(!matrix1.getDomainMap()->isSameAs(*matrix2.getDomainMap()))
        {
            if(output)
                out << std::endl;
            out << "Compare: DomainMap check failed." << std::endl;
            output = false;
        }

        auto lclmtx1 = matrix1.getLocalMatrixHost();
        auto lclmtx2 = matrix2.getLocalMatrixHost();

        auto rowptr1 = lclmtx1.graph.row_map;
        auto rowptr2 = lclmtx2.graph.row_map;

        auto colind1 = lclmtx1.graph.entries;
        auto colind2 = lclmtx2.graph.entries;

        auto rbo1 = lclmtx1.graph.row_block_offsets;
        auto rbo2 = lclmtx2.graph.row_block_offsets;

        if(rowptr1.extent(0) != rowptr2.extent(0))
        {
            if(output)
                out << std::endl;
            out << "Compare: Kokkos::StaticCrsGraph::rowptr extent check failed." << std::endl;
            output = false;
        }
        if(colind1.extent(0) != colind2.extent(0))
        {
            if(output)
                out << std::endl;
            out << "Compare: Kokkos::StaticCrsGraph::colind extent check failed." << std::endl;
            output = false;
        }
        if(rbo1.extent(0) != rbo2.extent(0))
        {
            if(output)
                out << std::endl;
            out << "Compare: Kokkos::StaticCrsGraph::row_block_offsets extents check failed." << std::endl;
            output = false;
        }

        bool success = true;
        out << std::endl;
        TEST_COMPARE_ARRAYS(rowptr1, rowptr2);
        if(!success)
        {
            out << "Compare: Kokkos::StaticCrsGraph::rowptr match failed." << std::endl;
            output = false;
        }

        TEST_COMPARE_ARRAYS(colind1, colind2);
        if(!success)
        {
            out << "Compare: Kokkos::StaticCrsGraph::colind match failed." << std::endl;
            output = false;
        }

        TEST_COMPARE_ARRAYS(rbo1, rbo2);
        if(!success)
        {
            out << "Compare: Kokkos::StaticCrsGraph::row_block_offsets match failed." << std::endl;
            output = false;
        }

        return output;
    }

};      // class crsMatrix_Swap_Tester

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
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsGraph, Swap, LO, GO, Scalar, Node)
{
    auto crsMatrixTester = Tpetra::crsMatrix_Swap_Tester<Scalar, LO, GO, Node>();
    crsMatrixTester.execute(out, success);
}


//
// INSTANTIATIONS
//



// Tests to build and run in both debug and release modes.  We will
// instantiate them over all enabled local ordinal (LO), global
// ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP_DEBUG_AND_RELEASE(SCALAR, LO, GO, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsGraph, Swap, LO, GO, SCALAR, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP_DEBUG_AND_RELEASE)

}      // namespace
