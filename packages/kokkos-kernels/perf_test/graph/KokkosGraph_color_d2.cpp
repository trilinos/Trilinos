/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

// EXERCISE 1 Goal:
//   Use Kokkos to parallelize the outer loop of <y,Ax> using Kokkos::parallel_reduce.
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <string>
#include <sys/time.h>

#include <Kokkos_Core.hpp>

#include <KokkosKernels_IOUtils.hpp>
#include <KokkosKernels_MyCRSMatrix.hpp>
#include <KokkosKernels_TestParameters.hpp>
#include <KokkosGraph_Distance2Color.hpp>


using namespace KokkosGraph;

#ifdef KOKKOSKERNELS_INST_DOUBLE
    typedef double kk_scalar_t;
#else
    #ifdef KOKKOSKERNELS_INST_FLOAT
        typedef float kk_scalar_t;
    #endif
#endif

#ifdef KOKKOSKERNELS_INST_OFFSET_INT
    typedef int kk_size_type;
#else
    #ifdef KOKKOSKERNELS_INST_OFFSET_SIZE_T
        typedef size_t kk_size_type;
    #endif
#endif

#ifdef KOKKOSKERNELS_INST_ORDINAL_INT
    typedef int kk_lno_t;
#else
    #ifdef KOKKOSKERNELS_INST_ORDINAL_INT64_T
        typedef int64_t kk_lno_t;
    #endif
#endif

// Toggle EXPERIMENTAL code for calculating
// a tight bound on distance-2 degree.
#define USE_EXPERIMENTAL_MAXD2DEGREE 0


using namespace KokkosGraph;



void print_options(std::ostream &os, const char *app_name, unsigned int indent = 0)
{
    std::string spaces(indent, ' ');
    os << "Usage:" << std::endl
       << spaces << "  " << app_name << " [parameters]" << std::endl
       << std::endl
       << spaces << "Parameters:" << std::endl
       << spaces << "  Parallelism (select one of the following):" << std::endl
       << spaces << "      --serial <N>        Execute serially." << std::endl
       << spaces << "      --threads <N>       Use N posix threads." << std::endl
       << spaces << "      --openmp <N>        Use OpenMP with N threads." << std::endl
       << spaces << "      --cuda              Use CUDA" << std::endl
       << std::endl
       << spaces << "  Required Parameters:" << std::endl
       << spaces << "      --amtx <filename>   Input file in Matrix Market format (.mtx)." << std::endl
       << std::endl
       << spaces << "      --algorithm <algorithm_name>   Set the algorithm to use.  Allowable values are:" << std::endl
       << spaces << "                 COLORING_D2_MATRIX_SQUARED  - Matrix-squared + Distance-1 method." << std::endl
       << spaces << "                 COLORING_D2_SERIAL          - Serial algorithm (must use with 'serial' mode)" << std::endl
       << spaces << "                 COLORING_D2_VB              - Vertex Based method using boolean forbidden array (Default)." << std::endl
       << spaces << "                 COLORING_D2_VB_BIT          - VB with Bitvector Forbidden Array" << std::endl
       << spaces << "                 COLORING_D2_VB_BIT_EF       - VB_BIT with Edge Filtering" << std::endl

       << std::endl
       << spaces << "  Optional Parameters:" << std::endl
       << spaces << "      --repeat <N>        Set number of test repetitions (Default: 1) " << std::endl
       << spaces << "      --verbose           Enable verbose mode (record and print timing + extra information)" << std::endl
       << spaces << "      --chunksize <N>     Set the chunk size." << std::endl
       << spaces << "      --dynamic           Use dynamic scheduling." << std::endl
       << spaces << "      --teamsize  <N>     Set the team size." << std::endl
       << spaces << "      --vectorsize <N>    Set the vector size." << std::endl
       << spaces << "      --help              Print out command line help." << std::endl
       << spaces << " " << std::endl;
}


int parse_inputs(KokkosKernels::Experiment::Parameters &params, int argc, char **argv)
{
    bool got_required_param_amtx      = false;
    bool got_required_param_algorithm = false;

    for(int i = 1; i < argc; ++i)
    {
        if(0 == strcasecmp(argv[i], "--threads"))
        {
            params.use_threads = atoi(argv[++i]);
        }
        else if(0 == strcasecmp(argv[i], "--serial"))
        {
            params.use_serial = atoi(argv[++i]);
        }
        else if(0 == strcasecmp(argv[i], "--openmp"))
        {
            params.use_openmp = atoi(argv[++i]);
            std::cout << "use_openmp = " << params.use_openmp << std::endl;
        }
        else if(0 == strcasecmp(argv[i], "--cuda"))
        {
            params.use_cuda = 1;
        }
        else if(0 == strcasecmp(argv[i], "--repeat"))
        {
            params.repeat = atoi(argv[++i]);
        }
        else if(0 == strcasecmp(argv[i], "--chunksize"))
        {
            params.chunk_size = atoi(argv[++i]);
        }
        else if(0 == strcasecmp(argv[i], "--teamsize"))
        {
            params.team_size = atoi(argv[++i]);
        }
        else if(0 == strcasecmp(argv[i], "--vectorsize"))
        {
            params.vector_size = atoi(argv[++i]);
        }
        else if(0 == strcasecmp(argv[i], "--amtx"))
        {
            got_required_param_amtx = true;
            params.a_mtx_bin_file   = argv[++i];
        }
        else if(0 == strcasecmp(argv[i], "--dynamic"))
        {
            params.use_dynamic_scheduling = 1;
        }
        else if(0 == strcasecmp(argv[i], "--verbose"))
        {
            params.verbose = 1;
        }
        else if(0 == strcasecmp(argv[i], "--algorithm"))
        {
            ++i;
            if(0 == strcasecmp(argv[i], "COLORING_D2_MATRIX_SQUARED"))
            {
                params.algorithm             = 1;
                got_required_param_algorithm = true;
            }
            else if(0 == strcasecmp(argv[i], "COLORING_D2_SERIAL"))
            {
                params.algorithm             = 2;
                got_required_param_algorithm = true;
            }
            else if(0 == strcasecmp(argv[i], "COLORING_D2_VB") || 0 == strcasecmp(argv[i], "COLORING_D2") )
            {
                params.algorithm             = 3;
                got_required_param_algorithm = true;
            }
            else if(0 == strcasecmp(argv[i], "COLORING_D2_VB_BIT"))
            {
                params.algorithm             = 4;
                got_required_param_algorithm = true;
            }
            else if(0 == strcasecmp(argv[i], "COLORING_D2_VB_BIT_EF"))
            {
                params.algorithm             = 5;
                got_required_param_algorithm = true;
            }
            else
            {
                std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
                print_options(std::cout, argv[0]);
                return 1;
            }
        }
        else if(0 == strcasecmp(argv[i], "--help") || 0 == strcasecmp(argv[i], "-h"))
        {
            print_options(std::cout, argv[0]);
            return 1;
        }
        else
        {
            std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
            print_options(std::cout, argv[0]);
            return 1;
        }
    }

    if(!got_required_param_amtx)
    {
        std::cout << "Missing required parameter amtx" << std::endl << std::endl;
        print_options(std::cout, argv[0]);
        return 1;
    }
    if(!got_required_param_algorithm)
    {
        std::cout << "Missing required parameter algorithm" << std::endl << std::endl;
        print_options(std::cout, argv[0]);
        return 1;
    }
    if(!params.use_serial && !params.use_threads && !params.use_openmp && !params.use_cuda)
    {
        print_options(std::cout, argv[0]);
        return 1;
    }
    return 0;
}


namespace KokkosKernels {
namespace Experiment {


std::string getCurrentDateTimeStr()
{
    // Note: This could be replaced with `std::put_time(&tm, "%FT%T%z")` but std::put_time isn't
    //       supported on the intel C++ compilers as of v. 17.0.x
    time_t now = time(0);
    char output[100];
    std::strftime(output, sizeof(output), "%FT%T%Z", std::localtime(&now));
    return output;
}


template<typename ExecSpace, typename crsGraph_t, typename crsGraph_t2, typename crsGraph_t3, typename TempMemSpace, typename PersistentMemSpace>
void run_experiment(crsGraph_t crsGraph, Parameters params)
{
    using namespace KokkosGraph;
    using namespace KokkosGraph::Experimental;

    int algorithm  = params.algorithm;
    int repeat     = params.repeat;
    int chunk_size = params.chunk_size;

    int shmemsize              = params.shmemsize;
    int team_size              = params.team_size;
    int use_dynamic_scheduling = params.use_dynamic_scheduling;
    int verbose                = params.verbose;

    // char spgemm_step = params.spgemm_step;
    int vector_size = params.vector_size;

    using lno_view_t = typename crsGraph_t3::row_map_type::non_const_type;
    using lno_nnz_view_t = typename crsGraph_t3::entries_type::non_const_type;

    using size_type = typename lno_view_t::non_const_value_type;
    using lno_t     = typename lno_nnz_view_t::non_const_value_type;

    typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, kk_scalar_t, ExecSpace, TempMemSpace, PersistentMemSpace> KernelHandle;

    // Get Date/Time stamps of start to use later when printing out summary data.
    //auto t  =  std::time(nullptr);
    //auto tm = *std::localtime(&t);

    // Note: crsGraph.numRows() == number of vertices in the 'graph'
    //       crsGraph.entries.extent(0) == number of edges in the 'graph'
    std::cout << "Num verts: " << crsGraph.numRows()         << std::endl
              << "Num edges: " << crsGraph.entries.extent(0) << std::endl;

    KernelHandle kh;
    kh.set_team_work_size(chunk_size);
    kh.set_shmem_size(shmemsize);
    kh.set_suggested_team_size(team_size);
    kh.set_suggested_vector_size(vector_size);

    if(use_dynamic_scheduling)
    {
        kh.set_dynamic_scheduling(true);
    }

    if(verbose)
    {
        kh.set_verbose(true);
    }

    // accumulators for average stats
    size_t total_colors = 0;
    size_t total_phases = 0;

    std::string label_algorithm;
    switch(algorithm)
    {
        case 1:
            kh.create_distance2_graph_coloring_handle(COLORING_D2_MATRIX_SQUARED);
            label_algorithm = "COLORING_D2_MATRIX_SQUARED";
            break;
        case 2:
            kh.create_distance2_graph_coloring_handle(COLORING_D2_SERIAL);
            label_algorithm = "COLORING_D2_SERIAL";
            break;
        case 3:
            kh.create_distance2_graph_coloring_handle(COLORING_D2_VB);
            label_algorithm = "COLORING_D2_VB";
            break;
        case 4:
            kh.create_distance2_graph_coloring_handle(COLORING_D2_VB_BIT);
            label_algorithm = "COLORING_D2_VB_BIT";
            break;
        case 5:
            kh.create_distance2_graph_coloring_handle(COLORING_D2_VB_BIT_EF);
            label_algorithm = "COLORING_D2_VB_BIT_EF";
            break;
        default:
            kh.create_distance2_graph_coloring_handle(COLORING_D2_VB);
            label_algorithm = "COLORING_D2_VB";
            break;
    }

    std::cout << std::endl << "Run Graph Color D2 (" << label_algorithm << ")" << std::endl;

    // If any of the runs have an invalid result, this will be set to false.
    bool all_results_valid = true;

    // Loop over # of experiments to run
    for(int i = 0; i < repeat; ++i)
    {
        graph_compute_distance2_color(&kh, crsGraph.numRows(), crsGraph.numCols(), crsGraph.row_map, crsGraph.entries, crsGraph.row_map, crsGraph.entries);

        total_colors += kh.get_distance2_graph_coloring_handle()->get_num_colors();
        total_phases += kh.get_distance2_graph_coloring_handle()->get_num_phases();

        std::cout << "Total Time: " << kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time() << std::endl
                  << "Num colors: " << kh.get_distance2_graph_coloring_handle()->get_num_colors() << std::endl
                  << "Num Phases: " << kh.get_distance2_graph_coloring_handle()->get_num_phases() << std::endl;

        std::cout << "\t";
        KokkosKernels::Impl::print_1Dview(kh.get_distance2_graph_coloring_handle()->get_vertex_colors());
        std::cout << std::endl;

        // If verbose mode is on and there the graph has fewer than 1500 verts, dump a GraphVIZ DOT file.
        if(verbose && repeat==i+1 && crsGraph.numRows() < 1500)
        {
            auto colors = kh.get_distance2_graph_coloring_handle()->get_vertex_colors();
            std::ofstream os("G.dot", std::ofstream::out);
            kh.get_distance2_graph_coloring_handle()->dump_graphviz(os, crsGraph.numRows(), crsGraph.row_map, crsGraph.entries, colors);
        }

        // ------------------------------------------
        // Verify correctness
        // ------------------------------------------
        bool d2_coloring_is_valid            = false;
        bool d2_coloring_validation_flags[4] = { false };

        d2_coloring_is_valid = KokkosGraph::Impl::graph_verify_distance2_color(&kh, crsGraph.numRows(), crsGraph.numCols(), crsGraph.row_map, crsGraph.entries, crsGraph.row_map, crsGraph.entries, d2_coloring_validation_flags);

        // Print out messages based on coloring validation check.
        if(d2_coloring_is_valid)
        {
            std::cout << std::endl << "Distance-2 Graph Coloring is VALID" << std::endl << std::endl;
        }
        else
        {
            all_results_valid = false;
            std::cout << std::endl
                      << "Distance-2 Graph Coloring is NOT VALID" << std::endl
                      << "  - Vert(s) left uncolored : " << d2_coloring_validation_flags[1] << std::endl
                      << "  - Invalid D2 Coloring    : " << d2_coloring_validation_flags[2] << std::endl
                      << std::endl;
        }
        if(d2_coloring_validation_flags[3])
        {
            std::cout << "Distance-2 Graph Coloring may have poor quality." << std::endl
                      << "  - Vert(s) have high color value : " << d2_coloring_validation_flags[3] << std::endl
                      << std::endl;
        }

        // ------------------------------------------
        // Print out the colors histogram
        // ------------------------------------------
        KokkosGraph::Impl::graph_print_distance2_color_histogram(&kh, crsGraph.numRows(), crsGraph.numCols(), crsGraph.row_map, crsGraph.entries, crsGraph.row_map, crsGraph.entries, false);

    } // for i...

    // ------------------------------------------
    // Compute Distance 2 Degree Stats
    // ------------------------------------------
    std::cout << "Compute Distance-2 Degree " << std::endl;

    Kokkos::Impl::Timer timer;

    #if defined(USE_EXPERIMENTAL_MAXD2DEGREE) && USE_EXPERIMENTAL_MAXD2DEGREE
    double time_d2_degree;
    timer.reset();

    typedef typename KernelHandle::GraphColoringHandleType::non_const_1d_size_type_view_t non_const_1d_size_type_view_t;
    non_const_1d_size_type_view_t degree_d2_dist = non_const_1d_size_type_view_t("degree d2", crsGraph.numRows());

    size_t degree_d2_max=0;
    KokkosGraph::Impl::graph_compute_distance2_degree(&kh, crsGraph.numRows(), crsGraph.numCols(),
                                   crsGraph.row_map, crsGraph.entries,
                                   crsGraph.row_map, crsGraph.entries,
                                   degree_d2_dist, degree_d2_max);
    time_d2_degree = timer.seconds();
    #endif

    double total_time                   = kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time();
    double total_time_color_greedy      = kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time_phase1();
    double total_time_find_conflicts    = kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time_phase2();
    double total_time_resolve_conflicts = kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time_phase3();
    double total_time_matrix_squared    = kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time_phase4();
    double total_time_matrix_squared_d1 = kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time_phase5();

    double avg_time                   = total_time / (double)repeat;
    double avg_time_color_greedy      = total_time_color_greedy / (double)repeat;
    double avg_time_find_conflicts    = total_time_find_conflicts / (double)repeat;
    double avg_time_resolve_conflicts = total_time_resolve_conflicts / (double)repeat;
    double avg_colors                 = total_colors / (double)repeat;
    double avg_phases                 = total_phases / (double)repeat;
    double avg_time_matrix_squared    = total_time_matrix_squared / (double)repeat;
    double avg_time_matrix_squared_d1 = total_time_matrix_squared_d1 / (double)repeat;

    std::string a_mtx_bin_file = params.a_mtx_bin_file;
    a_mtx_bin_file             = a_mtx_bin_file.substr(a_mtx_bin_file.find_last_of("/\\") + 1);


    int result;
    char hostname[100];
    char username[100];

    result = gethostname(hostname, 100);
    if(result)
    {
        perror("gethostname");
    }

    result = getlogin_r(username, 100);
    if(result)
    {
        perror("getlogin_r");
    }

    std::string all_results_valid_str = "PASSED";
    if(!all_results_valid)
        all_results_valid_str = "FAILED";

    std::string currentDateTimeStr = getCurrentDateTimeStr();

    std::cout << "Summary" << std::endl
              << "-------" << std::endl
              << "    Date/Time      : " << currentDateTimeStr << std::endl
              << "    KExecSName     : " << Kokkos::DefaultExecutionSpace::name() << std::endl
              << "    Filename       : " << a_mtx_bin_file << std::endl
              << "    Num Verts      : " << crsGraph.numRows() << std::endl
              << "    Num Edges      : " << crsGraph.entries.extent(0) << std::endl
              << "    Concurrency    : " << Kokkos::DefaultExecutionSpace::concurrency() << std::endl
              << "    Algorithm      : " << label_algorithm << std::endl
              #if defined(USE_EXPERIMENTAL_MAXD2DEGREE) && USE_EXPERIMENTAL_MAXD2DEGREE
              << "Graph Stats" << std::endl
              << "    Degree D2 Max  : " << degree_d2_max << std::endl
              << "    Degree D2 Time : " << time_d2_degree << std::endl
              #endif
              << "Overall Time/Stats" << std::endl
              << "    Total Time     : " << total_time << std::endl
              << "    Avg Time       : " << avg_time << std::endl
              << "VB Distance[1|2] Stats " << std::endl
              << "    Avg Time CG    : " << avg_time_color_greedy << std::endl
              << "    Avg Time FC    : " << avg_time_find_conflicts << std::endl
              << "    Avg Time RC    : " << avg_time_resolve_conflicts << std::endl
              << "Matrix-Squared + D1 Stats" << std::endl
              << "    Avg Time to M^2: " << avg_time_matrix_squared << std::endl
              << "    Avg Time to D1 : " << avg_time_matrix_squared_d1 << std::endl
              << "Coloring Stats" << std::endl
              << "    Avg colors     : " << avg_colors << std::endl
              << "    Avg Phases     : " << avg_phases << std::endl
              << "    Validation     : " << all_results_valid_str << std::endl
              << std::endl;

    std::cout << "CSVTIMEHDR"
              << "," << "Filename"
              << "," << "Host"
              << "," << "DateTime"
              << "," << "Num Rows"
              << "," << "Num Edges"
              << "," << "Execution Space"
              << "," << "Algorithm"
              << "," << "Concurrency"
              << "," << "Repetitions"
              << "," << "Total Time"
              << "," << "Total Time to M^2"
              << "," << "Total Time D1(M^2)"
              << "," << "Total Time CG"
              << "," << "Total Time FC"
              << "," << "Total Time RC"
              << "," << "Avg Colors"
              << "," << "Avg Num Phases"
              #if defined(USE_EXPERIMENTAL_MAXD2DEGREE) && USE_EXPERIMENTAL_MAXD2DEGREE
              << "," << "Time D2 Degree"
              << "," << "Degree D2 Max"
              #endif
              << "," << "Validation"
              << std::endl;

    std::cout << "CSVTIMEDATA"
              << "," << a_mtx_bin_file
              << "," << hostname
              << "," << currentDateTimeStr
              << "," << crsGraph.numRows()
              << "," << crsGraph.entries.extent(0)
              << "," << Kokkos::DefaultExecutionSpace::name()
              << "," << label_algorithm
              << "," << Kokkos::DefaultExecutionSpace::concurrency()
              << "," << repeat
              << "," << total_time
              << "," << total_time_matrix_squared
              << "," << total_time_matrix_squared_d1
              << "," << total_time_color_greedy
              << "," << total_time_find_conflicts
              << "," << total_time_resolve_conflicts

              << "," << avg_colors
              << "," << avg_phases
              #if defined(USE_EXPERIMENTAL_MAXD2DEGREE) && USE_EXPERIMENTAL_MAXD2DEGREE
              << "," << time_d2_degree
              << "," << degree_d2_max
              #endif
              << "," << all_results_valid_str
              << std::endl;

    std::cout << "CSVHISTHDR"
              << "," << "Filename"
              << "," << "Host"
              << "," << "DateTime"
              << "," << "Num Rows"
              << "," << "Num Edges"
              << "," << "Execution Space"
              << "," << "Algorithm"
              << "," << "Concurrency"
              << "," << "Histogram: 1 .. N"
              << std::endl;

    std::cout << "CSVHISTDATA"
              << "," << a_mtx_bin_file
              << "," << hostname
              << "," << currentDateTimeStr
              << "," << crsGraph.numRows()
              << "," << crsGraph.entries.extent(0)
              << "," << Kokkos::DefaultExecutionSpace::name()
              << "," << label_algorithm
              << "," << Kokkos::DefaultExecutionSpace::concurrency()
              << ",";
    KokkosGraph::Impl::graph_print_distance2_color_histogram(&kh, crsGraph.numRows(), crsGraph.numCols(), crsGraph.row_map, crsGraph.entries, crsGraph.row_map, crsGraph.entries, true);
    std::cout << std::endl;

    // Kokkos::print_configuration(std::cout);
}


template<typename size_type, typename lno_t, typename exec_space, typename hbm_mem_space>
void experiment_driver(Parameters params)
{
    using myExecSpace     = exec_space;
    using myFastDevice    = Kokkos::Device<exec_space, hbm_mem_space>;
    using fast_crstmat_t  = typename MyKokkosSparse::CrsMatrix<double, lno_t, myFastDevice, void, size_type>;
    using fast_graph_t    = typename fast_crstmat_t::StaticCrsGraphType;

    char *a_mat_file = params.a_mtx_bin_file;

    fast_graph_t a_fast_crsgraph, /*b_fast_crsgraph,*/ c_fast_crsgraph;

    if(params.a_mem_space == 1)
    {
        fast_crstmat_t a_fast_crsmat;
        a_fast_crsmat            = KokkosKernels::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
        a_fast_crsgraph          = a_fast_crsmat.graph;
        a_fast_crsgraph.num_cols = a_fast_crsmat.numCols();
    }

    if(params.a_mem_space == 1 && params.b_mem_space==1 && params.c_mem_space==1 && params.work_mem_space==1)
    {
        KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, fast_graph_t, fast_graph_t, hbm_mem_space, hbm_mem_space>
                (a_fast_crsgraph, /*b_fast_crsgraph,*/ params);
    }
    else
    {
        std::cout << ">>> unhandled memspace configuration flags:"      << std::endl
                  << ">>>   a_mem_space    = " << params.a_mem_space    << std::endl
                  << ">>>   b_mem_space    = " << params.a_mem_space    << std::endl
                  << ">>>   c_mem_space    = " << params.a_mem_space    << std::endl
                  << ">>>   work_mem_space = " << params.work_mem_space << std::endl;
    }
}


}      // namespace Experiment
}      // namespace KokkosKernels



int main(int argc, char *argv[])
{
    KokkosKernels::Experiment::Parameters params;

    // Override default repeats (default is 6)
    params.repeat = 1;

    if(parse_inputs(params, argc, argv))
    {
        return 1;
    }

    if(params.a_mtx_bin_file == NULL)
    {
        std::cerr << "Provide a matrix file" << std::endl;
        return 0;
    }

    std::cout << "Sizeof(kk_lno_t) : " << sizeof(kk_lno_t)     << std::endl
              << "Sizeof(size_type): " << sizeof(kk_size_type) << std::endl;

    const int num_threads = params.use_openmp;      // Assumption is that use_openmp variable is provided as number of threads
    const int device_id   = 0;
    Kokkos::initialize(Kokkos::InitArguments(num_threads, -1, device_id));

    // Print out verbose information about the configuration of the run.
    // Kokkos::print_configuration(std::cout);

    #if defined(KOKKOS_MULTI_MEM)
    const bool use_multi_mem = true;
    // todo: Throw an error or print a message if KOKKOS_MULTI_MEM is enabled for this test?  (WCMCLEN--SCAFFOLDING)
    #else
    const bool use_multi_mem = false;
    #endif

    #if defined(KOKKOS_ENABLE_OPENMP)
    if(params.use_openmp)
    {
        if(!use_multi_mem)
        {
            KokkosKernels::Experiment::experiment_driver<kk_size_type, kk_lno_t, Kokkos::OpenMP, Kokkos::OpenMP::memory_space>(params);
        }
    }
    #endif

    #if defined(KOKKOS_ENABLE_THREADS)
    if(params.use_threads)
    {
        if(!use_multi_mem)
        {
            KokkosKernels::Experiment::experiment_driver<kk_size_type, kk_lno_t, Kokkos::Threads, Kokkos::Threads::memory_space>(params);
        }
    }
    #endif

    #if defined(KOKKOS_ENABLE_CUDA)
    if(params.use_cuda)
    {
        if(!use_multi_mem)
        {
            KokkosKernels::Experiment::experiment_driver<kk_size_type, kk_lno_t, Kokkos::Cuda, Kokkos::Cuda::memory_space>(params);
        }
    }
    #endif

    #if defined(KOKKOS_ENABLE_SERIAL)
    if(params.use_serial)
    {
        if(!use_multi_mem)
        {
            KokkosKernels::Experiment::experiment_driver<kk_size_type, kk_lno_t, Kokkos::Serial, Kokkos::Serial::memory_space>(params);
        }
    }
    #endif

    Kokkos::finalize();

    return 0;
}
