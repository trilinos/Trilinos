/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact William McLendon (wcmclen@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef _KOKKOSGRAPH_DISTANCE2COLOR_IMPL_HPP
#define _KOKKOSGRAPH_DISTANCE2COLOR_IMPL_HPP

#include <iomanip>
#include <stdexcept>
#include <vector>

#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_UniqueToken.hpp>

#include <KokkosKernels_Uniform_Initialized_MemoryPool.hpp>
#include <KokkosKernels_HashmapAccumulator.hpp>
#include <KokkosSparse_spgemm.hpp>

#include <impl/Kokkos_Timer.hpp>

#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosGraph_Distance1ColorHandle.hpp"      // todo: remove this  (SCAFFOLDING - WCMCLEN)
#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosKernels_Handle.hpp"


namespace KokkosGraph {

namespace Impl {

#define VB_D2_COLORING_FORBIDDEN_SIZE 64
#define VBBIT_D2_COLORING_FORBIDDEN_SIZE 64



/*!
 * \brief Distance-2 Graph Coloring class
 *
 * This class supports direct methods for distance-2 graph coloring.
 *
 */
template<typename HandleType,
         typename lno_row_view_t_,
         typename lno_nnz_view_t_,
         typename clno_row_view_t_,
         typename clno_nnz_view_t_>
class GraphColorDistance2
{
    // static_assert( std::is_same< HandleType, KokkosGraph::GraphColorDistance2Handle>::value, "Incorrect HandleType for this
    // class.");

  public:
    using in_lno_row_view_type             = lno_row_view_t_;
    using in_lno_nnz_view_type             = lno_nnz_view_t_;
    using color_view_type                  = typename HandleType::color_view_type;
    using color_type                       = typename HandleType::color_type;
    using size_type                        = typename HandleType::size_type;
    using nnz_lno_type                     = typename HandleType::nnz_lno_type;
    using row_lno_host_view_type           = typename in_lno_row_view_type::HostMirror;
    using nnz_lno_host_view_type           = typename in_lno_nnz_view_type::HostMirror;
    using color_host_view_type             = typename HandleType::color_host_view_type;
    using my_exec_space                    = typename HandleType::HandleExecSpace;
    using MyTempMemorySpace                = typename HandleType::HandleTempMemorySpace;
    using const_size_type                  = typename HandleType::const_size_type;
    using row_lno_view_device_type         = typename lno_row_view_t_::device_type;
    using const_lno_row_view_type          = typename lno_row_view_t_::const_type;
    using const_lno_nnz_view_t             = typename lno_nnz_view_t_::const_type;
    using non_const_lno_nnz_view_type      = typename lno_nnz_view_t_::non_const_type;
    using const_clno_row_view_t            = typename clno_row_view_t_::const_type;
    using const_clno_nnz_view_t            = typename clno_nnz_view_t_::const_type;
    using non_const_clno_nnz_view_t        = typename clno_nnz_view_t_::non_const_type;
    using nnz_lno_temp_work_view_t         = typename HandleType::nnz_lno_temp_work_view_type;
    using single_dim_index_view_type       = typename Kokkos::View<nnz_lno_type, row_lno_view_device_type>;
    using range_policy_type                = Kokkos::RangePolicy<my_exec_space>;
    using team_policy_type                 = Kokkos::TeamPolicy<my_exec_space>;
    using team_member_type                 = typename team_policy_type::member_type;
    using non_const_1d_bool_view_t         = Kokkos::View<bool*>;
    using non_const_1d_size_type_view_type = typename HandleType::non_const_1d_size_type_view_type;
    using bit_64_forbidden_type            = long long int;

    // For HashmapAccumulator (EXPERIMENTAL)
    using pool_memory_space_t = typename KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, nnz_lno_type>;


  protected:
    nnz_lno_type            nv;          // num vertices
    nnz_lno_type            nr;          // num_rows  (# verts)
    nnz_lno_type            nc;          // num cols
    size_type               ne;          // # edges
    const_lno_row_view_type xadj;        // rowmap, transpose of rowmap
    const_lno_nnz_view_t    adj;         // entries, transpose of entries   (size = # edges)
    const_clno_row_view_t   t_xadj;      // rowmap, transpose of rowmap
    const_clno_nnz_view_t   t_adj;       // entries, transpose of entries

    HandleType* gc_handle;      // pointer to the graph coloring handle


  private:
    int  _chunkSize;      // the size of the minimum work unit assigned to threads.  Changes the convergence on GPUs
    int  _max_num_iterations;
    char _use_color_set;      // The VB Algorithm Type: 0: VB,  1: VBCS,  2: VBBIT  (1, 2 not implemented).
    bool _ticToc;             // if true print info in each step
    bool _verbose;            // if true, print verbose information


  public:
    /**
     * \brief GraphColorDistance2 constructor.
     * \param nv_: number of vertices in the graph
     * \param ne_: number of edges in the graph
     * \param row_map: the xadj array of the graph. Its size is nv_ +1
     * \param entries: adjacency array of the graph. Its size is ne_
     * \param handle: GraphColoringHandle object that holds the specification about the graph coloring,
     *    including parameters.
     */
    GraphColorDistance2(nnz_lno_type            nv_,
                        nnz_lno_type            nc_,
                        size_type               ne_,
                        const_lno_row_view_type row_map,
                        const_lno_nnz_view_t    entries,
                        const_clno_row_view_t   t_row_map,
                        const_clno_nnz_view_t   t_entries,
                        HandleType*             handle)
        : nv(nv_)
        , nr(nv_)
        , nc(nc_)
        , ne(ne_)
        , xadj(row_map)
        , adj(entries)
        , t_xadj(t_row_map)
        , t_adj(t_entries)
        , gc_handle(handle)
        , _chunkSize(handle->get_vb_chunk_size())
        , _max_num_iterations(handle->get_max_number_of_iterations())
        , _use_color_set(0)
        , _ticToc(handle->get_verbose())
        , _verbose(handle->get_verbose())
    {
    }


    /**
     *  \brief GraphColor destructor.
     */
    virtual ~GraphColorDistance2() {}




    // -----------------------------------------------------------------
    //
    // GraphColorDistance2::execute()
    //
    // -----------------------------------------------------------------

    /**
     * \brief Computes the distance-2 graph coloring.
     */
    virtual void compute_distance2_color()
    {
        color_view_type colors_out("Graph Colors", this->nv);

        // If the selected pass is using edge filtering we should set the flag.
        bool using_edge_filtering = false;
        if(COLORING_D2_VB_BIT_EF == this->gc_handle->get_coloring_algo_type())
        {
            using_edge_filtering = true;
        }

        // Data:
        // gc_handle = graph coloring handle
        // nr        = num_rows  (scalar)
        // nc        = num_cols  (scalar)
        // xadj      = row_map   (view 1 dimension - [num_verts+1] - entries index into adj )
        // adj       = entries   (view 1 dimension - [num_edges]   - adjacency list )
        if(this->_ticToc)
        {
            std::cout << "\tcolor_graph_d2 params:" << std::endl
                      << "\t  algorithm                : " << (int)this->_use_color_set << std::endl
                      << "\t  ticToc                   : " << this->_ticToc << std::endl
                      << "\t  max_num_iterations       : " << this->_max_num_iterations << std::endl
                      << "\t  chunkSize                : " << this->_chunkSize << std::endl
                      << "\t  use_color_set            : " << (int)this->_use_color_set << std::endl
                      << "\t  Edge Filtering Pass?     : " << (int)using_edge_filtering << std::endl
                      << "\tgraph information:" << std::endl
                      << "\t  nv                       : " << this->nv << std::endl
                      << "\t  ne                       : " << this->ne << std::endl;
            /*
            // For extra verbosity if debugging...
            prettyPrint1DView(this->xadj,   ">>> xadj   ", 500);
            prettyPrint1DView(this->adj,    ">>>  adj   ", 500);
            prettyPrint1DView(this->t_xadj, ">>> t_xadj ", 500);
            prettyPrint1DView(this->t_adj,  ">>> t_adj  ", 500);
            */
        }

        #if 0
        // (EXPERIMENTAL) Distance-2 Degree calculation
        // Compute Distance-2 Degree of the vertices.
        //  -- Keeping this around in case we need to use it later on as an example.
        size_t degree_d2_max=0;
        size_t degree_d2_sum=0;
        non_const_1d_size_type_view_type degree_d2 = non_const_1d_size_type_view_type("degree d2", this->nv);
        this->compute_distance2_degree(degree_d2, degree_d2_max, degree_d2_sum);
        if(true || _ticToc )
        {
            prettyPrint1DView(degree_d2, "Degree D2", 150);
            std::cout << ">>> Max D2 Degree: " << degree_d2_max << std::endl;
            std::cout << ">>> D2 Degree Sum: " << degree_d2_sum<< std::endl;
        }
        #endif

        // conflictlist - store conflicts that can happen when we're coloring in parallel.
        nnz_lno_temp_work_view_t current_vertexList =
          nnz_lno_temp_work_view_t(Kokkos::ViewAllocateWithoutInitializing("vertexList"), this->nv);

        // init conflictlist sequentially.
        Kokkos::parallel_for("InitList", range_policy_type(0, this->nv), functorInitList<nnz_lno_temp_work_view_t>(current_vertexList));

        // Next iteratons's conflictList
        nnz_lno_temp_work_view_t next_iteration_recolorList;

        // Size the next iteration conflictList
        single_dim_index_view_type next_iteration_recolorListLength;

        // Vertices to recolor.  Will swap with vertexList
        next_iteration_recolorList = nnz_lno_temp_work_view_t(Kokkos::ViewAllocateWithoutInitializing("recolorList"), this->nv);
        next_iteration_recolorListLength = single_dim_index_view_type("recolorListLength");

        nnz_lno_type numUncolored             = this->nv;
        nnz_lno_type numUncoloredPreviousIter = this->nv + 1;
        nnz_lno_type current_vertexListLength = this->nv;

        double              time;
        double              total_time = 0.0;
        Kokkos::Impl::Timer timer;

        int iter = 0;
        for(; (iter < _max_num_iterations) && (numUncolored > 0); iter++)
        {
            timer.reset();

            // Save the # of uncolored from the previous iteration
            numUncoloredPreviousIter = numUncolored;

            // ------------------------------------------
            // Do greedy color
            // ------------------------------------------
            if(using_edge_filtering)
            {
                // Temp copy of the adj array (mutable)
                // - This is required for edge-filtering passes to avoid
                //   side effects since edge filtering modifies the adj array.
                nnz_lno_temp_work_view_t adj_copy;
                adj_copy = nnz_lno_temp_work_view_t(Kokkos::ViewAllocateWithoutInitializing("adj copy"), this->ne);
                Kokkos::deep_copy(adj_copy, this->adj);

                non_const_clno_nnz_view_t t_adj_copy;
                t_adj_copy = non_const_clno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("t_adj copy"), this->ne);
                Kokkos::deep_copy(t_adj_copy, this->t_adj);

                // Debugging
                // prettyPrint1DView(t_adj_copy, "t_adj_copy", 100);
                // prettyPrint1DView(t_adj, "t_adj", 100);

                this->colorGreedyEF(
                  this->xadj, adj_copy, this->t_xadj, t_adj_copy, colors_out, current_vertexList, current_vertexListLength);
            }
            else
            {
                this->colorGreedy(
                  this->xadj, this->adj, this->t_xadj, this->t_adj, colors_out, current_vertexList, current_vertexListLength);
            }

            my_exec_space().fence();

            if(this->_ticToc)
            {
                time = timer.seconds();
                total_time += time;
                std::cout << "\tIteration: " << iter << std::endl
                          << "\t  - Time speculative greedy phase : " << time << std::endl
                          << "\t  - Num Uncolored (greedy-color)  : " << numUncolored << std::endl;

                gc_handle->add_to_overall_coloring_time_phase1(time);

                timer.reset();
            }

            // ------------------------------------------
            // Find conflicts
            // ------------------------------------------
            bool swap_work_arrays = true;      // NOTE: swap_work_arrays can go away in this example -- was only ever
                                               //       set false in the PPS code in the original D1 coloring...

            // NOTE: not using colorset algorithm in this so we don't include colorset data
            numUncolored = this->findConflicts(swap_work_arrays,
                                               this->xadj,
                                               this->adj,
                                               this->t_xadj,
                                               this->t_adj,
                                               colors_out,
                                               current_vertexList,
                                               current_vertexListLength,
                                               next_iteration_recolorList,
                                               next_iteration_recolorListLength);

            my_exec_space().fence();

            if(_ticToc)
            {
                time = timer.seconds();
                total_time += time;
                std::cout << "\t  - Time conflict detection       : " << time << std::endl;
                std::cout << "\t  - Num Uncolored (conflicts)     : " << numUncolored << std::endl;
                gc_handle->add_to_overall_coloring_time_phase2(time);
                timer.reset();
            }

            // Swap the work arrays (for conflictlist)
            if(swap_work_arrays)
            {
                // Swap Work Arrays
                if(iter + 1 < this->_max_num_iterations)
                {
                    nnz_lno_temp_work_view_t temp = current_vertexList;
                    current_vertexList            = next_iteration_recolorList;
                    next_iteration_recolorList    = temp;

                    current_vertexListLength         = numUncolored;
                    next_iteration_recolorListLength = single_dim_index_view_type("recolorListLength");
                }
            }

            // Bail out if we didn't make any progress since we're probably stuck and it's better to just clean up in serial.
            if(numUncolored == numUncoloredPreviousIter)
                break;

        }      // end for iterations && numUncolored > 0...

        // ------------------------------------------
        // clean up in serial (resolveConflictsSerial)
        // ------------------------------------------
        if(numUncolored > 0)
        {
            this->resolveConflictsSerial(this->nv,
                                         this->xadj,
                                         this->adj,
                                         this->t_xadj,
                                         this->t_adj,
                                         colors_out,
                                         current_vertexList,
                                         current_vertexListLength);
        }

        my_exec_space().fence();

        if(_ticToc)
        {
            time = timer.seconds();
            total_time += time;
            std::cout << "\tTime serial conflict resolution   : " << time << std::endl;
            gc_handle->add_to_overall_coloring_time_phase3(time);
        }

        // Save the number of phases and vertex colors to the graph coloring handle
        this->gc_handle->set_vertex_colors(colors_out);
        this->gc_handle->set_num_phases((double)iter);

    }      // color_graph_d2 (end)



    /**
     *  Validate Distance 2 Graph Coloring
     *
     *  @param[in] validation_flags Is an array of 3 booleans to capture if the graph coloring is invalid, and if so what is
     * invalid. validation_flags[0] : True IF the distance-2 coloring is invalid. validation_flags[1] : True IF the coloring is
     * bad because vertices are left uncolored. validation_flags[2] : True IF the coloring is bad because at least one pair of
     * vertices at distance=2 from each other has the same color. validation_flags[3] : True IF a vertex has a color that is
     * greater than number of vertices in the graph. This may not be an INVALID coloring, but can indicate poor quality in
     * coloring.
     *
     *  @return boolean that is TRUE if the Distance-2 coloring is valid. False if otherwise.
     */
    bool verify_coloring(const_lno_row_view_type xadj_,
                         const_lno_nnz_view_t    adj_,
                         const_clno_row_view_t   t_xadj_,
                         const_clno_nnz_view_t   t_adj_,
                         color_view_type         vertex_colors_,
                         bool                    validation_flags[])
    {
        bool output = true;

        validation_flags[ 0 ] = false;      // True if an invalid coloring.
        validation_flags[ 1 ] = false;      // True if uncolored vertices exist.
        validation_flags[ 2 ] = false;      // True if an invalid color exists.

        nnz_lno_type chunkSize_ = this->_chunkSize;
        if(nv < 100 * chunkSize_)
        {
            chunkSize_ = 1;
        }
        const size_t num_chunks = this->nv / chunkSize_ + 1;

        non_const_1d_bool_view_t d_flags("flags", 4);

        // Create host mirror view.
        non_const_1d_bool_view_t::HostMirror h_flags = Kokkos::create_mirror_view(d_flags);

        // Initialize flags on host
        for(int i = 0; i < 4; i++) h_flags(i) = false;

        // Deep copy initialized flags to device
        Kokkos::deep_copy(d_flags, h_flags);

        functorVerifyDistance2Coloring vr(this->nv, xadj_, adj_, t_xadj_, t_adj_, vertex_colors_, d_flags, chunkSize_);
        Kokkos::parallel_for("ValidateD2Coloring", range_policy_type(0, num_chunks), vr);

        // Deep copy flags back to host memory
        Kokkos::deep_copy(h_flags, d_flags);

        validation_flags[ 0 ] = h_flags[ 0 ];
        validation_flags[ 1 ] = h_flags[ 1 ];
        validation_flags[ 2 ] = h_flags[ 2 ];
        validation_flags[ 3 ] = h_flags[ 3 ];

        output = !h_flags[ 0 ];

        return output;
    }      // verify_coloring (end)



    /**
     * Generate a histogram of the distance-2 colors
     *
     * @param[out] histogram A Kokkos::View that will contain the histogram.
     *                       Must be of size num_colors+1
     *
     * @return None
     */
    void compute_color_histogram(nnz_lno_temp_work_view_t& histogram)
    {
        KokkosKernels::Impl::kk_get_histogram<typename HandleType::color_view_type, nnz_lno_temp_work_view_t, my_exec_space>
            (this->nv, this->gc_handle->get_vertex_colors(), histogram);

        my_exec_space().fence();
    }



    /**
     * Print out the histogram of colors in CSV format.
     */
    void print_color_histogram_csv()
    {
        nnz_lno_type num_colors = this->gc_handle->get_num_colors();

        nnz_lno_temp_work_view_t histogram("histogram", num_colors + 1);

        this->compute_color_histogram(histogram);

        auto h_histogram = Kokkos::create_mirror_view(histogram);
        Kokkos::deep_copy(h_histogram, histogram);

        size_t i = 0;
        for(i = 1; i < h_histogram.extent(0) - 1; i++)
        {
            std::cout << h_histogram(i) << ",";
        }
        std::cout << h_histogram(i);
    }



    /**
     * Print out the histogram of colors in a more human friendly format
     * This will not print out all the colors if there are many.
     */
    void print_color_histogram()
    {
        nnz_lno_type             num_colors = this->gc_handle->get_num_colors();
        nnz_lno_temp_work_view_t histogram("histogram", num_colors + 1);

        this->compute_color_histogram(histogram);
        auto histogram_slice = Kokkos::subview(histogram, std::make_pair((size_t)1, histogram.extent(0)));
        std::cout << "Distance-2 Color Histogram (1..N): " << std::endl;
        KokkosKernels::Impl::kk_print_1Dview(histogram_slice);
        std::cout << std::endl;
    }



    /**
     * Calculate the distance-2 degree of all the vertices in the graph.
     *
     * @param[out] degree_d2     A mutable view of size |V|
     * @param[out] degree_d2_max Saves the max distance-2 degree value in.
     *
     * @note This is EXPERIMENTAL development code.
     * @todo Uses HashMapAccumulator and still needs a lot of tweaking to make performant.
     */
    void compute_distance2_degree(non_const_1d_size_type_view_type& degree_d2, size_t& degree_d2_max)
    {
        // Vertex group chunking
        nnz_lno_type v_chunk_size = this->_chunkSize;
        nnz_lno_type v_num_chunks = this->nv / v_chunk_size + 1;

        std::cout << ">>> this->nv     = " << this->nv << std::endl;
        std::cout << ">>> v_chunk_size = " << v_chunk_size << std::endl;
        std::cout << ">>> v_num_chunks = " << v_num_chunks << std::endl;

        // Get the maximum distance-1 degree
        nnz_lno_type max_d1_degree = this->compute_max_distance1_degree();

        // Round up hash_size to next power of two
        nnz_lno_type hash_size = 1;
        while(hash_size < max_d1_degree * max_d1_degree) { hash_size *= 2; }

        // Max Nonzeros in D2 context can be max_degree(G)^2
        nnz_lno_type max_nonzeros = max_d1_degree * max_d1_degree;

        // Determine how much we'd need for UniformMemoryPool
        nnz_lno_type mem_chunk_size = hash_size;      // For hash indices
        mem_chunk_size += hash_size;                  // For hash begins
        mem_chunk_size += max_nonzeros;               // For hash nexts
        mem_chunk_size += max_nonzeros;               // For hash keys
        // nnz_lno_type mem_chunk_size += max_nonzeros;   // For hash_values

        /*
        std::cout << "num_verts      = " << this->nv << std::endl;
        std::cout << "v_chunk_size   = " << v_chunk_size << std::endl;
        std::cout << "v_num_chunks   = " << v_num_chunks << std::endl;
        std::cout << "max_d1_degree  = " << max_d1_degree << std::endl;
        std::cout << "mem_chunk_size = " << mem_chunk_size << std::endl;
        std::cout << "hash_size      = " << hash_size << std::endl;
        */

        // HashmapAccumulator requires a memory Pool
        KokkosKernels::Impl::PoolType my_pool_type = KokkosKernels::Impl::OneThread2OneChunk;
        pool_memory_space_t           m_space(v_num_chunks, mem_chunk_size, -1, my_pool_type);

        /*
        std::cout << std::endl << "Memory Pool:" << std::endl;
        m_space.print_memory_pool(true);
        */

        functorCalculateD2Degree calculateD2Degree(
          this->nv, this->xadj, this->adj, this->t_xadj, this->t_adj, v_chunk_size, degree_d2, m_space, hash_size, max_nonzeros);
        Kokkos::parallel_for("Compute Degree D2", range_policy_type(0, v_num_chunks), calculateD2Degree);

        // Compute maximum d2 degree
        size_t _degree_d2_max = 0;
        Kokkos::parallel_reduce(
          "Max D2 Degree",
          this->nv,
          KOKKOS_LAMBDA(const size_t& i, size_t & lmax) { lmax = degree_d2(i) > lmax ? degree_d2(i) : lmax; },
          Kokkos::Max<size_t>(_degree_d2_max));
        degree_d2_max = _degree_d2_max;

        // Tell the handle that coloring has been called.
        this->gc_handle->set_coloring_called();
    }



  private:
    // -----------------------------------------------------------------
    //
    // GraphColorDistance2::colorGreedy()
    //
    // -----------------------------------------------------------------
    void colorGreedy(const_lno_row_view_type  xadj_,
                     const_lno_nnz_view_t     adj_,
                     const_clno_row_view_t    t_xadj_,
                     const_clno_nnz_view_t    t_adj_,
                     color_view_type          vertex_colors_,
                     nnz_lno_temp_work_view_t current_vertexList_,
                     nnz_lno_type             current_vertexListLength_)
    {
        nnz_lno_type chunkSize_ = this->_chunkSize;

        if(current_vertexListLength_ < 100 * chunkSize_)
        {
            chunkSize_ = 1;
        }

        // Pick the right coloring algorithm to use based on which algorithm we're using
        switch(this->gc_handle->get_coloring_algo_type())
        {
            // Single level parallelism on chunks
            // 1. [P] loop over vertices
            // 2. [S] loop over color offset blocks
            // 3. [S] loop over vertex neighbors
            // 4. [S] loop over vertex neighbors of neighbors
            case COLORING_D2:
            case COLORING_D2_VB:
            {
                functorGreedyColorVB gc(
                  this->nv, xadj_, adj_, t_xadj_, t_adj_, vertex_colors_, current_vertexList_, current_vertexListLength_);
                Kokkos::parallel_for("LoopOverChunks", range_policy_type(0, this->nv), gc);
            }
            break;

            // One level Perallelism, BIT Array for coloring
            // 1. [P] loop over vertices
            // 2. [S] loop over color offset blocks
            // 3. [S] loop over vertex neighbors
            // 4. [S] loop over vertex neighbors of neighbors
            case COLORING_D2_VB_BIT:
            {
                functorGreedyColorVB_BIT gc(
                  this->nv, xadj_, adj_, t_xadj_, t_adj_, vertex_colors_, current_vertexList_, current_vertexListLength_);
                Kokkos::parallel_for("LoopOverChunks", range_policy_type(0, this->nv), gc);
            }
            break;

            default:
                throw std::invalid_argument("Unknown Distance-2 Algorithm Type or invalid for non Edge Filtering mode.");
        }

    }      // colorGreedy (end)



    // -----------------------------------------------------------------
    //
    // GraphColorDistance2::colorGreedyEF()
    //
    // -----------------------------------------------------------------
    void colorGreedyEF(const_lno_row_view_type   xadj_,
                       nnz_lno_temp_work_view_t  adj_copy_,
                       const_clno_row_view_t     t_xadj_,
                       non_const_clno_nnz_view_t t_adj_copy_,
                       color_view_type           vertex_colors_,
                       nnz_lno_temp_work_view_t  current_vertexList_,
                       nnz_lno_type              current_vertexListLength_)
    {
        // Pick the right coloring algorithm to use based on which algorithm we're using
        switch(this->gc_handle->get_coloring_algo_type())
        {
            // One level Perallelism, BIT Array for coloring + edge filtering
            // 1. [P] loop over vertices
            // 2. [S] loop over color offset blocks
            // 3. [S] loop over vertex neighbors
            // 4. [S] loop over vertex neighbors of neighbors
            case COLORING_D2_VB_BIT_EF:
            {
                functorGreedyColorVB_BIT_EF gc(this->nv,
                                               xadj_,
                                               adj_copy_,
                                               t_xadj_,
                                               t_adj_copy_,
                                               vertex_colors_,
                                               current_vertexList_,
                                               current_vertexListLength_);
                Kokkos::parallel_for("LoopOverChunks", range_policy_type(0, this->nv), gc);
                // prettyPrint1DView(vertex_colors_, "COLORS_GC_VB_BIT",500);
            }
            break;

            default:
                throw std::invalid_argument("Unknown Distance-2 Algorithm Type or algorithm does not use Edge Filtering.");
        }
    }      // colorGreedyEF (end)



    // -----------------------------------------------------------------
    //
    // GraphColorDistance2::findConflicts()
    //
    // -----------------------------------------------------------------
    template<typename adj_view_t>
    nnz_lno_type findConflicts(bool&                      swap_work_arrays,
                               const_lno_row_view_type    xadj_,
                               adj_view_t                 adj_,
                               const_clno_row_view_t      t_xadj_,
                               const_clno_nnz_view_t      t_adj_,
                               color_view_type            vertex_colors_,
                               nnz_lno_temp_work_view_t   current_vertexList_,
                               nnz_lno_type               current_vertexListLength_,
                               nnz_lno_temp_work_view_t   next_iteration_recolorList_,
                               single_dim_index_view_type next_iteration_recolorListLength_)
    {
        swap_work_arrays                 = true;
        nnz_lno_type output_numUncolored = 0;

        if(0 == this->_use_color_set)
        {
            functorFindConflicts_Atomic<adj_view_t> conf(this->nv,
                                                         xadj_,
                                                         adj_,
                                                         t_xadj_,
                                                         t_adj_,
                                                         vertex_colors_,
                                                         current_vertexList_,
                                                         next_iteration_recolorList_,
                                                         next_iteration_recolorListLength_);
            Kokkos::parallel_reduce("FindConflicts", range_policy_type(0, current_vertexListLength_), conf, output_numUncolored);
        }

        return output_numUncolored;
    }      // findConflicts (end)



    // -----------------------------------------------------------------
    //
    // GraphColorDistance2::resolveConflictsSerial()
    //
    // -----------------------------------------------------------------
    template<typename adj_view_t>
    void resolveConflictsSerial(nnz_lno_type             _nv,
                                const_lno_row_view_type  xadj_,
                                adj_view_t               adj_,
                                const_clno_row_view_t    t_xadj_,
                                const_clno_nnz_view_t    t_adj_,
                                color_view_type          vertex_colors_,
                                nnz_lno_temp_work_view_t current_vertexList_,
                                size_type                current_vertexListLength_)
    {
        color_type*  forbidden = new color_type[ _nv ];
        nnz_lno_type vid       = 0;
        nnz_lno_type end       = _nv;

        typename nnz_lno_temp_work_view_t::HostMirror h_recolor_list;

        end            = current_vertexListLength_;
        h_recolor_list = Kokkos::create_mirror_view(current_vertexList_);
        Kokkos::deep_copy(h_recolor_list, current_vertexList_);

        auto h_colors = Kokkos::create_mirror_view(vertex_colors_);

        auto h_idx = Kokkos::create_mirror_view(xadj_);
        auto h_adj = Kokkos::create_mirror_view(adj_);

        auto h_t_idx = Kokkos::create_mirror_view(t_xadj_);
        auto h_t_adj = Kokkos::create_mirror_view(t_adj_);

        Kokkos::deep_copy(h_colors, vertex_colors_);

        Kokkos::deep_copy(h_idx, xadj_);
        Kokkos::deep_copy(h_adj, adj_);

        Kokkos::deep_copy(h_t_idx, t_xadj_);
        Kokkos::deep_copy(h_t_adj, t_adj_);

        for(nnz_lno_type k = 0; k < end; k++)
        {
            vid = h_recolor_list(k);

            if(h_colors(vid) > 0)
                continue;

            // loop over distance-1 neighbors of vid
            for(size_type vid_d1_adj = h_idx(vid); vid_d1_adj < h_idx(vid + 1); vid_d1_adj++)
            {
                size_type vid_d1 = h_adj(vid_d1_adj);

                // loop over neighbors of vid_d1 (distance-2 from vid)
                for(size_type vid_d2_adj = h_t_idx(vid_d1); vid_d2_adj < h_t_idx(vid_d1 + 1); vid_d2_adj++)
                {
                    nnz_lno_type vid_d2 = h_t_adj(vid_d2_adj);

                    // skip over loops vid -- x -- vid
                    if(vid_d2 == vid)
                        continue;

                    forbidden[ h_colors(vid_d2) ] = vid;
                }
            }

            // color vertex vid with smallest available color
            int c = 1;
            while(forbidden[ c ] == vid) c++;

            h_colors(vid) = c;
        }
        Kokkos::deep_copy(vertex_colors_, h_colors);
        delete[] forbidden;
    }      // resolveConflictsSerial (end)



    // ------------------------------------------------------
    // Functions: Helpers
    // ------------------------------------------------------


  public:
    /**
     * Compute the maximum distance-1 degree.
     *
     * @return Maximum distance1 degree in the graph
     */
    nnz_lno_type compute_max_distance1_degree() const
    {
        nnz_lno_type max_degree = 0;
        Kokkos::parallel_reduce("Max D1 Degree",
                                this->nv,
                                KOKKOS_LAMBDA(const nnz_lno_type& vid, nnz_lno_type & lmax) {
                                    const nnz_lno_type degree = this->xadj(vid + 1) - this->xadj(vid);
                                    lmax                      = degree > lmax ? degree : lmax;
                                },
                                Kokkos::Max<nnz_lno_type>(max_degree));
        return max_degree;
    }



    // pretty-print a 1D View with label
    template<typename kokkos_view_t>
    void prettyPrint1DView(kokkos_view_t& view, const char* label, const size_t max_entries = 500) const
    {
        int max_per_line = 20;
        int line_count   = 1;
        std::cout << label << " = [ \n\t";
        for(size_t i = 0; i < view.extent(0); i++)
        {
            std::cout << std::setw(5) << view(i) << " ";
            if(line_count >= max_per_line)
            {
                std::cout << std::endl << "\t";
                line_count = 0;
            }
            line_count++;
            if(i >= max_entries - 1)
            {
                std::cout << "<snip>";
                break;
            }
        }
        if(line_count > 1)
            std::cout << std::endl;
        std::cout << "\t ]" << std::endl;
    }      // prettyPrint1DView (end)



    // ------------------------------------------------------
    // Functors: Distance-2 Graph Coloring
    // ------------------------------------------------------

    /**
     * Functor to init a list sequentialy, that is list[i] = i
     */
    template<typename view_type>
    struct functorInitList
    {
        view_type _vertexList;
        functorInitList(view_type vertexList) : _vertexList(vertexList) {}

        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type i) const
        {
            // Natural order
            _vertexList(i) = i;
        }
    };      // struct functorInitList (end)



    /**
     * Functor for VB algorithm speculative coloring without edge filtering.
     * Single level parallelism
     */
    struct functorGreedyColorVB
    {
        nnz_lno_type             nv;                     // num vertices
        const_lno_row_view_type  _idx;                   // vertex degree list
        const_lno_nnz_view_t     _adj;                   // vertex adjacency list
        const_clno_row_view_t    _t_idx;                 // transpose vertex degree list
        const_clno_nnz_view_t    _t_adj;                 // transpose vertex adjacency list
        color_view_type          _colors;                // vertex colors
        nnz_lno_temp_work_view_t _vertexList;            //
        nnz_lno_type             _vertexListLength;      //
        nnz_lno_type             _chunkSize;             //

        functorGreedyColorVB(nnz_lno_type             nv_,
                             const_lno_row_view_type  xadj_,
                             const_lno_nnz_view_t     adj_,
                             const_clno_row_view_t    t_xadj_,
                             const_clno_nnz_view_t    t_adj_,
                             color_view_type          colors,
                             nnz_lno_temp_work_view_t vertexList,
                             nnz_lno_type             vertexListLength)
            : nv(nv_)
            , _idx(xadj_)
            , _adj(adj_)
            , _t_idx(t_xadj_)
            , _t_adj(t_adj_)
            , _colors(colors)
            , _vertexList(vertexList)
            , _vertexListLength(vertexListLength)
        {
        }


        // Color vertex i with smallest available color.
        //
        // Each thread colors a chunk of vertices to prevent all vertices getting the same color.
        //
        // This version uses a bool array of size FORBIDDEN_SIZE.
        //
        // param: ii = vertex id
        //
        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type vid) const
        {
            // If vertex is not already colored...
            if(_colors(vid) <= 0)
            {
                bool foundColor = false;      // Have we found a valid color?

                const size_type vid_adj_begin = _idx(vid);
                const size_type vid_adj_end   = _idx(vid + 1);

                // Use forbidden array to find available color.
                // - should be small enough to fit into fast memory (use Kokkos memoryspace?)
                bool forbidden[ VB_D2_COLORING_FORBIDDEN_SIZE ];      // Forbidden Colors

                // Do multiple passes if the forbidden array is too small.
                // * The Distance-1 code used the knowledge of the degree of the vertex to cap the number of iterations
                //   but in distance-2 we'd need the total vertices at distance-2 which we don't easily have aprioi.
                //   This could be as big as all the vertices in the graph if diameter(G)=2...
                for(color_type offset = 0; !foundColor && offset < nv; offset += VB_D2_COLORING_FORBIDDEN_SIZE)
                {
                    // initialize
                    for(int i = 0; i < VB_D2_COLORING_FORBIDDEN_SIZE; i++) { forbidden[ i ] = false; }

                    // Colors start at 1.  0 is special in that it means a vertex is uncolored.
                    // For the range 0..63 we mark forbidden[0] as true to take color 0 out of
                    // consideration.
                    if(0 == offset)
                    {
                        forbidden[ 0 ] = true;
                    }

                    // Check neighbors, fill forbidden array.
                    for(size_type vid_adj = vid_adj_begin; vid_adj < vid_adj_end; vid_adj++)
                    {
                        const nnz_lno_type vid_d1           = _adj(vid_adj);
                        const size_type    vid_d1_adj_begin = _t_idx(vid_d1);
                        const size_type    vid_d1_adj_end   = _t_idx(vid_d1 + 1);

                        for(size_type vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; vid_d1_adj++)
                        {
                            const nnz_lno_type vid_d2 = _t_adj(vid_d1_adj);

                            // Skip distance-2-self-loops
                            if(vid_d2 != vid && vid_d2 < nv)
                            {
                                const color_type c = _colors(vid_d2);

                                if((c >= offset) && (c - offset < VB_D2_COLORING_FORBIDDEN_SIZE))
                                {
                                    forbidden[ c - offset ] = true;
                                }
                            }
                        }      // for vid_d1_adj...
                    }          // for vid_adj...

                    // color vertex i with smallest available color (firstFit)
                    for(int c = 0; c < VB_D2_COLORING_FORBIDDEN_SIZE; c++)
                    {
                        if(!forbidden[ c ])
                        {
                            _colors(vid) = offset + c;
                            foundColor   = true;
                            break;
                        }
                    }      // for c...
                }          // for !foundColor...
            }              // if _colors(vid) <= 0...
        }                  // operator() (end)
    };                     // struct functorGreedyColorVB (end)



    /**
     * Functor for VB_BIT algorithm coloring without edge filtering.
     * Single level parallelism
     */
    struct functorGreedyColorVB_BIT
    {
        nnz_lno_type             nv;                     // num vertices
        const_lno_row_view_type  _idx;                   // vertex degree list
        const_lno_nnz_view_t     _adj;                   // vertex adjacency list
        const_clno_row_view_t    _t_idx;                 // transpose vertex degree list
        const_clno_nnz_view_t    _t_adj;                 // transpose vertex adjacency list
        color_view_type          _colors;                // vertex colors
        nnz_lno_temp_work_view_t _vertexList;            //
        nnz_lno_type             _vertexListLength;      //

        functorGreedyColorVB_BIT(nnz_lno_type             nv_,
                                 const_lno_row_view_type  xadj_,
                                 const_lno_nnz_view_t     adj_,
                                 const_clno_row_view_t    t_xadj_,
                                 const_clno_nnz_view_t    t_adj_,
                                 color_view_type          colors,
                                 nnz_lno_temp_work_view_t vertexList,
                                 nnz_lno_type             vertexListLength)
            : nv(nv_)
            , _idx(xadj_)
            , _adj(adj_)
            , _t_idx(t_xadj_)
            , _t_adj(t_adj_)
            , _colors(colors)
            , _vertexList(vertexList)
            , _vertexListLength(vertexListLength)
        {
        }


        // Color vertex i with smallest available color.
        //
        // Each thread colors a chunk of vertices to prevent all vertices getting the same color.
        //
        // This version uses a bool array of size FORBIDDEN_SIZE.
        //
        // param: ii = vertex id
        //
        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type vid) const
        {
            // If vertex is not colored yet...
            if(_colors(vid) == 0)
            {
                const size_type vid_adj_begin = _idx(vid);
                const size_type vid_adj_end   = _idx(vid + 1);

                bool foundColor = false;
                for(color_type offset = 0; !foundColor && offset <= (nv + VBBIT_D2_COLORING_FORBIDDEN_SIZE);
                    offset += VBBIT_D2_COLORING_FORBIDDEN_SIZE)
                {
                    // Forbidden colors
                    // - single long int for forbidden colors
                    bit_64_forbidden_type forbidden = 0;

                    // If all available colors for this range are unavailable we can break out of the nested loops
                    bool break_out = false;

                    // Loop over distance-1 neighbors of vid
                    for(size_type vid_adj = vid_adj_begin; !break_out && vid_adj < vid_adj_end; ++vid_adj)
                    {
                        const nnz_lno_type vid_d1           = _adj(vid_adj);
                        const size_type    vid_d1_adj_begin = _t_idx(vid_d1);
                        const size_type    vid_d1_adj_end   = _t_idx(vid_d1 + 1);

                        // Loop over distance-2 neighbors of vid
                        for(size_type vid_d1_adj = vid_d1_adj_begin; !break_out && vid_d1_adj < vid_d1_adj_end; ++vid_d1_adj)
                        {
                            const nnz_lno_type vid_d2 = _t_adj(vid_d1_adj);

                            // Ignore Distance-2 Self Loops
                            if(vid_d2 != vid && vid_d2 < nv)
                            {
                                const color_type color        = _colors(vid_d2);
                                const color_type color_offset = color - offset;

                                // if color is within the current range, or if its color is in a previously traversed
                                // range
                                if(color && color_offset <= VBBIT_D2_COLORING_FORBIDDEN_SIZE)
                                {
                                    // if it is in the current range, then add the color to the banned colors
                                    if(color > offset)
                                    {
                                        // convert color to bit representation
                                        bit_64_forbidden_type ban_color_bit = 1;

                                        ban_color_bit = ban_color_bit << (color_offset - 1);

                                        // add it to forbidden colors
                                        forbidden = forbidden | ban_color_bit;

                                        // if there are no available colors in this range then exit early,
                                        // no need to traverse the rest.
                                        if(0 == ~forbidden)
                                        {
                                            break_out = true;
                                        }
                                    }      // if color > offset ...
                                }          // if color && color_offset ...
                            }              // if vid_d2 ...
                        }                  // for vid_d1_adj ...
                    }                      // for vid_adj ...
                    forbidden = ~(forbidden);

                    // check if an available color exists.
                    if(forbidden)
                    {
                        color_type val = 1;

                        // if there is an available color, choose the first color, using 2s complement.
                        bit_64_forbidden_type new_color = forbidden & (-forbidden);

                        // convert it back to decimal color.
                        while((new_color & 1) == 0)
                        {
                            ++val;
                            new_color = new_color >> 1;
                        }
                        _colors(vid) = val + offset;
                        foundColor   = true;
                        break;
                    }
                }      // for !foundColor
            }          // if _colors(vid)==0
        }              // operator() (end)
    };                 // struct functorGreedyColorVB_BIT (end)



    /**
     * Functor for VB_BIT_EF algorithm coloring without edge filtering.
     * Single level parallelism
     */
    struct functorGreedyColorVB_BIT_EF
    {
        nnz_lno_type              _nv;                    // num vertices
        const_lno_row_view_type   _idx;                   // vertex degree list
        nnz_lno_temp_work_view_t  _adj;                   // vertex adjacency list  (mutable)
        const_clno_row_view_t     _t_idx;                 // transpose vertex degree list
        non_const_clno_nnz_view_t _t_adj;                 // transpose vertex adjacency list (mutable)
        color_view_type           _colors;                // vertex colors
        nnz_lno_temp_work_view_t  _vertexList;            //
        nnz_lno_type              _vertexListLength;      //

        functorGreedyColorVB_BIT_EF(nnz_lno_type              nv_,
                                    const_lno_row_view_type   xadj_,
                                    nnz_lno_temp_work_view_t  adj_,
                                    const_clno_row_view_t     t_xadj_,
                                    non_const_clno_nnz_view_t t_adj_,
                                    color_view_type           colors,
                                    nnz_lno_temp_work_view_t  vertexList,
                                    nnz_lno_type              vertexListLength)
            : _nv(nv_)
            , _idx(xadj_)
            , _adj(adj_)
            , _t_idx(t_xadj_)
            , _t_adj(t_adj_)
            , _colors(colors)
            , _vertexList(vertexList)
            , _vertexListLength(vertexListLength)
        {
        }


        // Color vertex i with smallest available color.
        //
        // Each thread colors a chunk of vertices to prevent all vertices getting the same color.
        //
        // This version uses a bool array of size FORBIDDEN_SIZE.
        //
        // param: ii = vertex id
        //
        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type vid) const
        {
            // If vertex is not colored yet..
            if(_colors(vid) == 0)
            {
                size_type vid_adj_begin = _idx(vid);
                size_type vid_adj_end   = _idx(vid + 1);

                bool foundColor = false;
                for(color_type offset = 0; !foundColor && offset <= (_nv + VBBIT_D2_COLORING_FORBIDDEN_SIZE);
                    offset += VBBIT_D2_COLORING_FORBIDDEN_SIZE)
                {
                    // Forbidden colors
                    // - single long int for forbidden colors
                    bit_64_forbidden_type forbidden = 0;

                    // If all available colors for this range are unavailable we can break out of the nested loops
                    bool offset_colors_full = false;

                    // Loop over distance-1 neighbors of vid
                    for(size_type vid_adj = vid_adj_begin; !offset_colors_full && vid_adj < vid_adj_end; ++vid_adj)
                    {
                        const nnz_lno_type vid_d1 = _adj(vid_adj);

                        size_type       vid_d1_adj_begin            = _t_idx(vid_d1);
                        const size_type vid_d1_adj_end              = _t_idx(vid_d1 + 1);
                        const size_type degree_vid_d1               = vid_d1_adj_end - vid_d1_adj_begin;
                        size_type       num_vid_d2_colored_in_range = 0;

                        // Store the maximum color value found in the vertices adjacent to vid_d1
                        color_type max_color_adj_to_d1 = 0;

                        // Loop over distance-2 neighbors of vid
                        for(size_type vid_d1_adj = vid_d1_adj_begin; !offset_colors_full && vid_d1_adj < vid_d1_adj_end;
                            ++vid_d1_adj)
                        {
                            const nnz_lno_type vid_d2 = _t_adj(vid_d1_adj);

                            // Ignore Distance-2 Self Loops
                            if(vid_d2 != vid && vid_d2 < _nv)
                            {
                                color_type color = _colors(vid_d2);
                                color_type color_offset =
                                  color - offset;      // color_offset < 0 means color is from a previous offset.

                                // Update maximum color adjacent to vid_d1 found so far.
                                max_color_adj_to_d1 = color > max_color_adj_to_d1 ? color : max_color_adj_to_d1;

                                // if color is within the current range, or if its color is in a previously traversed
                                // range
                                if(color && color_offset <= VBBIT_D2_COLORING_FORBIDDEN_SIZE)
                                {
                                    num_vid_d2_colored_in_range++;

                                    // if it is in the current range, then add the color to the banned colors
                                    if(color > offset)
                                    {
                                        // convert color to bit representation
                                        bit_64_forbidden_type ban_color_bit = 1;

                                        ban_color_bit = ban_color_bit << (color_offset - 1);

                                        // add it to forbidden colors
                                        forbidden = forbidden | ban_color_bit;

                                        // if there are no available colors in this range then exit early,
                                        // no need to traverse the rest b/c they contribute no new information
                                        // at this offset.
                                        if(0 == ~forbidden)
                                        {
                                            offset_colors_full = true;
                                            // Note: with edge-filtering, this can short-circuit the loop over all
                                            //       neighbors of VID and will reduce the number of filtered edges.
                                        }
                                    }      // if color > offset
                                }          // if color && color_offset
                            }              // if vid_d2 != vid ...
                            else
                            {
                                // If there's a self-loop then we should increment our 'colored in range' so we don't
                                // block filtering since we know there must be a (v2,v1) edge
                                num_vid_d2_colored_in_range++;
                            }
                        }      // for vid_d1_adj ...

                        // Edge filtering on the neighbors of vid.  We can only do this if ALL neighbors of vid_d1
                        // have been visited and if all were colored in current offset range or lower.
                        if(degree_vid_d1 == num_vid_d2_colored_in_range)
                        {
                            if(vid_adj_begin > vid_adj)
                            {
                                _adj(vid_adj)       = _adj(vid_adj_begin);
                                _adj(vid_adj_begin) = vid_d1;
                            }
                            vid_adj_begin++;
                        }

                    }      // for vid_adj
                    forbidden = ~(forbidden);

                    // check if an available color exists.
                    if(forbidden)
                    {
                        // if there is an available color, choose the first color, using 2s complement.
                        bit_64_forbidden_type new_color = forbidden & (-forbidden);
                        color_type         val       = 1;

                        // convert it back to decimal color.
                        while((new_color & 1) == 0)
                        {
                            ++val;
                            new_color = new_color >> 1;
                        }
                        _colors(vid) = val + offset;
                        foundColor   = true;
                        break;
                    }
                }      // for offset=0...
            }          // if _colors(vid)==0
        }              // operator() (end)
    };                 // struct functorGreedyColorVB_BIT_EF (end)



    template<typename adj_view_t>
    struct functorFindConflicts_Atomic
    {
        nnz_lno_type               nv;      // num verts
        const_lno_row_view_type    _idx;
        adj_view_t                 _adj;
        const_clno_row_view_t      _t_idx;
        const_clno_nnz_view_t      _t_adj;
        color_view_type            _colors;
        nnz_lno_temp_work_view_t   _vertexList;
        nnz_lno_temp_work_view_t   _recolorList;
        single_dim_index_view_type _recolorListLength;


        functorFindConflicts_Atomic(nnz_lno_type               nv_,
                                    const_lno_row_view_type    xadj_,
                                    adj_view_t                 adj_,
                                    const_clno_row_view_t      t_xadj_,
                                    const_clno_nnz_view_t      t_adj_,
                                    color_view_type            colors,
                                    nnz_lno_temp_work_view_t   vertexList,
                                    nnz_lno_temp_work_view_t   recolorList,
                                    single_dim_index_view_type recolorListLength)
            : nv(nv_)
            , _idx(xadj_)
            , _adj(adj_)
            , _t_idx(t_xadj_)
            , _t_adj(t_adj_)
            , _colors(colors)
            , _vertexList(vertexList)
            , _recolorList(recolorList)
            , _recolorListLength(recolorListLength)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type vid_, nnz_lno_type& numConflicts) const
        {
            using atomic_incr_type = typename std::remove_reference<decltype(_recolorListLength())>::type;

            const nnz_lno_type vid      = _vertexList(vid_);
            const color_type   my_color = _colors(vid);

            size_type       vid_d1_adj     = _idx(vid);
            const size_type vid_d1_adj_end = _idx(vid + 1);

            bool break_out = false;

            for(; !break_out && vid_d1_adj < vid_d1_adj_end; vid_d1_adj++)
            {
                nnz_lno_type vid_d1 = _adj(vid_d1_adj);

                for(size_type vid_d2_adj = _t_idx(vid_d1); !break_out && vid_d2_adj < _t_idx(vid_d1 + 1); vid_d2_adj++)
                {
                    const nnz_lno_type vid_d2 = _t_adj(vid_d2_adj);

                    if(vid != vid_d2 && vid_d2 < nv)
                    {
                        if(_colors(vid_d2) == my_color)
                        {
                            _colors(vid) = 0;      // uncolor vertex

                            // Atomically add vertex to recolorList
                            const nnz_lno_type k = Kokkos::atomic_fetch_add(&_recolorListLength(), atomic_incr_type(1));
                            _recolorList(k)      = vid;
                            numConflicts += 1;
                            break_out = true;
                            // break;      // Can exit if vertex gets marked as a conflict.
                        }
                    }      // if vid != vid_d2 ...
                }          // for vid_d2_adj ...
            }              // for vid_d1_adj ...
        }                  // operator() (end)
    };                     // struct functorFindConflicts_Atomic (end)



    /**
     * functorVerifyDistance2Coloring
     *  - Validate correctness of the distance-2 coloring
     */
    struct functorVerifyDistance2Coloring
    {
        nnz_lno_type             nv;           // num vertices
        const_lno_row_view_type  _idx;         // vertex degree list
        const_lno_nnz_view_t     _adj;         // vertex adjacency list
        const_clno_row_view_t    _t_idx;       // transpose vertex degree list
        const_clno_nnz_view_t    _t_adj;       // transpose vertex adjacency list
        color_view_type          _colors;      // vertex colors
        non_const_1d_bool_view_t _flags;       // [0]: valid or not?  [1]=true means something is uncolored. [2]=true means
                                               // something has distance-2 neighbor of same color
        nnz_lno_type _chunkSize;               //

        functorVerifyDistance2Coloring(nnz_lno_type             nv_,
                                       const_lno_row_view_type  xadj_,
                                       const_lno_nnz_view_t     adj_,
                                       const_clno_row_view_t    t_xadj_,
                                       const_clno_nnz_view_t    t_adj_,
                                       color_view_type          colors,
                                       non_const_1d_bool_view_t flags,
                                       nnz_lno_type             chunkSize)
            : nv(nv_)
            , _idx(xadj_)
            , _adj(adj_)
            , _t_idx(t_xadj_)
            , _t_adj(t_adj_)
            , _colors(colors)
            , _flags(flags)
            , _chunkSize(chunkSize)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type chunk_id) const
        {
            bool has_uncolored_vertex            = false;
            bool has_invalid_color               = false;
            bool has_color_bigger_than_num_verts = false;

            // Loop over vertices in chunks
            for(nnz_lno_type chunk_idx = 0; chunk_idx < _chunkSize; chunk_idx++)
            {
                bool         break_out = false;
                nnz_lno_type vid       = chunk_id * _chunkSize + chunk_idx;

                if(vid < nv)
                {
                    const color_type color_vid        = _colors(vid);
                    const size_type  vid_d1_adj_begin = _idx(vid);
                    const size_type  vid_d1_adj_end   = _idx(vid + 1);

                    if(color_vid == 0)
                    {
                        has_uncolored_vertex = true;
                    }

                    if(color_vid > nv)
                    {
                        has_color_bigger_than_num_verts = true;
                    }

                    // Loop over neighbors of vid (distance-1 from vid)
                    for(size_type vid_d1_adj = vid_d1_adj_begin; !break_out && vid_d1_adj < vid_d1_adj_end; ++vid_d1_adj)
                    {
                        const nnz_lno_type vid_d1           = _adj(vid_d1_adj);
                        const size_type    vid_d2_adj_begin = _t_idx(vid_d1);
                        const size_type    vid_d2_adj_end   = _t_idx(vid_d1 + 1);

                        // Loop over neighbors of vid_d1 (distance-2 from vid)
                        for(size_type vid_d2_adj = vid_d2_adj_begin; !break_out && vid_d2_adj < vid_d2_adj_end; ++vid_d2_adj)
                        {
                            const nnz_lno_type vid_d2 = _t_adj(vid_d2_adj);

                            // Ignore Distance-2 self loops
                            if(vid != vid_d2)
                            {
                                const color_type color_d2 = _colors(vid_d2);
                                // If distance-2 neighbor of vid has the same color as vid, then the coloring is invalid
                                if(color_vid == color_d2)
                                {
                                    has_invalid_color = true;
                                    break_out         = true;
                                }
                            }
                        }
                    }
                }
            }
            if(has_uncolored_vertex || has_invalid_color)
                Kokkos::atomic_fetch_or(&_flags[ 0 ], has_uncolored_vertex || has_invalid_color);
            if(has_uncolored_vertex)
                Kokkos::atomic_fetch_or(&_flags[ 1 ], has_uncolored_vertex);
            if(has_invalid_color)
                Kokkos::atomic_fetch_or(&_flags[ 2 ], has_invalid_color);
            if(has_color_bigger_than_num_verts)
                Kokkos::atomic_fetch_or(&_flags[ 3 ], has_color_bigger_than_num_verts);
        }      // operator()
    };         // struct functorVerifyDistance2Coloring (end)




    /**
     * functorCalculateD2Degree
     */
    struct functorCalculateD2Degree
    {
        nnz_lno_type                     _num_verts;      // num vertices
        const_lno_row_view_type          _idx;            // vertex degree list
        const_lno_nnz_view_t             _adj;            // vertex adjacency list
        const_clno_row_view_t            _t_idx;          // transpose vertex degree list
        const_clno_nnz_view_t            _t_adj;          // transpose vertex adjacency list
        nnz_lno_type                     _chunk_size;
        non_const_1d_size_type_view_type _degree_d2;      // Distance-2 Degree (assumes all are initialized to 0)

        // EXPERIMENTAL BEGIN (for HashMapAccumulator)
        pool_memory_space_t _m_space;
        const nnz_lno_type  _hash_size;
        const nnz_lno_type  _max_nonzeros;

        Kokkos::Experimental::UniqueToken<my_exec_space, Kokkos::Experimental::UniqueTokenScope::Global> tokens;
        // EXPERIMENTAL END

        functorCalculateD2Degree(nnz_lno_type                      num_verts,
                                 const_lno_row_view_type           xadj,
                                 const_lno_nnz_view_t              adj,
                                 const_clno_row_view_t             t_xadj,
                                 const_clno_nnz_view_t             t_adj,
                                 nnz_lno_type                      chunk_size,
                                 non_const_1d_size_type_view_type& degree_d2,
                                 pool_memory_space_t               m_space,
                                 const nnz_lno_type                hash_size,
                                 const nnz_lno_type                max_nonzeros)
            : _num_verts(num_verts)
            , _idx(xadj)
            , _adj(adj)
            , _t_idx(t_xadj)
            , _t_adj(t_adj)
            , _chunk_size(chunk_size)
            , _degree_d2(degree_d2)
            , _m_space(m_space)
            , _hash_size(hash_size)
            , _max_nonzeros(max_nonzeros)
        {
        }


        KOKKOS_INLINE_FUNCTION
        void operator()(const nnz_lno_type chunk_id) const
        {
            // BEGIN EXPERIMENTAL
            // Get a unique token since we aren't using TeamPolicy (for UniformMemoryPool, that the HashmapAccumulator requires)
            auto tid = tokens.acquire();

            volatile nnz_lno_type* tmp = NULL;
            while(NULL == tmp) { tmp = (volatile nnz_lno_type*)(_m_space.allocate_chunk(tid)); }

            KokkosKernels::Experimental::HashmapAccumulator<nnz_lno_type, nnz_lno_type, nnz_lno_type> hash_map;

            // What is globally_used_hash_indices?
            nnz_lno_type* globally_used_hash_indices = (nnz_lno_type*)tmp;
            tmp += _hash_size;

            // set up hash_begins
            hash_map.hash_begins = (nnz_lno_type*)(tmp);
            tmp += _hash_size;

            // set up hash_nexts
            hash_map.hash_nexts = (nnz_lno_type*)(tmp);
            tmp += _max_nonzeros;

            // set up hash_keys
            hash_map.keys = (nnz_lno_type*)(tmp);

            // Set key and value sizes
            hash_map.hash_key_size  = _max_nonzeros;      // Max # of nonzeros that can be added into the hash table (unused?)
            hash_map.max_value_size = _max_nonzeros;      // Max # of nonzeros that can be added into the hash table

            nnz_lno_type pow2_hash_func = _hash_size - 1;
            // END EXPERIMENTAL

            for(nnz_lno_type ichunk = 0; ichunk < _chunk_size; ichunk++)
            {
                nnz_lno_type vid = chunk_id * _chunk_size + ichunk;

                if(vid >= _num_verts)
                    continue;

                nnz_lno_type used_hash_size           = 0;
                nnz_lno_type globally_used_hash_count = 0;

                const size_type vid_d1_adj_begin = _idx(vid);
                const size_type vid_d1_adj_end   = _idx(vid + 1);

                // Loop over neighbors of vid (distance-1 from vid)
                for(size_type vid_d1_adj = vid_d1_adj_begin; vid_d1_adj < vid_d1_adj_end; ++vid_d1_adj)
                {
                    const nnz_lno_type vid_d1           = _adj(vid_d1_adj);
                    const size_type    vid_d2_adj_begin = _t_idx(vid_d1);
                    const size_type    vid_d2_adj_end   = _t_idx(vid_d1 + 1);

                    // EXPERIMENTAL BEGIN
                    for(size_type vid_d2_adj = vid_d2_adj_begin; vid_d2_adj < vid_d2_adj_end; ++vid_d2_adj)
                    {

                        const nnz_lno_type vid_d2 = _t_adj(vid_d2_adj);

                        nnz_lno_type hash = vid_d2 & pow2_hash_func;

                        int r = hash_map.sequential_insert_into_hash_TrackHashes(hash,
                                                                                 vid_d2,
                                                                                 &used_hash_size,
                                                                                 hash_map.max_value_size,
                                                                                 &globally_used_hash_count,
                                                                                 globally_used_hash_indices);

                        // TODO: Check r to make sure the insertion was successful.
                        if(r != 0)
                        {
                            // Do something if we couldn't insert...
                        }
                    }
                    // EXPERIMENTAL END
                }      // for vid_d1_adj ...

                // EXPERIMENTAL BEGIN
                // std::cout << "::: " << vid << " -> " << used_hash_size << std::endl;
                _degree_d2(vid) = used_hash_size;

                // Clear the Begins values.
                for(nnz_lno_type i = 0; i < globally_used_hash_count; i++)
                {
                    nnz_lno_type dirty_hash            = globally_used_hash_indices[ i ];
                    hash_map.hash_begins[ dirty_hash ] = -1;
                }
                // EXPERIMENTAL END

            }      // for ichunk ...
            _m_space.release_chunk(globally_used_hash_indices);
        }      // operator() (end)
    };         // struct functorCalculateD2Degree (end)


};      // end class GraphColorDistance2




/**
 * Compute Distance-2 Degree Stats
 *
 * Distance-2 Degree of a vertex, v, is the sum of the unique paths from u to v
 * where v is 2 hops from u. (i.e., paths like: `u -> * -> v`)
 *
 * This function calculates the distance-2 degree of all the vertices in the graph,
 * the maximum distance-2 degree.
 *
 * If the graph is symmetric, give the same value for col_map and row_map,
 * and for row_entries and col_entries.
 *
 * @param[in]  handle         The Kernel Handle
 * @param[in]  num_rows       Number of rows in the matrix (number of vertices)
 * @param[in]  num_cols       Number of columns in the matrix
 * @param[in]  row_map        Row map
 * @param[in]  row_entries    Row entries
 * @param[in]  col_map        Column map
 * @param[in]  col_entries    Column entries
 * @param[out] degree_d2_dist View to fill with distance-2 degree information.
 * @param[out] degree_d2_max  Maximum distance-2 degree found.
 *
 * @return Nothing
 *
 * Note: This is currently EXPERIMENTAL.
 */
template<class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename lno_col_view_t_, typename lno_colnnz_view_t_>
void graph_compute_distance2_degree(KernelHandle *handle,
                                    typename KernelHandle::nnz_lno_t num_rows,
                                    typename KernelHandle::nnz_lno_t num_cols,
                                    lno_row_view_t_ row_map,
                                    lno_nnz_view_t_ row_entries,
                                    lno_col_view_t_ col_map,
                                    lno_colnnz_view_t_ col_entries,
                                    typename KernelHandle::GraphColoringHandleType::non_const_1d_size_type_view_t& degree_d2_dist,
                                    size_t& degree_d2_max)
{
    typename KernelHandle::GraphColorDistance2HandleType *gch_d2 = handle->get_distance2_graph_coloring_handle();

    Impl::GraphColorDistance2<typename KernelHandle::GraphColorDistance2HandleType, lno_row_view_t_, lno_nnz_view_t_, lno_col_view_t_, lno_colnnz_view_t_>
        gc(num_rows, num_cols, row_entries.extent(0), row_map, row_entries, col_map, col_entries, gch_d2);

    gc.compute_distance2_degree(degree_d2_dist, degree_d2_max);
}




/**
 * Validate Distance 2 Graph Coloring
 *
 * If the graph is symmetric, give the same value for col_map and row_map,
 * and for row_entries and col_entries.
 *
 * @param[in]  handle           The kernel handle
 * @param[in]  num_rows         Number of rows in the matrix (number of vertices)
 * @param[in]  num_cols         Number of columns in the matrix
 * @param[in]  row_map          The row map
 * @param[in]  row_entries      The row entries
 * @param[in]  col_map          The column map
 * @param[in]  col_entries      The column entries
 * @param[out] validation_flags An array of 4 booleans.
 *                              validation_flags[0] : True IF the distance-2 coloring is invalid.
 *                              validation_flags[1] : True IF the coloring is bad because vertices are left uncolored.
 *                              validation_flags[2] : True IF the coloring is bad because at least one pair of vertices
 *                                                    at distance=2 from each other has the same color.
 *                              validation_flags[3] : True IF a vertex has a color greater than number of vertices in the graph.
 *                                                    May not be an INVALID coloring, but can indicate poor quality in coloring.
 *
 * @return boolean that is TRUE if the Distance-2 coloring is valid. False if otherwise.
 */
template<class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename lno_col_view_t_, typename lno_colnnz_view_t_>
bool graph_verify_distance2_color(KernelHandle *handle,
                                  typename KernelHandle::nnz_lno_t num_rows,
                                  typename KernelHandle::nnz_lno_t num_cols,
                                  lno_row_view_t_ row_map,
                                  lno_nnz_view_t_ row_entries,
                                  // If graph is symmetric, simply give same for col_map and row_map, and row_entries and col_entries.
                                  lno_col_view_t_ col_map,
                                  lno_colnnz_view_t_ col_entries,
                                  bool validation_flags[])
{
    bool output = true;

    typename KernelHandle::GraphColorDistance2HandleType *gch_d2 = handle->get_distance2_graph_coloring_handle();

    Impl::GraphColorDistance2<typename KernelHandle::GraphColorDistance2HandleType, lno_row_view_t_, lno_nnz_view_t_, lno_col_view_t_, lno_colnnz_view_t_>
        gc(num_rows, num_cols, row_entries.extent(0), row_map, row_entries, col_map, col_entries, gch_d2);

    output = gc.verify_coloring(row_map, row_entries, col_map, col_entries, gch_d2->get_vertex_colors(), validation_flags);

    return output;
}




/**
 * Prints out a histogram of graph colors for Distance-2 Graph Coloring
 *
 * If the graph is symmetric, give the same value for col_map and row_map,
 * and for row_entries and col_entries.
 *
 * @param[in]  handle           The kernel handle
 * @param[in]  num_rows         Number of rows in the matrix (number of vertices)
 * @param[in]  num_cols         Number of columns in the matrix
 * @param[in]  row_map          The row map
 * @param[in]  row_entries      The row entries
 * @param[in]  col_map          The column map
 * @param[in]  col_entries      The column entries
 * @param[out] validation_flags An array of 4 booleans.
 *                              validation_flags[0] : True IF the distance-2 coloring is invalid.
 *                              validation_flags[1] : True IF the coloring is bad because vertices are left uncolored.
 *                              validation_flags[2] : True IF the coloring is bad because at least one pair of vertices
 *                                                    at distance=2 from each other has the same color.
 *                              validation_flags[3] : True IF a vertex has a color greater than number of vertices in the graph.
 *                                                    May not be an INVALID coloring, but can indicate poor quality in coloring.
 * @param[in] csv               Output in CSV format? Default: false
 *
 * @return nothing
 */
template<class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename lno_col_view_t_, typename lno_colnnz_view_t_>
void graph_print_distance2_color_histogram(KernelHandle *handle,
                                           typename KernelHandle::nnz_lno_t num_rows,
                                           typename KernelHandle::nnz_lno_t num_cols,
                                           lno_row_view_t_ row_map,
                                           lno_nnz_view_t_ row_entries,
                                           // If graph is symmetric, simply give same for col_map and row_map, and row_entries and col_entries.
                                           lno_col_view_t_ col_map,
                                           lno_colnnz_view_t_ col_entries,
                                           bool csv=false)
{
    typename KernelHandle::GraphColorDistance2HandleType *gch_d2 = handle->get_distance2_graph_coloring_handle();

    Impl::GraphColorDistance2<typename KernelHandle::GraphColorDistance2HandleType, lno_row_view_t_, lno_nnz_view_t_, lno_col_view_t_, lno_colnnz_view_t_>
        gc(num_rows, num_cols, row_entries.extent(0), row_map, row_entries, col_map, col_entries, gch_d2);

    if(csv)
    {
        gc.print_color_histogram_csv();
    }
    else
    {
        gc.print_color_histogram();
    }
}




}      // namespace Impl
}      // namespace KokkosGraph


#endif      // _KOKKOSGRAPH_DISTANCE2COLOR_IMPL_HPP
