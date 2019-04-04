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
#include <iomanip>
#include <vector>

#include <KokkosSparse_spgemm.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "KokkosGraph_Distance1ColorHandle.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosKernels_Handle.hpp"

#ifndef _KOKKOSCOLORINGD2MATRIXSQUAREDIMP_HPP
#define _KOKKOSCOLORINGD2MATRIXSQUAREDIMP_HPP


namespace KokkosGraph {

namespace Impl {


/*! \brief Base class for graph coloring purposes.
 *  Each color represents the set of the vertices that are independent,
 *  e.g. no vertex having same color shares an edge.
 *  General aim is to find the minimum number of colors, minimum number of independent sets.
 */
template<typename HandleType,
         typename lno_row_view_t_,
         typename lno_nnz_view_t_,
         typename clno_row_view_t_,
         typename clno_nnz_view_t_>
class GraphColorDistance2MatrixSquared
{
  public:
    using in_lno_row_view_type              = lno_row_view_t_;
    using in_lno_nnz_view_type              = lno_nnz_view_t_;
    using color_type                        = typename HandleType::GraphColorDistance2HandleType::color_type;
    using color_view_type                   = typename HandleType::GraphColorDistance2HandleType::color_view_type;
    using size_type                         = typename HandleType::size_type;
    using const_size_type                   = typename HandleType::const_size_type;
    using nnz_lno_type                      = typename HandleType::nnz_lno_t;
    using MyExecSpace                       = typename HandleType::HandleExecSpace;
    using MyTempMemorySpace                 = typename HandleType::HandleTempMemorySpace;
    using row_lno_view_device_type          = typename lno_row_view_t_::device_type;
    using const_lno_row_view_type           = typename lno_row_view_t_::const_type;
    using const_lno_nnz_view_type           = typename lno_nnz_view_t_::const_type;
    using non_const_lno_nnz_view_type       = typename lno_nnz_view_t_::non_const_type;
    using const_clno_row_view_type          = typename clno_row_view_t_::const_type;
    using const_clno_nnz_view_type          = typename clno_nnz_view_t_::const_type;
    using non_const_clno_nnz_view_type      = typename clno_nnz_view_t_::non_const_type;
    using size_type_temp_work_view_type     = typename HandleType::size_type_temp_work_view_t;
    using scalar_temp_work_view_type        = typename HandleType::scalar_temp_work_view_t;
    using nnz_lno_persistent_work_view_type = typename HandleType::nnz_lno_persistent_work_view_t;
    using nnz_lno_temp_work_view_type       = typename HandleType::nnz_lno_temp_work_view_t;
    using single_dim_index_view_type        = typename Kokkos::View<nnz_lno_type, row_lno_view_device_type>;
    using my_exec_space                     = Kokkos::RangePolicy<MyExecSpace>;


  protected:
    nnz_lno_type             nr;          // num_rows  (# verts)
    nnz_lno_type             nc;          // num cols
    size_type                ne;          // # edges
    const_lno_row_view_type  xadj;        // rowmap, transpose of rowmap
    const_lno_nnz_view_type  adj;         // entries, transpose of entries   (size = # edges)
    const_clno_row_view_type t_xadj;      // rowmap, transpose of rowmap
    const_clno_nnz_view_type t_adj;       // entries, transpose of entries
    nnz_lno_type             nv;          // num vertices
    HandleType*              handle;      // the handle.

    bool verbose;

  public:
    /**
     * \brief GraphColor constructor.
     *
     * \param nv_: number of vertices in the graph
     * \param ne_: number of edges in the graph
     * \param row_map: the xadj array of the graph. Its size is nv_ +1
     * \param entries: adjacency array of the graph. Its size is ne_
     * \param handle: GraphColoringHandle object that holds the specification about the graph coloring,
     *    including parameters.
     */
    GraphColorDistance2MatrixSquared(nnz_lno_type             nr_,
                                     nnz_lno_type             nc_,
                                     size_type                ne_,
                                     const_lno_row_view_type  row_map,
                                     const_lno_nnz_view_type  entries,
                                     const_clno_row_view_type t_row_map,
                                     const_clno_nnz_view_type t_entries,
                                     HandleType*              handle_)
        : nr(nr_)
        , nc(nc_)
        , ne(ne_)
        , xadj(row_map)
        , adj(entries)
        , t_xadj(t_row_map)
        , t_adj(t_entries)
        , nv(nr_)
        , handle(handle_)
        , verbose(handle_->get_verbose())
    {
    }


    /** \brief GraphColor destructor.
     */
    virtual ~GraphColorDistance2MatrixSquared() {}


    /**
     * \brief Function to color the vertices of the graphs. This is the base class,
     * therefore, it only performs sequential coloring on the host device, ignoring the execution space.
     *
     * \param colors is the output array corresponding the color of each vertex.Size is this->nv.
     *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
     *   it should be all initialized with zeros. Any positive value in the given array, will make the
     *   algorithm to assume that the color is fixed for the corresponding vertex.
     * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
     */
    virtual void compute_distance2_color()
    {
        double              time = 0;
        Kokkos::Impl::Timer timer;
        timer.reset();

        std::string algName = "SPGEMM_KK_MEMSPEED";
        handle->create_spgemm_handle(KokkosSparse::StringToSPGEMMAlgorithm(algName));

        size_type_temp_work_view_type cRowptrs("cRowptrs", nr + 1);

        // Call symbolic multiplication of graph with itself (no transposes, and A and B are the same)
        KokkosSparse::Experimental::spgemm_symbolic(handle, nr, nc, nr, xadj, adj, false, t_xadj, t_adj, false, cRowptrs);

        // Get num nz in C
        auto Cnnz = handle->get_spgemm_handle()->get_c_nnz();

        // Must create placeholder value views for A and C (values are meaningless)
        // Said above that the scalar view type is the same as the colinds view type
        scalar_temp_work_view_type aFakeValues("A/B placeholder values (meaningless)", adj.size());

        // Allocate C entries array, and placeholder values
        nnz_lno_persistent_work_view_type cColinds("C colinds", Cnnz);
        scalar_temp_work_view_type        cFakeValues("C placeholder values (meaningless)", Cnnz);

        // Run the numeric kernel
        KokkosSparse::Experimental::spgemm_numeric(
          handle, nr, nc, nr, xadj, adj, aFakeValues, false, t_xadj, t_adj, aFakeValues, false, cRowptrs, cColinds, cFakeValues);

        if(this->verbose)
        {
            time = timer.seconds();
            std::cout << "\tTime Phase Square Matrix : " << time << std::endl << std::endl;
            this->handle->get_distance2_graph_coloring_handle()->add_to_overall_coloring_time_phase4(time);
            timer.reset();
        }

        // done with spgemm
        handle->destroy_spgemm_handle();

        // Now run distance-1 graph coloring on C, use LocalOrdinal for storing colors
        handle->create_graph_coloring_handle();
        KokkosGraph::Experimental::graph_color(
          handle, nr, nr, /*(const_rowptrs_view)*/ cRowptrs, /*(const_colinds_view)*/ cColinds);

        if(this->verbose)
        {
            time = timer.seconds();
            std::cout << "\tTime Phase Graph Coloring : " << time << std::endl << std::endl;
            this->handle->get_distance2_graph_coloring_handle()->add_to_overall_coloring_time_phase5(time);
            timer.reset();
        }

        // extract the colors
        // auto coloringHandle = handle->get_graph_coloring_handle();
        // color_view_type colorsDevice = coloringHandle->get_vertex_colors();

        // std::cout << "Num phases: " << handle->get_graph_coloring_handle()->get_num_phases() << std::endl;
        this->handle->get_distance2_graph_coloring_handle()->set_num_phases(
          this->handle->get_graph_coloring_handle()->get_num_phases());
        this->handle->get_distance2_graph_coloring_handle()->set_vertex_colors(
          this->handle->get_graph_coloring_handle()->get_vertex_colors());

        // clean up coloring handle
        handle->destroy_graph_coloring_handle();
    }

};      // GraphColorDistance2MatrixSquared (end)



}      // namespace Impl
}      // namespace KokkosGraph


#endif      // _KOKKOSCOLORINGD2MATRIXSQUAREDIMP_HPP
