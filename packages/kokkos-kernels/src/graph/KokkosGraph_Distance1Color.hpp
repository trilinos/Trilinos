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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef _KOKKOSGRAPH_DISTANCE1_COLOR_HPP
#define _KOKKOSGRAPH_DISTANCE1_COLOR_HPP

#include <sstream>

#include "KokkosGraph_Distance1ColorHandle.hpp"
#include "KokkosGraph_Distance1Color_impl.hpp"
#include "KokkosKernels_Utils.hpp"

namespace KokkosGraph{

namespace Experimental{



template <class KernelHandle,typename lno_row_view_t_, typename lno_nnz_view_t_>
void graph_color_symbolic(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t num_rows,
    typename KernelHandle::nnz_lno_t /* num_cols */,
    lno_row_view_t_ row_map,
    lno_nnz_view_t_ entries,
    bool /* is_symmetric */ = true){

  Kokkos::Impl::Timer timer;

  typename KernelHandle::GraphColoringHandleType *gch = handle->get_graph_coloring_handle();

  ColoringAlgorithm algorithm = gch->get_coloring_algo_type();

  typedef typename KernelHandle::GraphColoringHandleType::color_view_t color_view_type;

  gch->set_tictoc(handle->get_verbose());

  color_view_type colors_out = color_view_type("Graph Colors", num_rows);


  typedef typename Impl::GraphColor
      <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> BaseGraphColoring;
  BaseGraphColoring *gc = NULL;

  switch (algorithm){
  case COLORING_SERIAL:
    gc = new BaseGraphColoring(num_rows, entries.extent(0), row_map, entries, gch);
    break;

  case COLORING_VB:
  case COLORING_VBBIT:
  case COLORING_VBCS:
    typedef typename Impl::GraphColor_VB <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> VBGraphColoring;
    gc = new VBGraphColoring(num_rows, entries.extent(0), row_map, entries, gch);
    break;

  case COLORING_VBD:
  case COLORING_VBDBIT:
    typedef typename Impl::GraphColor_VBD <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> VBDGraphColoring;
    gc = new VBDGraphColoring(num_rows, entries.extent(0), row_map, entries, gch);
    break;

  case COLORING_EB:
    typedef typename Impl::GraphColor_EB <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> EBGraphColoring;
    gc = new EBGraphColoring(num_rows, entries.extent(0),row_map, entries, gch);
    break;

  case COLORING_DEFAULT:
    break;

  default:
    break;
  }

  int num_phases = 0;
  gc->color_graph(colors_out, num_phases);

  delete gc;
  double coloring_time = timer.seconds();
  gch->add_to_overall_coloring_time(coloring_time);
  gch->set_coloring_time(coloring_time);
  gch->set_num_phases(num_phases);
  gch->set_vertex_colors(colors_out);
}



template <class KernelHandle,typename lno_row_view_t_, typename lno_nnz_view_t_>
void graph_color(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t num_rows,
    typename KernelHandle::nnz_lno_t num_cols,
    lno_row_view_t_ row_map,
    lno_nnz_view_t_ entries,
    bool is_symmetric = true)
{
  graph_color_symbolic(handle, num_rows, num_cols, row_map, entries, is_symmetric);
}


// Distance 2 graph coloring -- serial only --
// todo: move this over to the Distance2 handle
template <class KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_,
          typename lno_col_view_t_, typename lno_colnnz_view_t_>
void graph_compute_distance2_color_serial(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t num_rows,
    typename KernelHandle::nnz_lno_t num_cols,
    lno_row_view_t_ row_map,
    lno_nnz_view_t_ row_entries,
    //if graph is symmetric, simply give same for col_map and row_map, and row_entries and col_entries.
    lno_col_view_t_ col_map,
    lno_colnnz_view_t_ col_entries)
{

  Kokkos::Impl::Timer timer;
  typename KernelHandle::GraphColoringHandleType *gch = handle->get_graph_coloring_handle();

  ColoringAlgorithm algorithm = gch->get_coloring_algo_type();

  typedef typename KernelHandle::GraphColoringHandleType::color_view_t color_view_type;

  // typedef typename Impl::GraphColor<typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_>BaseGraphColoring;

  int num_phases = 0;

  switch (algorithm)
  {
    #if 0   // SERIAL shouldn't be in the d2 coloring function since you can get it from graph_color()
    case COLORING_SERIAL:
    {
      color_view_type colors_out = color_view_type("Graph Colors", num_rows);
      BaseGraphColoring gc(num_rows, row_entries.extent(0), row_map, row_entries, gch);
      gc.d2_color_graph/*<lno_col_view_t_,lno_colnnz_view_t_>*/(colors_out, num_phases, num_cols, col_map, col_entries);
      gch->set_num_phases(num_phases);
      gch->set_vertex_colors(colors_out);
      break;
    }
    #endif

    // todo: Remove SERIAL Distance-2 Graph Coloring from this interface
    //
    case COLORING_SERIAL2:
    {
      color_view_type colors_out = color_view_type("Graph Colors", num_rows);
      Impl::GraphColor2<typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_>
          gc(num_rows, row_entries.extent(0), row_map, row_entries, gch);
      gc.d2_color_graph(colors_out, num_phases, num_cols, col_map, col_entries);
      gch->set_num_phases(num_phases);
      gch->set_vertex_colors(colors_out);
      break;
    }

    default:
    {
      color_view_type colors_out = color_view_type("Graph Colors", num_rows);
      Impl::GraphColor2<typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_>
          gc(num_rows, row_entries.extent(0), row_map, row_entries, gch);
      gc.d2_color_graph(colors_out, num_phases, num_cols, col_map, col_entries);
      gch->set_num_phases(num_phases);
      gch->set_vertex_colors(colors_out);
      break;
    }
  }

  double coloring_time = timer.seconds();
  gch->add_to_overall_coloring_time(coloring_time);
  gch->set_coloring_time(coloring_time);
}


#if defined(KOKKOS_ENABLE_DEPRECATED_CODE)
template <class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename lno_col_view_t_, typename lno_colnnz_view_t_>
void d2_graph_color(KernelHandle *handle,
                    typename KernelHandle::nnz_lno_t num_rows,
                    typename KernelHandle::nnz_lno_t num_cols,
                    lno_row_view_t_ row_map,
                    lno_nnz_view_t_ row_entries,
                    //if graph is symmetric, simply give same for col_map and row_map, and row_entries and col_entries.
                    lno_col_view_t_ col_map,
                    lno_colnnz_view_t_ col_entries)
{
    graph_compute_distance2_color_serial(handle, num_rows, num_cols, row_map, row_entries, col_map, col_entries);
    std::cout << "Deprecation warning: d2_graph_color() will be replaced with graph_compute_distance2_color_serial()" << std::endl;
}
#endif


}  // end namespace Experimental
}  // end namespace KokkosGraph

#endif   // _KOKKOSGRAPH_DISTANCE1_COLOR_HPP

