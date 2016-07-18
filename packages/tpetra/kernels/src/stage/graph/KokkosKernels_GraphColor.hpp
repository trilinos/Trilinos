#ifndef _KOKKOS_GRAPH_COLOR_HPP
#define _KOKKOS_GRAPH_COLOR_HPP

#include "KokkosKernels_GraphColor_impl.hpp"
#include "KokkosKernels_GraphColorHandle.hpp"
#include "KokkosKernels_Utils.hpp"
namespace KokkosKernels{

namespace Experimental{


namespace Graph{


template <class KernelHandle,typename lno_row_view_t_, typename lno_nnz_view_t_>
void graph_color_symbolic(
    KernelHandle *handle,
    typename KernelHandle::row_lno_t num_rows,
    typename KernelHandle::row_lno_t num_cols,
    lno_row_view_t_ row_map,
    lno_nnz_view_t_ entries,
    bool is_symmetric = true){

  Kokkos::Impl::Timer timer;

  typename KernelHandle::GraphColoringHandleType *gch = handle->get_graph_coloring_handle();

  ColoringAlgorithm algorithm = gch->get_coloring_algo_type();

  typedef typename KernelHandle::GraphColoringHandleType::color_view_t color_view_type;

  color_view_type colors_out = color_view_type("Graph Colors", num_rows);


  typedef typename Impl::GraphColor
      <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> BaseGraphColoring;
  BaseGraphColoring *gc = NULL;


  switch (algorithm){
  case COLORING_SERIAL:

    gc = new BaseGraphColoring(
        num_rows, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case COLORING_SERIAL2:

    gc = new Impl::GraphColor2<typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_>(
        num_rows, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case COLORING_VB:
  case COLORING_VBBIT:
  case COLORING_VBCS:

    typedef typename Impl::GraphColor_VB
        <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> VBGraphColoring;
    gc = new VBGraphColoring(
        num_rows, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case COLORING_EB:

    typedef typename Impl::GraphColor_EB
        <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> EBGraphColoring;

    gc = new EBGraphColoring(num_rows, entries.dimension_0(),row_map, entries, gch);
    break;
  case COLORING_DEFAULT:
    break;

  }

  int num_phases = 0;
  gc->color_graph(colors_out, num_phases);
  /*
  switch (gch->get_coloring_type()){
    default:
    case KokkosKernels::Experimental::Graph::Distance1:

      break;
    case KokkosKernels::Experimental::Graph::Distance2:
      gc->d2_color_graph(colors_out, num_phases);
      break;
  }
  */

  delete gc;
  double coloring_time = timer.seconds();
  gch->add_to_overall_coloring_time(coloring_time);
  gch->set_coloring_time(coloring_time);
  gch->set_num_phases(num_phases);
  gch->set_vertex_colors(colors_out);
}

template <class KernelHandle,
          typename lno_row_view_t_, typename lno_nnz_view_t_,
          typename lno_col_view_t_, typename lno_colnnz_view_t_>
void d2_graph_color(
    KernelHandle *handle,
    typename KernelHandle::row_lno_t num_rows,
    typename KernelHandle::row_lno_t num_cols,
    lno_row_view_t_ row_map,
    lno_nnz_view_t_ row_entries,
    lno_col_view_t_ col_map, //if graph is symmetric, simply give same for col_map and row_map, and row_entries and col_entries.
    lno_colnnz_view_t_ col_entries){

  Kokkos::Impl::Timer timer;
  typename KernelHandle::GraphColoringHandleType *gch = handle->get_graph_coloring_handle();

  ColoringAlgorithm algorithm = gch->get_coloring_algo_type();

  typedef typename KernelHandle::GraphColoringHandleType::color_view_t color_view_type;

  color_view_type colors_out = color_view_type("Graph Colors", num_rows);


  typedef typename Impl::GraphColor <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> BaseGraphColoring;
  //BaseGraphColoring *gc = NULL;

  int num_phases = 0;

  switch (algorithm){
  case COLORING_SERIAL:
  {
    BaseGraphColoring gc(
        num_rows, row_entries.dimension_0(),
        row_map, row_entries, gch);
    gc.d2_color_graph/*<lno_col_view_t_,lno_colnnz_view_t_>*/(colors_out, num_phases, num_cols, col_map, col_entries);
    break;
  }
  case COLORING_SERIAL2:
  {

    Impl::GraphColor2<typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> gc(
        num_rows, row_entries.dimension_0(),
        row_map, row_entries, gch);
    gc.d2_color_graph(colors_out, num_phases, num_cols, col_map, col_entries);
    break;
  }
  case COLORING_VB:
  case COLORING_VBBIT:
  case COLORING_VBCS:
  {
    typedef typename Impl::GraphColor_VB
        <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> VBGraphColoring;
    VBGraphColoring gc(
        num_rows, row_entries.dimension_0(),
        row_map, row_entries, gch);
    gc.d2_color_graph(colors_out, num_phases, num_cols, col_map, col_entries);

  }
  break;
  case COLORING_EB:
  {
    typedef typename Impl::GraphColor_EB
        <typename KernelHandle::GraphColoringHandleType, lno_row_view_t_, lno_nnz_view_t_> EBGraphColoring;

    EBGraphColoring gc(num_rows, row_entries.dimension_0(),row_map, row_entries, gch);
    gc.d2_color_graph(colors_out, num_phases, num_cols, col_map, col_entries);

  }
  break;
  case COLORING_DEFAULT:
    break;

  }



  double coloring_time = timer.seconds();
  gch->add_to_overall_coloring_time(coloring_time);
  gch->set_coloring_time(coloring_time);
  gch->set_num_phases(num_phases);
  gch->set_vertex_colors(colors_out);
}
}
}
}

#endif//_KOKKOS_GRAPH_COLOR_HPP
