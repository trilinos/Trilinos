#ifndef _KOKKOS_GRAPH_COLOR_HPP
#define _KOKKOS_GRAPH_COLOR_HPP

#include "GraphColor_impl.hpp"
#include "GraphColoringHandle.hpp"
namespace Experimental{

namespace KokkosKernels{
namespace Graph{

template <class KernelHandle>
void graph_color_symbolic(
    KernelHandle *handle,
    typename KernelHandle::idx_array_type row_map,
    typename KernelHandle::idx_edge_array_type entries){

  Kokkos::Impl::Timer timer;

  //typename KernelHandle::idx_array_type row_map = handle->get_row_map();
  //typename KernelHandle::idx_edge_array_type entries = handle->get_entries();

  typename KernelHandle::GraphColoringHandleType *gch = handle->get_graph_coloring_handle();

  Experimental::KokkosKernels::Graph::ColoringAlgorithm algorithm = gch->get_coloring_type();

  typedef typename KernelHandle::GraphColoringHandleType::color_array_type color_view_type;
  color_view_type colors_out = color_view_type("Graph Colors", row_map.dimension_0() - 1);

  typedef typename Experimental::KokkosKernels::Graph::Impl::GraphColor
      <typename KernelHandle::GraphColoringHandleType> BaseGraphColoring;
  BaseGraphColoring *gc = NULL;



  switch (algorithm){
  case Experimental::KokkosKernels::Graph::COLORING_SERIAL:

    gc = new BaseGraphColoring(
        row_map.dimension_0() - 1, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case Experimental::KokkosKernels::Graph::COLORING_VB:
  case Experimental::KokkosKernels::Graph::COLORING_VBBIT:
  case Experimental::KokkosKernels::Graph::COLORING_VBCS:

    typedef typename Experimental::KokkosKernels::Graph::Impl::GraphColor_VB
        <typename KernelHandle::GraphColoringHandleType> VBGraphColoring;
    gc = new VBGraphColoring(
        row_map.dimension_0() - 1, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case Experimental::KokkosKernels::Graph::COLORING_EB:

    typedef typename Experimental::KokkosKernels::Graph::Impl::GraphColor_EB
        <typename KernelHandle::GraphColoringHandleType> EBGraphColoring;

    gc = new EBGraphColoring(row_map.dimension_0() - 1, entries.dimension_0(),row_map, entries, gch);
    break;
  case Experimental::KokkosKernels::Graph::COLORING_DEFAULT:
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

/*
template <class KernelHandle>
void graph_color_numeric(
    KernelHandle *handle){
  if (!handle->get_graph_coloring_handle()->is_coloring_called()){
    graph_color_symbolic<KernelHandle>(
        handle);
  }
  else{
    handle->get_graph_coloring_handle()->add_to_overall_coloring_time(0);
    handle->get_graph_coloring_handle()->set_coloring_time(0);
  }

}

template <class KernelHandle>
void graph_color_solve(
    KernelHandle *handle
){

  if (!handle->get_graph_coloring_handle()->is_coloring_called()){
    graph_color_symbolic<KernelHandle>(
        handle);
  }  else{
    handle->get_graph_coloring_handle()->add_to_overall_coloring_time(0);
    handle->get_graph_coloring_handle()->set_coloring_time(0);
  }
}
*/

}
}
}

#endif//_KOKKOS_GRAPH_COLOR_HPP
