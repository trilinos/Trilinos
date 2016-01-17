#ifndef _KOKKOS_GRAPH_COLOR_HPP
#define _KOKKOS_GRAPH_COLOR_HPP

#include "GraphColor_impl.hpp"
#include "GraphColoringHandle.hpp"
#include "KokkosKernelsUtils.hpp"
namespace KokkosKernels{

namespace Experimental{


namespace Graph{

template <class KernelHandle>
void graph_color_symbolic(
    KernelHandle *handle,
    typename KernelHandle::idx num_rows,
    typename KernelHandle::idx num_cols,
    typename KernelHandle::idx_array_type row_map,
    typename KernelHandle::idx_edge_array_type entries,
    bool is_symmetric = true){

  Kokkos::Impl::Timer timer;

  //typename KernelHandle::idx_array_type row_map = handle->get_row_map();
  //typename KernelHandle::idx_edge_array_type entries = handle->get_entries();

  typename KernelHandle::GraphColoringHandleType *gch = handle->get_graph_coloring_handle();

  ColoringAlgorithm algorithm = gch->get_coloring_type();

  typedef typename KernelHandle::GraphColoringHandleType::color_array_type color_view_type;
  color_view_type colors_out = color_view_type("Graph Colors", num_rows);

  typedef typename Impl::GraphColor
      <typename KernelHandle::GraphColoringHandleType> BaseGraphColoring;
  BaseGraphColoring *gc = NULL;


  switch (algorithm){
  case COLORING_SERIAL:

    gc = new BaseGraphColoring(
        num_rows, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case COLORING_VB:
  case COLORING_VBBIT:
  case COLORING_VBCS:

    typedef typename Impl::GraphColor_VB
        <typename KernelHandle::GraphColoringHandleType> VBGraphColoring;
    gc = new VBGraphColoring(
        num_rows, entries.dimension_0(),
        row_map, entries, gch);
    break;
  case COLORING_EB:

    typedef typename Impl::GraphColor_EB
        <typename KernelHandle::GraphColoringHandleType> EBGraphColoring;

    gc = new EBGraphColoring(num_rows, entries.dimension_0(),row_map, entries, gch);
    break;
  case COLORING_DEFAULT:
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

}
}
}

#endif//_KOKKOS_GRAPH_COLOR_HPP
