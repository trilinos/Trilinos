KokkosGraph::graph_color_distance2
##################################

Defined in header: :code:`KokkosGraph_Distance2Color.hpp`

.. code:: cppkokkos

  template <class KernelHandle, typename InRowmap, typename InEntries>
  void graph_color_distance2(KernelHandle *handle, typename KernelHandle::nnz_lno_t num_verts, InRowmap row_map,
                             InEntries row_entries);

Colors the vertices of a graph such that every vertex, its neighbors and its neighbors' neighbors have distict colors.

The graph must be symmetric, but it is not required to have diagonal entries. The coloring will not have distance-1 or distance-2 conflicts.

A view of length num_vertices, containing the colors will be returned through the handle: `handle->get_distance2_graph_coloring_handle()->get_vertex_colors()`

Parameters
==========

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` that stores algorithm parameters and the output colors.

:num_verts: the number of vertices in the graph.

:row_map: the graph row map.

:row_entries: the graph column indices.

Type Requirements
=================

No type requirements will be asserted.

..
   .. note::

      Obviously we should probably look at improving this "No requirement asserted"...

Example
=======
.. code-block:: c++

   crsMat A =
     KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat>(numVerts, numVerts, nnz, row_size_variance, bandwidth);
   auto G = A.graph;
   // Symmetrize the graph
   rowmap_t symRowmap;
   entries_t symEntries;
   KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<c_rowmap_t, c_entries_t, rowmap_t, entries_t, execution_space>(
      numVerts, G.row_map, G.entries, symRowmap, symEntries);
   std::vector<GraphColoringAlgorithmDistance2> algos = {COLORING_D2_DEFAULT, COLORING_D2_SERIAL,    COLORING_D2_VB,
                                                         COLORING_D2_VB_BIT,  COLORING_D2_VB_BIT_EF, COLORING_D2_NB_BIT};
   for (auto algo : algos) {
     KernelHandle kh;
     kh.create_distance2_graph_coloring_handle(algo);
     // Compute the Distance-2 graph coloring.
     graph_color_distance2<KernelHandle, c_rowmap_t, c_entries_t>(&kh, numVerts, symRowmap, symEntries);
     execution_space().fence();
     auto coloring_handle = kh.get_distance2_graph_coloring_handle();
     auto colors          = coloring_handle->get_vertex_colors();
     ...
   }

