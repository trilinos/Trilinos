KokkosGraph::graph_color_symbolic
#################################

Defined in header: :code:`KokkosGraph_Distance1Color.hpp`

.. code:: cppkokkos

  template <class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void graph_color_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t num_rows,
                            typename KernelHandle::nnz_lno_t /* num_cols */, lno_row_view_t_ row_map,
                            lno_nnz_view_t_ entries, bool /* is_symmetric */ = true);

Colors the vertices of a graph such that every vertex and its neighbors have distinc colors.

.. math::

   \text{Given a graph}\ G=(\mathcal{V}, \mathcal{E})\\
   \forall v\in\mathcal{V}, \forall w\in neigh(v),\ color(v) != color(w)

..
   .. note::

      We could consider adding a stream interface and a staticcrsgraph interface as well.
      Should we also deprecate the ``graph_color`` interface as well as the current interface and move to new API without ``num_cols`` and ``is_symmetric``? Also technically ``num_rows`` is a bit redundant as ``row_map.extent(0) - 1`` should be the same?

Parameters
==========

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` that stores algorithm parameters and the output colors.

:num_rows: the number of vertices in the graph.

:row_map: the graph row map.

:entries: the graph column indices.

:num_cols, is_symmetric: these two parameters are ignored and are only present for backward compatibility purposes.

Type Requirements
=================

No type requirements will be asserted.

..
   .. note::

      Obviously we should probably look at improving this "No requirement asserted"...

Example
=======

.. literalinclude:: ../../../../example/wiki/graph/KokkosGraph_wiki_coloring.cpp
  :language: c++
  :lines: 16-
