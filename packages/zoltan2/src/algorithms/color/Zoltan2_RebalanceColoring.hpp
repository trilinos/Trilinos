// Sketch: Does not compile nor work yet!

// Given a valid coloring of a graph, rebalance the colors such that
// every color class has at least min_size vertices.
// This function can be called as a post-processing after any initial
// coloring.
   void rebalanceColoring(
     const lno_t nVtx,
     ArrayView<const lno_t> edgeIds,
     ArrayView<const lno_t> offsets,
     ArrayRCP<int> colors,
     lno_t thres
     )
   {
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL
    Z2_THROW_EXPERIMENTAL("Zoltan2 rebalanceColoring is experimental and not tested");
#else
   // Count size of each color class.
   // For each color i with < min_size vertices, pair it with
   // a color j with > min_size vertices. Recolor vertices with color j to i.

#endif
   }
