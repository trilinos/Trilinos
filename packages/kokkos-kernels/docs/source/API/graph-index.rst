API: Graph
##########

.. toctree::
   :maxdepth: 2
   :hidden:

   graph/distance1_color
   graph/distance2_color
   graph/rcb
   graph/merge
   graph/merge_path
   graph/load_balance

Graph
=====

These algorithms provide generic graph coloring capabilities used by Gauss-Seidel, multigrid aggregation, etc... We distinguish to main categories of algorithms, distance 1 coloring and distance 2 coloring. We also provide a coloring handle that allows users to easily control the behavior of the algorithms. We also include an implementation of the recursive partitioning bisection (RCB) algorithm.

Graph coloring handle
=====================

The graph coloring handle is extending the KokkosKernels handle to provide options that are specific to coloring algorithms.

Distance 1 graph coloring
=========================

Distance 1 coloring algorithms will ensure that each node has a different color than all of its neighbors.

- :doc:`Graph Coloring <graph/distance1_color>`

Distance 2 and One-sided Bipartite graph coloring
=================================================

Distance 2 coloring algorithms will ensure that each node has a different color than its neighbors and its neighbors' neighbors.

- :doc:`Distance-2 Graph Coloring <graph/distance2_color>`

Recursive Coordinate Bisection (RCB)
====================================

RCB performs recursive partitioning on a set of coordinates of the mesh points.

- :doc:`Recursive Coordinate Bisection <graph/rcb>`

Ordered Sequence Utilities
==========================

These public interfaces operate on ordered rank-1 inputs and are used to build
merge-based graph and sparse algorithms.

- :doc:`KokkosGraph merge routines <graph/merge>`
- :doc:`KokkosGraph merge-path routines <graph/merge_path>`
- :doc:`KokkosGraph load-balance routines <graph/load_balance>`
