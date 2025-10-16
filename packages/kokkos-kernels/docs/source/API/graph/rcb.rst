KokkosGraph::Experimental::recursive_coordinate_bisection
#########################################################

Defined in header: :code:`KokkosGraph_RCB.hpp`

.. code:: cppkokkos

  template <typename coors_view_type, typename perm_view_type>
  std::vector<typename perm_view_type::value_type> recursive_coordinate_bisection(coors_view_type &coordinates, perm_view_type &perm,
                                                                                  perm_view_type &reverse_perm, const int &n_levels);						 

Performs the recursive coordinate bisection (RCB) algorithm to partition a graph according to the coordinates of the mesh points.

This function will

1. return a vector containing sizes of partitions,
2. reorder the coordinate list to the RCB order,
3. update the permutation array describing the mapping from the original order to RCB order, and
4. update the reverse permutation array describing the mapping from the RCB order to the original order.

The function will throw a runtime exception if any of the following conditions are not met:

- ``coordinates`` are not 1-D, 2-D, or 3-D coordinates (i.e., ``coordinates.extent(1)`` is greater than 3)
- ``n_levels`` is smaller than 2

Parameters
==========

:coordinates: 1/2/3-D coordinates of the mesh points.

:perm: 1-D array describing the mapping from the original order to RCB order.

:reverse_perm: 1-D array describing the mapping from the RCB order to the original order.

:n_levels: the number of bisection levels.

Type Requirements
=================

- ``coors_view_type`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2

- ``perm_view_type`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1

Example
=======

Will be added later.

