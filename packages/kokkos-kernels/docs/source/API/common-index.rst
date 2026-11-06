API: Common
###########

.. toctree::
   :maxdepth: 2
   :hidden:

   common/lower_bound
   common/upper_bound
   common/exclusive_parallel_prefix_sum
   common/inclusive_parallel_prefix_sum

Common
======

The common API provides low-level utilities that are reused across multiple
Kokkos Kernels components. The current public interfaces in this section expose
portable ordered-search operations on rank-1 view-like inputs.

Lower Bound Search
==================

Lower-bound search returns the first index whose entry is not ordered before a
target value.

- :doc:`KokkosKernels::lower_bound_thread and KokkosKernels::lower_bound_team <common/lower_bound>`

Upper Bound Search
==================

Upper-bound search returns the first index whose entry is ordered after a
target value.

- :doc:`KokkosKernels::upper_bound_thread and KokkosKernels::upper_bound_team <common/upper_bound>`

Prefix Sums
===========

These functions perform in-place exclusive or inclusive parallel prefix sums. There are two
reasons to use these instead of ``Kokkos::Experimental::exclusive_scan`` or ``Kokkos::Experimental::inclusive_scan``. First,
``KokkosKernels::exclusive_parallel_prefix_sum`` can optionally compute the "final sum", or the sum of all elements on input.
Second, all other versions (those not computing the final sum) are non-blocking.

- :doc:`KokkosKernels::exclusive_parallel_prefix_sum <common/exclusive_parallel_prefix_sum>`
- :doc:`KokkosKernels::inclusive_parallel_prefix_sum <common/inclusive_parallel_prefix_sum>`
