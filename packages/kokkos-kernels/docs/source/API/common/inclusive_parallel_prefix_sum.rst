KokkosKernels::inclusive_parallel_prefix_sum
############################################

Defined in header: :code:`KokkosKernels_SimpleUtils.hpp`

.. code:: c++

  template <typename ExecSpace, typename view_t>
  void inclusive_parallel_prefix_sum(const ExecSpace& exec, const view_t& arr);

  template <typename view_t>
  void inclusive_parallel_prefix_sum(const view_t& arr);

These functions perform an in-place inclusive parallel prefix sum on the rank-1 ``Kokkos::View`` ``arr``.
That is, element i will be replaced with the sum of elements 0...i, including i itself.
If ``exec`` is provided, the underlying parallel scan is executed on that instance.

Parameters
==========

:arr: rank-1 ``Kokkos::View`` to prefix sum.

:exec: execution space instance on which to run the prefix sum.

Type Requirements
-----------------

- ``view_t`` must be a rank-1 ``Kokkos::View`` type.

- ``ExecSpace`` must be a Kokkos execution space.

Notes
=====

- Both versions of this function are non-blocking.
- If no ``exec`` is provided, the kernel is run on the default instance of ``view_t::execution_space``.
