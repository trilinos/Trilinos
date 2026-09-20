KokkosKernels::exclusive_parallel_prefix_sum
############################################

Defined in header: :code:`KokkosKernels_SimpleUtils.hpp`

.. code:: c++

  template <typename ExecSpace, typename view_t>
  void exclusive_parallel_prefix_sum(const ExecSpace &exec, const view_t& arr);

  template <typename view_t>
  void exclusive_parallel_prefix_sum(const view_t& arr);

  template <typename ExecSpace, typename view_t>
  void exclusive_parallel_prefix_sum(const ExecSpace &exec, const view_t& arr,
                                     typename view_t::non_const_value_type &finalSum);

  template <typename view_t>
  void exclusive_parallel_prefix_sum(const view_t& arr,
                                     typename view_t::non_const_value_type &finalSum);

These functions perform an in-place exclusive parallel prefix sum on the rank-1 ``Kokkos::View`` ``arr``.
That is, element i will be replaced with the sum of elements 0...i-1 (not including i itself in the sum).
If ``exec`` is provided, the underlying parallel scan is executed on that instance.
If ``finalSum`` (output parameter) is provided, the total sum of all elements of ``arr`` on input is computed and stored here.

Parameters
==========

:arr: rank-1 ``Kokkos::View`` to prefix sum.

:exec: execution space instance on which to run the prefix sum.

:finalSum: output parameter for computing the sum of all elements originally contained in ``arr``.

Type Requirements
-----------------

- ``view_t`` must be a rank-1 ``Kokkos::View`` type.

- ``ExecSpace`` must be a Kokkos execution space.

Notes
=====

- Versions 1 and 2 (those without the ``finalSum`` parameter) are non-blocking.
- Versions 3 and 4 are blocking; ``finalSum`` is valid immediately when the function returns.
- If no ``exec`` is provided, the kernel is run on the default instance of ``view_t::execution_space``.
