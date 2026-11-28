KokkosBlas::sum
###############

Defined in header: :code:`KokkosBlas1_sum.hpp`

.. code:: c++

  template <class execution_space, class XVector>
  requires Kokkos::is_execution_space_v<execution_space>
  typename XVector::non_const_value_type sum(const execution_space& space, const XVector& x);

  template <class XVector>
  typename XVector::non_const_value_type sum(const XVector& x);

  template <class execution_space, class RV, class XMV>
  void sum(const execution_space& space, const RV& R, const XMV& X);

  template <class RV, class XMV>
  requires Kokkos::is_view_v<RV>
  void sum(const RV& R, const XMV& X);

computes the sum of the elements in vector(s) `X` and stores the result in `R`. When the sum is returned by the function, the function is blocking and will complete before returning.

.. math::

   R_i = \sum_j X_{ij}

1. Computes the sum of the entries in `X` using ``space`` instance and returns them.
2. Computes the sum of the entries in `X` using the default instance of ``typename XVector::execution_space`` and returns them.
3. Computes the sum of the entries in `X` column wise and store the results in `R` using ``space`` instance.
4. Computes the sum of the entries in `X` column wise and store the results in `R` using the default instance of ``typename XVector::execution_space``.

3. and 4. will throw a runtime exception if ``X.extent(1) != R.extent(0)`` when `X` has rank 2.

Parameters
==========

:space: execution space instance

:R: output vector storing the sum along the column(s) of `X`

:X: input vector(s)


Type Requirements
=================

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2, that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``

- `RV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 or 1, that satisfies

  - ``!std::is_const_v<typename RV::value_type>``
  - ``RV::rank == XMV::rank - 1``

Example
=======

.. literalinclude:: ../../../../example/docs/blas/KokkosBlas1_docs_sum.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   orignal values:
      x: {0, -1, -2, 3, -4}
   result:
      -4
