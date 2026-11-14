KokkosBlas::mult
################

Defined in header: :code:`KokkosBlas1_mult.hpp`

.. code:: c++

  template <class execution_space, class YMV, class AV, class XMV>
  void mult(const execution_space& space, typename YMV::const_value_type& beta, const YMV& Y,
            typename AV::const_value_type& alpha, const AV& A, const XMV& X);

  template <class YMV, class AV, class XMV>
  void mult(typename YMV::const_value_type& beta, const YMV& Y, typename AV::const_value_type& alpha, const AV& A,
            const XMV& X);

computes an element wise multiplication of `A` and `X` scaled by `alpha`, then accumulate into `Y` scaled by `beta`.

.. math::

   Y = \beta * Y + \alpha * A * X

1. Perform the entry wise multipliaction and addition using ``space`` instance.
2. Perform the entry wise multipliaction and addition using the default instance of ``typename XMV::execution_space``.

The function will throw a runtime exception if
1. the vectors do not have same length: ``Y.extent(0) != A.extent(0) || Y.extent(0) != X.extent(0)``,
2. the number of vectors in `X` and `Y` does not match: ``Y.extent(2) != X.extent(1)`` for the rank-2 case.

Parameters
==========

:space: execution space instance

:beta: scaling factor for `Y`

:Y: output vector(s)

:alpha: scaling factor for the `A`, `X` product

:A: input vector to multiply element wise with `X`

:X: input vector(s)


Type Requirements
=================

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `YMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2, that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``
  - ``std::is_same_v<typename YMV::value_type, typename YMV::non_const_value_type>``

- `AV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2, that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``
  - ``XMV::rank == YMV::rank``

Example
=======

.. literalinclude:: ../../../../example/docs/blas/KokkosBlas1_docs_mult.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   orignal values:
      beta: 2.0,  alpha: -1.0
      y: {3, 2, 3, -2, 4}
      a: {-3, -1, 2, -2, -3}
      x: {0, -1, -2, 3, -4}
   final values:
      y: {6, 3, 10, 2, -4}
