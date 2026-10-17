KokkosBlas::scal
################

Defined in header: :code:`KokkosBlas1_scal.hpp`

.. code:: c++

  template <class execution_space, class RMV, class AV, class XMV>
  void scal(const execution_space& space, const RMV& R, const AV& a, const XMV& X);

  template <class RMV, class AV, class XMV>
  void scal(const RMV& R, const AV& a, const XMV& X);

Scales the vector(s) in `X` by `a` and stores the result(s) in `R`.

1. store X scaled by a in R and fences the ``space`` instance
2. store X scaled by a in R and fences the default instance of ``typename RMV::execution_space``

The function will throw a runtime exception if ``X.extent(0) != R.extent(0) || X.extent(1) != R.extent(1)``

Parameters
==========

:space: execution space instance

:R: scaled version of the vector(s) of ``X``

:a: scaling factor(s) for the vector(s) of ``X``

:X: vector(s) to scale

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible == true``
  - ``RMV::rank == XMV::rank``

- `RMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename RMV::memory_space>::accessible == true``
  - ``RMV::rank == XMV::rank``

- `AV` must be one of:

  - A **scalar** value (e.g. ``double``). A single value scales all elements of ``X``.

  - A **rank-0** Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ (e.g. ``Kokkos::View<const double>``).

  - A **rank-1** Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ with ``AV::rank == 1``. Only valid when ``RMV::rank == XMV::rank == 2`` (multivector). Column ``j`` of ``X`` is scaled by ``a(j)``.

Example
=======

.. literalinclude:: ../../../../example/docs/blas/KokkosBlas1_docs_scal.cpp
  :language: c++

output:

.. code::

   orignal values:
      {-3, -3, -3, -3, -3}
   final values:
      {6, 6, 6, 6, 6}
