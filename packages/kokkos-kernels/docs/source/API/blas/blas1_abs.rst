KokkosBlas::abs
###############

Defined in header: :code:`KokkosBlas1_abs.hpp`

.. code:: c++

  template <class execution_space, class RMV, class XMV>
  void abs(const execution_space& space, const RMV& R, const XMV& X);

  template <class RMV, class XMV>
  void abs(const RMV& R, const XMV& X);

Replaces the entries of `R` by the absolute value of the corresponding entries in `X`.

1. Overwrites the entries of R with the magnitude of the corresponding entries in X and fences the ``space`` instance.
2. Overwrites the entries of R with the magnitude of the corresponding entries in X and fences the default instance of ``typename RMV::execution_space``.

The function will throw a runtime exception if ``R.extent(0) != X.extent(0) || R.extent(1) != X.extent(1)``

Parameters
==========

:space: execution space instance

:R: output vector(s) storing the magnitude of the entries of ``X``

:X: input vector(s)

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `RMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``
  - ``std::is_same_v<typename RMV::value_type, typename RMV::non_const_value_type> == true``

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``
  - ``XMV::rank == RMV::rank``
