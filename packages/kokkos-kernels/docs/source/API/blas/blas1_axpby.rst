KokkosBlas::axpby
#################

Defined in header: :code:`KokkosBlas1_axpby.hpp`

.. code:: c++

  template <class execution_space, class AV, class XMV, class BV, class YMV>
  void axpby(const execution_space& exec_space, const AV& a, const XMV& X, const BV& b, const YMV& Y) {

  template <class AV, class XMV, class BV, class YMV>
  void axpby(const AV& a, const XMV& X, const BV& b, const YMV& Y) {

Scale ``Y`` by coefficient ``b`` and then add entries of ``X`` scaled by coefficient ``a``: ``Y := a*X + b*Y``

1. Compute ``Y := a*X + b*Y`` on the provided ``exec_space``.
2. Compute ``Y := a*X + b*Y`` on the default execution space instance of type ``XMV::execution_space``.

``a`` and ``b`` may each be any of the following:
  - a scalar value
  - a rank-0 host-accessible ``Kokkos::View``
  - a rank-0 device-accessible ``Kokkos::View``
  - a rank-1 device-accessible ``Kokkos::View`` with extent 1 (the coefficient will be applied to all columns)
  - a rank-1 device-accessible ``Kokkos::View`` with extent ``Y.extent(1)`` (one coefficient per column)

The function will throw a runtime exception if any of the following conditions are **not** met:
  - ``Y.extent(0) == X.extent(0) && Y.extent(1) == X.extent(1)``

If ``a`` is a rank-1 ``View``:
  - ``a.extent(0) == 1 || a.extent(0) == X.extent(1)``

If ``b`` is a rank-1 ``View``:
  - ``b.extent(0) == 1 || b.extent(0) == Y.extent(1)``

Note: if ``b`` is a single coefficient to apply to all columns of ``Y``, and ``b`` is zero,
then the values of ``Y`` on input are treated as zeros. This means that the ``axpby`` will
overwrite any NaN or infinity values in ``Y``.

Parameters
==========

:space: execution space instance

:a: scaling factor(s) applied to ``X``

:X: input vector or multi-vector

:b: scaling factor(s) applied to ``Y``

:Y: output vector or multi-vector

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``

- `YMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``YMV::rank == XMV::rank``
  - ``!std::is_const_v<typename YMV::value_type>``
  - ``Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible``

- `AV` must be a scalar or a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies one of the following:

  - ``AV::rank == 0`` 
  - ``AV::rank == 1 && Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible``

- `BV` must be a scalar or a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies one of the following:

  - ``BV::rank == 0`` 
  - ``BV::rank == 1 && Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible``

