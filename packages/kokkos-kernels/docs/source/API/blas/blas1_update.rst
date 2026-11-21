KokkosBlas::update
##################

Defined in header: :code:`KokkosBlas1_update.hpp`

.. code:: c++

  template <class execution_space, class XMV, class YMV, class ZMV>
  void update(const execution_space& space, const typename XMV::non_const_value_type& alpha, const XMV& X,
              const typename YMV::non_const_value_type& beta, const YMV& Y,
              const typename ZMV::non_const_value_type& gamma, const ZMV& Z);

  template <class XMV, class YMV, class ZMV>
  void update(const typename XMV::non_const_value_type& alpha, const XMV& X,
              const typename YMV::non_const_value_type& beta, const YMV& Y,
              const typename ZMV::non_const_value_type& gamma, const ZMV& Z);

updates the values of `Z` using the scales sum of the values of the vector(s) `X`, `Y` and `Z` itself.

.. math::

   Z = \alpha X + \beta Y + \gamma Z

1. Adds the vectors `X`, `Y` and `Z`, scaled by :math:`\alpha`, :math:`\beta` and :math:`\gamma` resp., together into `Z`.
2. Computes the sum of the entries in `X` using the default instance of ``typename ZMV::execution_space`` and returns them.

The functions will throw at runtime unless the first and second extents of `X`, `Y` and `Z` are equal.

Parameters
==========

:space: execution space instance

:alpha: scaling factor for the values of `X`

:X: input vector(s)

:beta: scaling factor for the values of `Y`

:Y: input vector(s)

:gamma: scaling factor for the values of `Z`

:Z: output vector(s)


Type Requirements
=================

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2, that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``

- `YMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2, that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``

- `ZMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2, that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``
  - ``!std::is_cont_v<typename ZMV::value_type>``
  - ``ZMV::rank == XMV::rank && ZMV::rank == YMV::rank``

Example
=======

.. literalinclude:: ../../../../example/docs/blas/KokkosBlas1_docs_update.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   orignal values:
     alpha: -1.0,  beta: 2.0,  gamma:
     x: {0, -1, -2, 3, -4}
     y: {3, 2, 3, -2, 4}
     z: {-3, -1, 2, -2, -3}

   final values:
     z: {3, 4, 10, -9, 9}
