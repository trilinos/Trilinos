KokkosBlas::swap
################

Defined in header: :code:`KokkosBlas1_swap.hpp`

.. code:: c++

  template <class execution_space, class XVector, class YVector>
  void swap(execution_space const& space, XVector const& x, YVector const& y);

  template <class XVector, class YVector>
  void swap(const XVector& x, const YVector& y);

Exchanges the values of `x` with corresponding values of `y`.

1. iterates over the extents of ``x``, exchange entries of ``x`` and ``y`` on the ``space`` instance
2. iterates over the extents of ``x``, exchange entries of ``x`` and ``y`` on the default instance of ``typename XVector::execution_space``

The function will throw a runtime exception if ``x.extent(0) != y.extent(0) || x.extent(1) != y.extent(1)``

Parameters
==========

:space: execution space instance

:x: vector(s) to swap with ``y``

:y: vector(s) to sawp with ``x``

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XVector` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``std::is_same_v<typename XVector::value_type, typename XVector::non_const_value_type == true``
  - ``Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible == true``

- `YVector` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``std::is_same_v<typename YVector::value_type, typename YVector::non_const_value_type == true``
  - ``Kokkos::SpaceAccessibility<execution_space, typename YVector::memory_space>::accessible == true``
