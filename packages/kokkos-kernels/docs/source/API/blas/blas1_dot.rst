KokkosBlas::dot
###############

Defined in header: :code:`KokkosBlas1_dot.hpp`

.. code:: c++

  template <class execution_space, class XVector, class YVector>
  typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
  dot(const execution_space& space, const XVector& x, const YVector& y);

  template <class XVector, class YVector>
  typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
  dot(const XVector& x, const YVector& y);

  template <class execution_space, class RV, class XMV, class YMV>
  void dot(const execution_space& space, const RV& R, const XMV& X, const YMV& Y);

  template <class RV, class XMV, class YMV>
  void dot(const RV& R, const XMV& X, const YMV& Y);

Computes the inner product of vector(s) ``X`` with vector(s) ``Y``.

1. Iterate over the entries of ``X`` and multiply them (or their conjugate) by the corresponding entry in ``Y``, accumulate the results and fence ``space``
2. Iterate over the entries of ``X`` and multiply them (or their conjugate) by the corresponding entry in ``Y``, accumulate the results and fence the default instance of ``typename XVector::execution_space``
3. Iterate over the entries of ``X`` and multiply them (or their conjugate) by the corresponding entry in ``Y``, accumulate the results in ``R`` using the resources specified by ``space``
4. Iterate over the entries of ``X`` and multiply them (or their conjugate) by the corresponding entry in ``Y``, accumulate the results in ``R`` using the resources specified by the default instance of ``typename XVector::execution_space``

The function will throw a runtime exception if any of the following conditions are **not** met
  - ``X.extent(0) == Y.extent(0)``
  - ``X.extent(1) = 1 || Y.extent(1) || X.extent(1) == Y.extent(1)``
  - ``R.extent(0) = max(X.extent(1), Y.extent(1))``

The result (returned or stored value) is undefined if ``X`` has no entries.

..
   .. note::

     We should simplify our API by using unique template parameter names for all the overloads, similarly to the actual input parameter names. The input parameter names should use unique capitalization across overloads.

     The current overloads do not really allow for the interface of `dsdot <https://netlib.org/lapack/explore-html/d1/dcc/group__dot_ga17a1bc70455f422325f92943b48c7240.html#ga17a1bc70455f422325f92943b48c7240>`_ and `sdsdot <https://netlib.org/lapack/explore-html/d1/dcc/group__dot_ga5e09e98ca27006a197d7c5fa49a9da4b.html#ga5e09e98ca27006a197d7c5fa49a9da4b>`_, it also does not distinguish between `cdotc` and `cdotu`.

Parameters
==========

:space: execution space instance that specifies the resources and stream/queue used to execute this kernel

:R: computed inner products

:X, Y: vectors used to compute the inner products

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XVector` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible == true``

- `YVector` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename YVector::memory_space>::accessible == true``

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible == true``

- `YMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible == true``

- `RV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 or 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename RV::memory_space>::accessible == true``
  - ``std::is_same_v<typename RV::value_type, typename RV::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas1_wiki_dot.cpp
  :language: c++

output:

.. code::

  X_dot_Y: 600 Expected: 600 Diff: 0
