KokkosBlas::nrminf
##################

Defined in header: :code:`KokkosBlas1_nrminf.hpp`

.. code:: c++

  template <class execution_space, class XVector, typename std::enable_if<Kokkos::is_execution_space_v<execution_space>, int>::type = 0>
  typename KokkosKernels::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type
  nrminf(const execution_space& space, const XVector& x);

  template <class XVector>
  typename KokkosKernels::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type
  nrminf(const XVector& x);

  template <class execution_space, class RV, class XMV>
  void nrminf(const execution_space& space, const RV& R, const XMV& X,
              typename std::enable_if<Kokkos::is_view_v<RV>, int>::type = 0);

  template <class RV, class XMV>
  void nrminf(const RV& R, const XMV& X, typename std::enable_if<Kokkos::is_view_v<RV>, int>::type = 0)

Computes the norm inf (maximum absolute value) of the X vector(s) and returns or stores the result in R.

1. returns the largest magnitude in ``X`` and fences the ``space`` instance.
2. returns the largest magnitude in ``X`` and fences the default instance of ``typename XVector::execution_space``.
3. returns the largest magnitude in each column of ``X``, stores it in the corresponding entry of ``R`` and fences the ``space`` instance.
4. returns the largest magnitude in each column of ``X``, stores it in the corresponding entry of ``R`` and fences the default instance of ``typename XVector::execution_space``.

The result (returned or stored value) is undefined if ``X`` has no entries.

Note: For 3. and 4. the return view is expected to be templated on `Kokkos::HostSpace`, the function will fail with a segfault at runtime otherwise.

Parameters
==========

:space: execution space instance

:X: vector(s) to compute the norminf

:R: computed norm(s)

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XVector` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible == true``

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible == true``

- `RV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``RV::rank == XMV::rank - 1``
  - ``std::is_same_v<typename RV::value_type, typename RV::non_const_value_type> == true``
  - ``std::is_same_v<typename RV::value_type, typename KokkosKernels::Details::InnerProductSpaceTraits<typename XMV::non_const_value_type>::mag_type> == true``
  - ``RV::memory_space == Kokkos::HostSpace``

Example
=======

.. literalinclude:: ../../../../example/docs/blas/KokkosBlas1_docs_nrminf.cpp
  :language: c++

output:

.. code::

   X_nrm: 50.011 Expected: 50.011
