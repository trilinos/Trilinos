KokkosBlas::nrm1
################

Defined in header: :code:`KokkosBlas1_nrm1.hpp`

.. code:: c++

  template <class execution_space, class XVector>
  typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type
  nrm1(const execution_space& space, const XVector& x);

  template <class XVector>
  typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type
  nrm1(const XVector& x)

  template <class execution_space, class RV, class XMV>
  void nrm1(const execution_space& space, const RV& R, const XMV& X)

  template <class RV, class XMV>
  void nrm1(const RV& R, const XMV& X);

Computes the norm 1 (sum of absolute values) of the X vector(s) and returns or stores the result in R. Also known as ``asum``.

1. accumulates the sum of absolute values of ``X`` and fences the ``space`` instance
2. accumulates the sum of absolute values of ``X`` and fences the default instance of ``typename XVector::execution_space``
3. accumulates the sum of absolute values of each column of ``X``, stores it in the corresponding entry of ``R`` and fences the ``space`` instance
4. accumulates the sum of absolute values of each column of ``X``, stores it in the corresponding entry of ``R`` and fences the default instance of ``typename XVector::execution_space``

The result (returned or stored value) is undefined if ``X`` has no entries.

Parameters
==========

:space: execution space instance

:X: vector(s) to compute the norm1

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
  - ``std::is_same_v<typename RV::value_type, typename Kokkos::Details::InnerProductSpaceTraits<typename XMV::non_const_value_type>::mag_type> == true``

Example
=======

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas1_wiki_nrm1.cpp
  :language: c++

output:

.. code::

   X_nrm: 300 Expected: 300
