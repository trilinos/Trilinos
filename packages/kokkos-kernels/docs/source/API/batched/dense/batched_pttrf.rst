KokkosBatched::Pttrf
####################

Defined in header: :code:`KokkosBatched_Pttrf.hpp`

.. code:: c++

    template <typename ArgAlgo>
    struct SerialPttrf {
      template <typename DViewType, typename EViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const DViewType& d,
             const EViewType& e);
    };

Computes a factorization of a tridiagonal matrix A.

1. For a symmetric positive definite tridiagonal matrix A, this computes the :math:`A = L \cdot D \cdot L^T` factorization.
This operation is equivalent to the LAPACK routine ``SPTTRF`` or ``DPTTRF`` for single or double precision.

2. For a complex Hermitian positive definite tridiagonal matrix A, this computes the :math:`A = L \cdot D \cdot L^H` (or :math:`A = U^H \cdot D \cdot U`) factorization.
This operation is equivalent to the LAPACK routine ``CPTTRF`` or ``ZPTTRF`` for single or double precision.

Parameters
==========

:d: Input view containing the diagonal elements of D

:e: Input view containing the subdiagonal (superdiagonal) elements of L (U)

Type Requirements
-----------------

- ``ArgAlgo`` must be ``KokkosBatched::Algo::Pttrf::Unblocked`` for the unblocked algorithm

- ``DViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the diagonal elements (length n)

- ``EViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the superdiagonal or subdiagonal elements (length n-1)

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_pttrs.cpp
  :language: c++

output:

.. code::

   pttrf/pttrs works correctly!