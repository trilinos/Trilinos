KokkosBatched::Pttrs
####################

Defined in header: :code:`KokkosBatched_Pttrs.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgAlgo>
    struct SerialPttrs {
      template <typename DViewType, typename EViewType, typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const DViewType& d,
             const EViewType& e,
             const BViewType& b);
    };

Solves a tridiagonal system of the form :math:`A \cdot X = B` using the factorization computed by ``Pttrf``.

1. For a symmetric positive definite tridiagonal matrix A, this solves the linear system of equation :math:`A \cdot X = B` for X using the :math:`A = L \cdot D \cdot L^T` factorization computed by ``Pttrf``.
This operation is equivalent to the LAPACK routine ``SPTTRS`` or ``DPTTRS`` for single or double precision.

2. For a complex Hermitian positive definite tridiagonal matrix A, this solves :math:`A \cdot X = B` for X using the :math:`A = L \cdot D \cdot L^H` or (:math:`A = U^H \cdot D \cdot U`) factorization computed by ``Pttrf``.
This operation is equivalent to the LAPACK routine ``CPTTRS`` or ``ZPTTRS`` for single or double precision.

Parameters
==========

:d: Input view containing the diagonal elements of D from the factorization

:e: Input view containing the subdiagonal (superdiagonal) elements of L (U) from the factorization

:b: Input/output view containing the right-hand side on input and the solution on output

Type Requirements
-----------------

- ``ArgUplo`` must be one of the following:

    - ``KokkosBatched::Uplo::Upper`` if vector `e` specifies the superdiagonal of a unit bidiagonal matrix U
    - ``KokkosBatched::Uplo::Lower`` if vector `e` specifies the subdiagonal of a unit bidiagonal matrix L
    - This parameter is unused for real matrices

- ``ArgAlgo`` must be ``KokkosBatched::Algo::Pttrs::Unblocked`` for the unblocked algorithm

- ``DViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the diagonal elements (length n)

- ``EViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the superdiagonal or subdiagonal elements (length n-1)

- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the right-hand side that satisfies

  - ``std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_pttrs.cpp
  :language: c++

output:

.. code::

   pttrf/pttrs works correctly!