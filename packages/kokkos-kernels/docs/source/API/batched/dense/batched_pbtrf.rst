KokkosBatched::Pbtrf
####################

Defined in header: :code:`KokkosBatched_Pbtrf.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgAlgo>
    struct SerialPbtrf {
      template <typename ABViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ABViewType &ab);
    };

Computes the Cholesky factorization of a real symmetric (or complex Hermitian) positive definite band matrix A.

1. For a real symmetric positive definite band matrix A, this computes the :math:`A = U^T \cdot U \: \text{(if ArgUplo == Uplo::Upper)}` or :math:`A = L \cdot L^T \: \text{(if ArgUplo == Uplo::Lower)}` factorization,
where U is an upper triangular matrix and L is a lower triangular matrix. This operation is equivalent to the LAPACK routine ``SPBTRF`` or ``DPBTRF`` for single or double precision.

2. For a complex Hermitian positive definite band matrix A, this computes the :math:`A = U^H \cdot U \: \text{(if ArgUplo == Uplo::Upper)}` or :math:`A = L \cdot L^H \: \text{(if ArgUplo == Uplo::Lower)}` factorization,
where U is an upper triangular matrix and L is a lower triangular matrix. This operation is equivalent to the LAPACK routine ``CPBTRF`` or ``ZPBTRF`` for single or double precision.

Parameters
==========

:ab: On input, the upper or lower triangle of the symmetric band matrix A. On output, the triangular factor U or L from the Cholesky factorization. See `LAPACK reference <https://www.netlib.org/lapack/lug/node124.html>`_ for the band storage format.

Type Requirements
-----------------

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` for upper triangular factorization
   - ``KokkosBatched::Uplo::Lower`` for lower triangular factorization

- ``ArgAlgo`` must be ``KokkosBatched::Algo::Pbtrf::Unblocked`` for the unblocked algorithm
- ``ABViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the band matrix A

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_pbtrs.cpp
  :language: c++

output:

.. code::

   pbtrf/pbtrs works correctly!
