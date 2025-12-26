KokkosBatched::Pbtrs
####################

Defined in header: :code:`KokkosBatched_Pbtrs.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgAlgo>
    struct SerialPbtrs {
      template <typename ABViewType, typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ABViewType &ab,
             const BViewType& b);
    };

Solves a system of the linear equations :math:`A \cdot X = B` with a symmetric (or Hermitian) positive definite band matrix A using the Cholesky factorization computed by ``Pbtrf``.

1. For a real symmetric positive definite band matrix A, this solves a system of the linear equations :math:`A \cdot X = B` for X using the :math:`A = U^T \cdot U \: \text{(if ArgUplo == Uplo::Upper)}` or :math:`A = L \cdot L^T \: \text{(if ArgUplo == Uplo::Lower)}` factorization 
computed by ``Pbtrf``. Here, U is an upper triangular matrix and L is a lower triangular matrix. This operation is equivalent to the LAPACK routine ``SPBTRS`` or ``DPBTRS`` for single or double precision.

2. For a complex Hermitian positive definite band matrix A, this solves a system of the linear equations :math:`A \cdot X = B` for X using the :math:`A = U^H \cdot U \: \text{(if ArgUplo == Uplo::Upper)}` or :math:`A = L \cdot L^H \: \text{(if ArgUplo == Uplo::Lower)}` factorization 
computed by ``Pbtrs``. Here, U is an upper triangular matrix and L is a lower triangular matrix. This operation is equivalent to the LAPACK routine ``CPBTRS`` or ``ZPBTRS`` for single or double precision.

Parameters
==========

:ab: Input view containing the triangular factor U or L from the Cholesky factorization. See `LAPACK reference <https://www.netlib.org/lapack/lug/node124.html>`_ for the band storage format.

:b: Input/output view containing the right-hand side on input and the solution on output.


Type Requirements
-----------------

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` for upper triangular factorization
   - ``KokkosBatched::Uplo::Lower`` for lower triangular factorization

- ``ArgAlgo`` must be ``KokkosBatched::Algo::Pbtrs::Unblocked`` for the unblocked algorithm
- ``ABViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the band matrix A
- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the right-hand side that satisfies
  - ``std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_pbtrs.cpp
  :language: c++

output:

.. code::

   pbtrf/pbtrs works correctly!
