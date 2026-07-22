KokkosBatched::Tbsv
###################

Defined in header: :code:`KokkosBatched_Tbsv.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
    struct SerialTbsv {
      template <typename AViewType, typename XViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &X, const int k);
    };


Solves a system of the linear equations :math:`A \cdot X = B` or :math:`A^T \cdot X = B` or :math:`A^H \cdot X = B` where :math:`A` is an n-by-n unit or non-unit, upper or lower triangular band matrix with :math:`(k + 1)` diagonals.

1. For a real band matrix :math:`A`, this solves a system of the linear equations :math:`A \cdot X = B` or :math:`A^T \cdot X = B`.
   This operation is equivalent to the BLAS routine ``STBSV`` or ``DTBSV`` for single or double precision.

2. For a complex band matrix :math:`A`, this solves a system of the linear equations :math:`A \cdot X = B` or :math:`A^T \cdot X = B` or :math:`A^H \cdot X = B`.
   This operation is equivalent to the BLAS routine ``CTBSV`` or ``ZTBSV`` for single or double precision.

.. note::

  No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.

Parameters
==========

:A: Input view containing the upper or lower triangular band matrix. See `LAPACK reference <https://www.netlib.org/lapack/lug/node124.html>`_ for the band storage format.
:X: Input/output view containing the right-hand side on input and the solution on output.
:k: The number of superdiagonals or subdiagonals within the band of :math:`A`. :math:`k >= 0`


Type Requirements
-----------------

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` for upper triangular solve
   - ``KokkosBatched::Uplo::Lower`` for lower triangular solve

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` to solve a system :math:`A \cdot X = B`
   - ``KokkosBatched::Trans::Transpose`` to solve a system :math:`A^T \cdot X = B`
   - ``KokkosBatched::Trans::ConjTranspose`` to solve a system :math:`A^H \cdot X = B`

- ``ArgDiag`` must be one of the following:
   - ``KokkosBatched::Diag::Unit`` for the unit triangular matrix :math:`A`
   - ``KokkosBatched::Diag::NonUnit`` for the non-unit triangular matrix :math:`A`

- ``ArgAlgo`` must be ``KokkosBatched::Algo::tbsv::Unblocked`` for the unblocked algorithm
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the band matrix A
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the right-hand side that satisfies
  - ``std::is_same_v<typename XViewType::value_type, typename XViewType::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_tbsv.cpp
  :language: c++

output:

.. code::

   tbsv works correctly!
