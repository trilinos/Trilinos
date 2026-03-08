KokkosBatched::Trsm
###################

Defined in header: :code:`KokkosBatched_Trsm_Decl.hpp`

.. code:: c++

    template <typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
    struct SerialTrsm {
      template <typename ScalarType, typename AViewType, typename BViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B);
    };

    template <typename MemberType, typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag,
              typename ArgAlgo>
    struct TeamTrsm {
      template <typename ScalarType, typename AViewType, typename BViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                               const BViewType &B);
    };
    
    template <typename MemberType, typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag,
              typename ArgAlgo>
    struct TeamVectorTrsm {
      template <typename ScalarType, typename AViewType, typename BViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                               const BViewType &B);
    };


Solves a system of the linear equations :math:`op(A) \cdot X = \alpha B` or :math:`X \cdot op(A) = \alpha B` where :math:`\alpha` is a scalar, :math:`X` and :math:`B` are m-by-n matrices, :math:`A` is a unit or non-unit, upper or lower triangular matrix and 
:math:`op(A)` is one of :math:`A`, :math:`A^T`, or :math:`A^H`. The matrix :math:`X` is overwritten on :math:`B`.

1. For a real matrix :math:`A`, :math:`op(A)` is one of :math:`A` or :math:`A^T`.
   This operation is equivalent to the BLAS routine ``STRSM`` or ``DTRSM`` for single or double precision.

2. For a complex matrix :math:`A`, :math:`op(A)` is one of :math:`A`, :math:`A^T`, or :math:`A^H`
   This operation is equivalent to the BLAS routine ``CTRSM`` or ``ZTRSM`` for single or double precision.

Parameters
==========

:member: Kokkos team member handle (only for ``TeamTrsm`` and ``TeamVectorTrsm``).
:alpha: Scalar multiplier for :math:`B`.
:A: Input view containing the upper or lower triangular matrix.
:B: Input/output view containing the right-hand side on input and the solution on output.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamTrsm`` and ``TeamVectorTrsm``).

- ``ArgSide`` must be one of the following:
   - ``KokkosBatched::Side::Left`` to solve a system :math:`op(A) \cdot X = \alpha B`
   - ``KokkosBatched::Side::Right`` to solve a system :math:`X \cdot op(A) = \alpha B`

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` for upper triangular solve
   - ``KokkosBatched::Uplo::Lower`` for lower triangular solve

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` for :math:`op(A) = A`
   - ``KokkosBatched::Trans::Transpose`` for :math:`op(A) = A^T`
   - ``KokkosBatched::Trans::ConjTranspose`` for :math:`op(A) = A^H`

- ``ArgDiag`` must be one of the following:
   - ``KokkosBatched::Diag::Unit`` for the unit triangular matrix :math:`A`
   - ``KokkosBatched::Diag::NonUnit`` for the non-unit triangular matrix :math:`A`

- ``ArgAlgo`` must be one of the following:
   - ``KokkosBatched::Algo::trsm::Blocked`` for the blocked algorithm
   - ``KokkosBatched::Algo::trsm::Unblocked`` for the unblocked algorithm

- ``ScalarType`` must be a built-in arithmetic type like ``float``, ``double``, ``Kokkos::complex<float>``, or ``Kokkos::complex<double>``.
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the band matrix A
- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the right-hand side that satisfies
  - ``std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type> == true``

.. note::

  Some combinations of template parameters may not be supported yet.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_trsm.cpp
  :language: c++

output:

.. code::

   trsm works correctly!
