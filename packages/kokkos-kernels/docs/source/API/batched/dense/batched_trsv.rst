KokkosBatched::Trsv
###################

Defined in header: :code:`KokkosBatched_Trsv_Decl.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
    struct SerialTrsv {
      template <typename ScalarType, typename AViewType, typename bViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A,
                                               const bViewType &b);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
    struct TeamTrsv {
      template <typename ScalarType, typename AViewType, typename bViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha,
                                               const AViewType &A, const bViewType &b);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
    struct TeamVectorTrsv {
      template <typename ScalarType, typename AViewType, typename bViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha,
                                               const AViewType &A, const bViewType &b);
    };


Solves a system of the linear equations :math:`A \cdot X = \alpha B` or :math:`A^T \cdot X = \alpha B` or :math:`A^H \cdot X = \alpha B` where :math:`A` is an n-by-n unit or non-unit, upper or lower triangular matrix.

1. For a real matrix :math:`A`, this solves a system of the linear equations :math:`A \cdot X = \alpha B` or :math:`A^T \cdot X = \alpha B`.
   This operation is equivalent to the BLAS routine ``STRSV`` or ``DTRSV`` for single or double precision.

2. For a complex matrix :math:`A`, this solves a system of the linear equations :math:`A \cdot X = \alpha B` or :math:`A^T \cdot X = \alpha B` or :math:`A^H \cdot X = \alpha B`.
   This operation is equivalent to the BLAS routine ``CTRSV`` or ``ZTRSV`` for single or double precision.

.. note::

  No test for singularity or near-singularity is included in this routine. Such tests must be performed before calling this routine.

Parameters
==========

:member: Kokkos team member handle (only for ``TeamTrsv`` and ``TeamVectorTrsv``).
:alpha: Scalar multiplier for :math:`b` (usually set to one).
:A: Input view containing the upper or lower triangular matrix.
:b: Input/output view containing the right-hand side on input and the solution on output.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamTrsv`` and ``TeamVectorTrsv``).

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` for upper triangular solve
   - ``KokkosBatched::Uplo::Lower`` for lower triangular solve

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` to solve a system :math:`A \cdot X = \alpha B`
   - ``KokkosBatched::Trans::Transpose`` to solve a system :math:`A^T \cdot X = \alpha B`
   - ``KokkosBatched::Trans::ConjTranspose`` to solve a system :math:`A^H \cdot X = \alpha B`

- ``ArgDiag`` must be one of the following:
   - ``KokkosBatched::Diag::Unit`` for the unit triangular matrix :math:`A`
   - ``KokkosBatched::Diag::NonUnit`` for the non-unit triangular matrix :math:`A`

- ``ArgAlgo`` must be one of the following:
   - ``KokkosBatched::Algo::tbsv::Blocked`` for the blocked algorithm
   - ``KokkosBatched::Algo::tbsv::Unblocked`` for the unblocked algorithm

- ``ScalarType`` must be a built-in arithmetic type like ``float``, ``double``, ``Kokkos::complex<float>``, or ``Kokkos::complex<double>``.
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the band matrix A
- ``bViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the right-hand side that satisfies
  - ``std::is_same_v<typename bViewType::value_type, typename bViewType::non_const_value_type> == true``

.. note::

  Some combinations of template parameters may not be supported yet.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_trsv.cpp
  :language: c++

output:

.. code::

   trsv works correctly!
