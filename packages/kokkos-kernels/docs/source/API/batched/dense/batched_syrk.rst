KokkosBatched::Syrk
###################

Defined in header: :code:`KokkosBatched_Syrk.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans>
    struct SerialSyrk {
      template <typename ScalarType, typename AViewType, typename CViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const ScalarType beta,
                                               const CViewType &C);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamSyrk {
      template <typename ScalarType, typename AViewType, typename CViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                               const ScalarType beta, const CViewType &C);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamVectorSyrk {
      template <typename ScalarType, typename AViewType, typename CViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                               const ScalarType beta, const CViewType &C);
    };

Performs one of the symmetric or hermitian rank-k update

.. math::

   \begin{align}
   C &= \alpha A * A^T + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::NoTranspose)} \\
   C &= \alpha A^T * A + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   C &= \alpha A * A^H + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::ConjNoTranspose)} \\
   C &= \alpha A^H * A + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. If ``ArgTrans == KokkosBatched::Trans::Transpose`` or ``ArgTrans == KokkosBatched::Trans::NoTranspose``, this operation is equivalent to the BLAS routine `SSYRK <https://www.netlib.org/blas/ssyrk.f>`_ (`CSYRK <https://www.netlib.org/blas/csyrk.f>`_) or `DSYRK <https://www.netlib.org/blas/dsyrk.f>`_ (`ZSYRK <https://www.netlib.org/blas/zsyrk.f>`_) for single or double precision for real (complex) symmetric matrix.

2. If ``ArgTrans == KokkosBatched::Trans::ConjTranspose`` or ``ArgTrans == KokkosBatched::Trans::ConjNoTranspose``, this operation is equivalent to the BLAS routine `CHERK <https://www.netlib.org/blas/cherk.f>`_ or `ZHERK <https://www.netlib.org/blas/zherk.f>`_ for single or double precision for complex hermitian matrix.

Parameters
==========

:alpha: Scaling factor
:A: :math:`A` is an n by k matrix if ``ArgTrans == KokkosBatched::Trans::NoTranspose`` or ``ArgTrans == KokkosBatched::Trans::ConjNoTranspose``. Otherwise, :math:`A` is a k by n matrix.
:C: :math:`C` is a n by n symmetric (or hermitian) matrix. Only the upper or lower triangular part of :math:`C` is referenced based on ``ArgUplo``.
:beta: Scaling factor

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamSyrk`` and ``TeamVectorSyrk``)
- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` to update the upper triangular part of :math:`C`
   - ``KokkosBatched::Uplo::Lower`` to update the lower triangular part of :math:`C`
- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` to perform :math:`C = \alpha A * A^T + \beta C` on the symmetric matrix :math:`C`
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`C = \alpha A^T * A + \beta C` on the symmetric matrix :math:`C`
   - ``KokkosBatched::Trans::ConjNoTranspose`` to perform :math:`C = \alpha A * A^H + \beta C` on the hermitian matrix :math:`C`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`C = \alpha A^H * A + \beta C` on the hermitian matrix :math:`C`
- ``ScalarType`` must be a built-in floating point type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``)
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`A`
- ``CViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`C` that satisfies ``std::is_same_v<typename CViewType::value_type, typename CViewType::non_const_value_type>``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_syrk.cpp
  :language: c++

output:

.. code::

   syrk works correctly!
