KokkosBatched::Syr2k
####################

Defined in header: :code:`KokkosBatched_Syr2k.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans>
    struct SerialSyr2k {
      template <typename AlphaType, typename AViewType, typename BViewType, typename BetaType, typename CViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const AlphaType alpha, const AViewType &A, const BViewType &B, const BetaType beta,
                                               const CViewType &C);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamSyr2k {
      template <typename AlphaType, typename AViewType, typename BViewType, typename BetaType, typename CViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AlphaType alpha, const AViewType &A,
                                               const BViewType &B, const BetaType beta, const CViewType &C);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamVectorSyr2k {
      template <typename AlphaType, typename AViewType, typename BViewType, typename BetaType, typename CViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AlphaType alpha, const AViewType &A,
                                               const BViewType &B, const BetaType beta, const CViewType &C);
    };

Performs one of the symmetric or hermitian rank 2k operations

.. math::

   \begin{align}
   C &= \alpha A * B^T + \alpha B * A^T + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::NoTranspose)} \\
   C &= \alpha A^T * B + \alpha B^T * A + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   C &= \alpha A * B^H + \bar{\alpha} B * A^H + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::ConjNoTranspose)} \\
   C &= \alpha A^H * B + \bar{\alpha} B^H * A + \beta C \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. If ``ArgTrans == KokkosBatched::Trans::Transpose`` or ``ArgTrans == KokkosBatched::Trans::NoTranspose``, this operation is equivalent to the BLAS routine `SSYR2K <https://www.netlib.org/blas/ssyr2k.f>`_ (`CSYR2K <https://www.netlib.org/blas/csyr2k.f>`_) or `DSYR2K <https://www.netlib.org/blas/dsyr2k.f>`_ (`ZSYR2K <https://www.netlib.org/blas/zsyr2k.f>`_) for single or double precision for real (complex) symmetric matrix.

2. If ``ArgTrans == KokkosBatched::Trans::ConjTranspose`` or ``ArgTrans == KokkosBatched::Trans::ConjNoTranspose``, this operation is equivalent to the BLAS routine `CHER2K <https://www.netlib.org/blas/cher2k.f>`_ or `ZHER2K <https://www.netlib.org/blas/zher2k.f>`_ for single or double precision for complex hermitian matrix.

Parameters
==========

:alpha: Scaling factor
:A: :math:`A` is an n by k matrix if ``ArgTrans == KokkosBatched::Trans::NoTranspose`` or ``ArgTrans == KokkosBatched::Trans::ConjNoTranspose``. Otherwise, :math:`A` is a k by n matrix.
:B: :math:`B` is an n by k matrix if ``ArgTrans == KokkosBatched::Trans::NoTranspose`` or ``ArgTrans == KokkosBatched::Trans::ConjNoTranspose``. Otherwise, :math:`B` is a k by n matrix.
:beta: Scaling factor
:C: :math:`C` is a n by n symmetric (or hermitian) matrix. Only the upper or lower triangular part of :math:`C` is referenced based on ``ArgUplo``.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamSyr2k`` and ``TeamVectorSyr2k``)
- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` to update the upper triangular part of :math:`C`
   - ``KokkosBatched::Uplo::Lower`` to update the lower triangular part of :math:`C`
- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` to perform :math:`C = \alpha A * B^T + \alpha B * A^T + \beta C` on the symmetric matrix :math:`C`
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`C = \alpha A^T * B + \alpha B^T * A + \beta C` on the symmetric matrix :math:`C`
   - ``KokkosBatched::Trans::ConjNoTranspose`` to perform :math:`C = \alpha A * B^H + \bar{\alpha} B * A^H + \beta C` on the hermitian matrix :math:`C`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`C = \alpha A^H * B + \bar{\alpha} B^H * A + \beta C` on the hermitian matrix :math:`C`
- ``ScalarType`` must be a built-in floating point type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``)
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`A`
- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`B`
- ``CViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`C` that satisfies ``std::is_same_v<typename CViewType::value_type, typename CViewType::non_const_value_type>``

.. note::

   If ``ArgTrans == KokkosBatched::Trans::ConjTranspose`` or ``ArgTrans == KokkosBatched::Trans::ConjNoTranspose``, the coefficient :math:`\beta` must be real to ensure Hermitian symmetry of the updated matrix :math:`C`.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_syr2k.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   syr2k works correctly!
