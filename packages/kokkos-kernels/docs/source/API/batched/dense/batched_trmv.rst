KokkosBatched::Trmv
###################

Defined in header: :code:`KokkosBatched_Trmv.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans, typename ArgDiag>
    struct SerialTrmv {
      template <typename ScalarType, typename AViewType, typename XViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const XViewType &x);
    };

Performs the matrix-vector operations :math:`A` by vector :math:`x` with scaling factor :math:`alpha`

.. math::

   \begin{align}
   x &= \alpha A * x \: \text{(if ArgTrans == KokkosBatched::Trans::NoTranspose)} \\
   x &= \alpha A^T * x \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   x &= \alpha A^H * x \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. This operation is equivalent to the BLAS routine `STRMV <https://www.netlib.org/blas/strmv.f>`_ (`CTRMV <https://www.netlib.org/blas/ctrmv.f>`_) or `DTRMV <https://www.netlib.org/blas/dtrmv.f>`_ (`ZTRMV <https://www.netlib.org/blas/ztrmv.f>`_) for single or double precision for real (complex) triangular matrix.

Parameters
==========

:alpha: Scalar multiplier for :math:`x` (usually set to one).
:A: :math:`A` is an n by n upper or lower triangular matrix.
:x: On entry, :math:`x` is a length n vector. On exit, :math:`x` contains the result of the matrix-vector operation.

Type Requirements
-----------------

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` for upper triangular solve
   - ``KokkosBatched::Uplo::Lower`` for lower triangular solve

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` to perform :math:`x = \alpha A * x`
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`x = \alpha A^T * x`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`x = \alpha A^H * x`

- ``ArgDiag`` must be one of the following:
   - ``KokkosBatched::Diag::Unit`` for the unit triangular matrix :math:`A`
   - ``KokkosBatched::Diag::NonUnit`` for the non-unit triangular matrix :math:`A`

- ``ScalarType`` must be a built-in arithmetic type like ``float``, ``double``, ``Kokkos::complex<float>``, or ``Kokkos::complex<double>``.
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the triangular matrix A
- ``xViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the right-hand side that satisfies
  - ``std::is_same_v<typename xViewType::value_type, typename xViewType::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_trmv.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   trmv works correctly!
