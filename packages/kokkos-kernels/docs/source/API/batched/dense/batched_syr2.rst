KokkosBatched::Syr2
###################

Defined in header: :code:`KokkosBatched_Syr2.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans>
    struct SerialSyr2 {
      template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                               const AViewType &a);
    };

Perform a symmetric rank-2 update of matrix :math:`A` by vectors :math:`x` and :math:`y` with scaling factor :math:`alpha`

.. math::

   \begin{align}
   A &= A + \alpha (x \cdot y^T) + \alpha (y \cdot x^T) \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   A &= A + \alpha (x \cdot y^H) + \bar{\alpha} (y \cdot x^H) \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. If ``ArgTrans == KokkosBatched::Trans::Transpose``, this operation is equivalent to the BLAS routine ``SSYR2`` (``CSYR2``) or ``DSYR2`` (``ZSYR2``) for single or double precision for real (complex) matrix.
2. If ``ArgTrans == KokkosBatched::Trans::ConjTranspose``, this operation is equivalent to the BLAS routine ``CHER2`` or ``ZHER2`` for single or double precision for complex matrix.

Parameters
==========

:alpha: Scaling factor of the rank-2 update
:x: On input, :math:`x` is a length n vector.
:y: On input, :math:`y` is a length n vector.
:A: With ``ArgUplo == KokkosBatched::Uplo::Upper``, the upper triangular part of :math:`A` is referenced and updated. With ``ArgUplo == KokkosBatched::Uplo::Lower``, the lower triangular part of :math:`A` is referenced and updated.

Type Requirements
-----------------

- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` to update the upper triangular part of :math:`A`
   - ``KokkosBatched::Uplo::Lower`` to update the lower triangular part of :math:`A`
- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`A = A + \alpha (x \cdot y^T) + \alpha (y \cdot x^T)`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`A = A + \alpha (x \cdot y^H) + \bar{\alpha} (y \cdot x^H)`
- ``ScalarType`` must be a built-in floating point type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``)
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`X`
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`Y`
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`A` that satisfies ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_syr2.cpp
  :language: c++

output:

.. code::

   syr2 works correctly!
