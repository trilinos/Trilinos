KokkosBatched::Ger
##################

Defined in header: :code:`KokkosBatched_Ger.hpp`

.. code:: c++

    template <typename ArgTrans>
    struct SerialGer {
      template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                               const AViewType &a);
    };

Perform a rank-1 update of matrix :math:`A` by vectors :math:`x` and :math:`y` with scaling factor :math:`alpha`

.. math::

   \begin{align}
   A &= A + \alpha (x * y^T) \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   A &= A + \alpha (x * y^H) \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. For real vectors :math:`x` and :math:`y`, this operation is equivalent to the BLAS routine ``SGER`` or ``DGER`` for single or double precision.

2. For complex vectors :math:`x` and :math:`y` with ``ArgTrans == KokkosBatched::Trans::Transpose``, this operation is equivalent to the BLAS routine ``CGERU`` or ``ZGERU`` for single or double precision.

3. For complex vectors :math:`x` and :math:`y` with ``ArgTrans == KokkosBatched::Trans::ConjTranspose``, this operation is equivalent to the BLAS routine ``CGERC`` or ``ZGERC`` for single or double precision.

Parameters
==========

:alpha: Scaling factor of the rank-1 update
:x: On input, :math:`x` is a length m vector.
:y: On input, :math:`y` is a length n vector.
:A: On input, :math:`A` is a m by n matrix being updated.

Type Requirements
-----------------

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`A = A + \alpha (x * y^T)`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`A = A + \alpha (x * y^H)`
- ``ScalarType`` must be a built-in floating point type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``)
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`X`
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`Y`
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`A` that satisfies ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_ger.cpp
  :language: c++

output:

.. code::

   ger works correctly!
