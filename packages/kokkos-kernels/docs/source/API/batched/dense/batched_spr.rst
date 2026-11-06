KokkosBatched::Spr
##################

Defined in header: :code:`KokkosBatched_Spr.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans>
    struct SerialSpr {
      template <typename ScalarType, typename XViewType, typename APViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const APViewType &ap);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamSpr {
      template <typename ScalarType, typename XViewType, typename APViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const XViewType &x,
                                               const APViewType &ap);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamVectorSpr {
      template <typename ScalarType, typename XViewType, typename APViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const XViewType &x,
                                               const APViewType &ap);
    };

Performs one of the symmetric or hermitian rank-1 operation of matrix :math:`A` by vector :math:`x` with scaling factor :math:`alpha`, where :math:`A` is stored in packed format.

.. math::

   \begin{align}
   A &= A + \alpha (x * x^T) \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   A &= A + \alpha (x * x^H) \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. If ``ArgTrans == KokkosBatched::Trans::Transpose``, this operation is equivalent to the BLAS routine `SSPR <https://www.netlib.org/blas/sspr.f>`_ (`CSPR <https://www.netlib.org/blas/cspr.f>`_) or `DSPR <https://www.netlib.org/blas/dspr.f>`_ (`ZSPR <https://www.netlib.org/blas/zspr.f>`_) for single or double precision for real (complex) symmetric matrix.
   
2. If ``ArgTrans == KokkosBatched::Trans::ConjTranspose``, this operation is equivalent to the BLAS routine `CHPR <https://www.netlib.org/blas/chpr.f>`_ or `ZHPR <https://www.netlib.org/blas/zhpr.f>`_ for single or double precision for complex hermitian matrix.

Parameters
==========

:alpha: Scaling factor
:x: On input, :math:`x` is a length :math:`n` vector.
:ap: With ``ArgUplo == KokkosBatched::Uplo::Upper``, the upper triangular part of :math:`A` is referenced and updated. With ``ArgUplo == KokkosBatched::Uplo::Lower``, the lower triangular part of :math:`A` is referenced and updated. The matrix :math:`A` is supplied in the packed format as a length :math:`n*(n+1)/2` vector, where :math:`n` is the size of the matrix.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamSpr`` and ``TeamVectorSpr``)
- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` to update the upper triangular part of :math:`A`
   - ``KokkosBatched::Uplo::Lower`` to update the lower triangular part of :math:`A`
- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`A = A + \alpha (x * x^T)`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`A = A + \alpha (x * x^H)`
- ``ScalarType`` must be a built-in floating point type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``)
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`X`
- ``APViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a matrix :math:`A` in packed format  that satisfies ``std::is_same_v<typename APViewType::value_type, typename APViewType::non_const_value_type>``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_spr.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   spr works correctly!
