KokkosBatched::Symv
###################

Defined in header: :code:`KokkosBatched_Symv.hpp`

.. code:: c++

    template <typename ArgUplo, typename ArgTrans>
    struct SerialSymv {
      template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const XViewType &x,
                                               const ScalarType beta, const YViewType &y);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamSymv {
      template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                               const XViewType &x, const ScalarType beta, const YViewType &y);
    };

    template <typename MemberType, typename ArgUplo, typename ArgTrans>
    struct TeamVectorSymv {
      template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                               const XViewType &x, const ScalarType beta, const YViewType &y);
    };

Performs one of the symmetric or hermitian matrix-vector multiplication

.. math::

   \begin{align}
   Y &= \alpha A * x + \beta Y
   \end{align}

1. If ``ArgTrans == KokkosBatched::Trans::Transpose``, this operation is equivalent to the BLAS routine `SSYMV <https://www.netlib.org/blas/ssymv.f>`_ (`CSYMV <https://www.netlib.org/blas/csymv.f>`_) or `DSYMV <https://www.netlib.org/blas/dsymv.f>`_ (`ZSYMV <https://www.netlib.org/blas/zsymv.f>`_) for single or double precision for real (complex) symmetric matrix.

2. If ``ArgTrans == KokkosBatched::Trans::ConjTranspose``, this operation is equivalent to the BLAS routine `CHEMV <https://www.netlib.org/blas/chemv.f>`_ or `ZHEMV <https://www.netlib.org/blas/zhemv.f>`_ for single or double precision for complex hermitian matrix.

Parameters
==========

:alpha: Scaling factor
:A: :math:`A` is a symmetric (or hermitian) matrix. Only the upper or lower triangular part of :math:`A` is referenced based on ``ArgUplo``.
:x: :math:`x` is a length n vector.
:beta: Scaling factor
:y: :math:`y` is a length n vector. Before entry, y must contain the vector y. On exit, y is overwritten by the result :math:`\alpha A * x + \beta Y`

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamSymv`` and ``TeamVectorSymv``)
- ``ArgUplo`` must be one of the following:
   - ``KokkosBatched::Uplo::Upper`` to update the upper triangular part of :math:`A`
   - ``KokkosBatched::Uplo::Lower`` to update the lower triangular part of :math:`A`
- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::Transpose`` to perform :math:`y = \alpha A * x + \beta y` on the symmetric matrix :math:`A`
   - ``KokkosBatched::Trans::ConjTranspose`` to perform :math:`y = \alpha A * x + \beta y` on the hermitian matrix :math:`A`
- ``ScalarType`` must be a built-in floating point type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``)
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing a matrix :math:`A`
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`X`
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`Y` that satisfies ``std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_symv.cpp
  :language: c++

output:

.. code::

   symv works correctly!
