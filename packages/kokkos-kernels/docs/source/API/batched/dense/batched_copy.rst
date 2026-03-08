KokkosBatched::Copy
###################

Defined in header: :code:`KokkosBatched_Copy_Decl.hpp`

.. code:: c++

  template <typename ArgTrans = Trans::NoTranspose>
  struct SerialCopy {
    template <typename AViewType, typename BViewType>
    KOKKOS_INLINE_FUNCTION static invoke(const AViewType &A, const BViewType &B);
  };

  template <typename MemberType, typename ArgTrans = Trans::NoTranspose>
  struct TeamCopy {
    template <typename AViewType, typename BViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B);
  };

  template <typename MemberType, typename ArgTrans = Trans::NoTranspose>
  struct TeamVectorCopy {
    template <typename AViewType, typename BViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B);
  };

Performs :math:`B = op(A)` where :math:`op(A)` is one of :math:`A`, :math:`A^T`, or :math:`A^H`.

1. For real vectors :math:`A` and :math:`B`, this operation is equivalent to the BLAS routine ``SCOPY`` or ``DCOPY`` for single or double precision.
2. For complex vectors :math:`A` and :math:`B`, this operation is equivalent to the BLAS routine ``CCOPY`` or ``ZCOPY`` for single or double precision.

Parameters
==========

:A: On input, :math:`A` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix.
:B: On input, :math:`B` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix. On output, :math:`B` is overwritten by the updated vector or matrix.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamCopy`` and ``TeamVectorCopy``)

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` for :math:`op(A) = A`
   - ``KokkosBatched::Trans::Transpose`` for :math:`op(A) = A^T`
   - ``KokkosBatched::Trans::ConjTranspose`` for :math:`op(A) = A^H`

- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector :math:`A`
- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector :math:`B` that satisfies ``std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type>``

.. note::

  This kernel supports both vector and matrix operations. When the input views :math:`A` and :math:`B` are of rank 1, the kernel performs a vector operation (BLAS copy).
  When the input views :math:`A` and :math:`B` are of rank 2, the kernel performs a matrix operation where the matrix :math:`A` is copied to :math:`B` with an optional transpose or conjugate transpose.
  The template argument to specify the rank of the input views is deprecated from 5.1.0.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_copy.cpp
  :language: c++

output:

.. code::

   copy works correctly!
