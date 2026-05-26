KokkosBatched::Axpy
###################

Defined in header: :code:`KokkosBatched_Axpy.hpp`

.. code:: c++

  struct SerialAxpy {
    template <typename XViewType, typename YViewType, typename alphaViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const alphaViewType &alpha, const XViewType &X, const YViewType &Y);
  };

  template <typename MemberType>
  struct TeamAxpy {
    template <typename XViewType, typename YViewType, typename alphaViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha, const XViewType &X,
                                             const YViewType &Y);
  };

  template <typename MemberType>
  struct TeamVectorAxpy {
    template <typename XViewType, typename YViewType, typename alphaViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha, const XViewType &X,
                                             const YViewType &Y);
  };

Add entries of :math:`X` scaled by :math:`\alpha` to entries of :math:`Y`

.. math::

   \begin{align}
   Y &= Y + \alpha X
   \end{align}

1. For real vectors :math:`X` and :math:`Y`, this operation is equivalent to the BLAS routine ``SAXPY`` or ``DAXPY`` for single or double precision.
2. For complex vectors :math:`X` and :math:`Y`, this operation is equivalent to the BLAS routine ``CAXPY`` or ``ZAXPY`` for single or double precision.

Parameters
==========

:alpha: A scalar scaling factor or a length :math:`m` vector of scaling factors.
:x: On input, :math:`X` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix.
:y: On input, :math:`Y` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix. On output, :math:`Y` is overwritten by the updated vector or matrix.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamAxpy`` and ``TeamVectorAxpy``).
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector :math:`X`
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector :math:`Y` that satisfies ``std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>``
- ``alphaViewType`` must be either a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a matrix :math:`\alpha` or a built-in arithmetic type like ``float``, ``double``, ``Kokkos::complex<float>``, or ``Kokkos::complex<double>``

.. note::

  This kernel supports both vector and matrix operations. When the input views :math:`X` and :math:`Y` are of rank 1, the kernel performs a vector operation (BLAS axpy).
  When the input views :math:`X` and :math:`Y` are of rank 2, the kernel performs a matrix operation where each row of the matrices is treated as a separate vector.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_axpy.cpp
  :language: c++

output:

.. code::

   axpy works correctly!
