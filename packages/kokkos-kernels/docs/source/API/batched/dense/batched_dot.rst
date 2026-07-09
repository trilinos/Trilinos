KokkosBatched::Dot
##################

Defined in header: :code:`KokkosBatched_Dot.hpp`

.. code:: c++

  template <typename ArgTrans, int Axis>
  struct SerialDot {
    template <typename XViewType, typename YViewType, typename NormViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X, const YViewType &Y, const NormViewType &dot);
  };

  template <typename MemberType, typename ArgTrans, int Axis>
  struct TeamDot {
    template <typename XViewType, typename YViewType, typename NormViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y, const NormViewType &dot);
  };

  template <typename MemberType, typename ArgTrans, int Axis>
  struct TeamVectorDot {
    template <typename XViewType, typename YViewType, typename NormViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y, const NormViewType &dot);
  };

Performs the dot product of two vectors :math:`X` and :math:`Y`.

.. math::

   \begin{align}
   dot &= X^T * Y \: \text{(if ArgTrans == KokkosBatched::Trans::Transpose)} \\
   dot &= X^H * Y \: \text{(if ArgTrans == KokkosBatched::Trans::ConjTranspose)}
   \end{align}

1. If ``ArgTrans == KokkosBatched::Trans::Transpose``, this operation is equivalent to the BLAS routine `SDOT <https://www.netlib.org/blas/sdot.f>`_ (`CDOTU <https://www.netlib.org/blas/cdotu.f>`_) or `DDOT <https://www.netlib.org/blas/ddot.f>`_ (`ZDOTU <https://www.netlib.org/blas/zdotu.f>`_) for single or double precision for real (complex) vectors.
2. If ``ArgTrans == KokkosBatched::Trans::ConjTranspose``, this operation is equivalent to the BLAS routine `CDOTC <https://www.netlib.org/blas/cdotc.f>`_ or `ZDOTC <https://www.netlib.org/blas/zdotc.f>`_ for single or double precision for complex vectors.

.. note::
  
  This kernel does not support the BLAS routine `SDSDOT <https://www.netlib.org/blas/sdsdot.f>`_ which returns the single precision dot product of two single precision vectors with dot product accumulated in double precision.
  For `DSDOT <https://www.netlib.org/blas/dsdot.f>`_, provide :math:`X` and :math:`Y` in single precision and provide the output :math:`dot` in double precision.

Parameters
==========

:X: On input, :math:`X` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix.
:Y: On input, :math:`Y` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix.
:dot: On output, :math:`dot` is the computed dot product if :math:`X` and :math:`Y` are vectors, or the computed dot products along the specified axis if :math:`X` and :math:`Y` are matrices.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamDot`` and ``TeamVectorDot``)

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::Transpose`` for :math:`dot = X^T * Y`
   - ``KokkosBatched::Trans::ConjTranspose`` for :math:`dot = X^H * Y`

- ``Axis`` must be one of the following:
   - ``0`` to perform the operation along the first dimension (columns) when :math:`X` and :math:`Y` are matrices
   - ``1`` to perform the operation along the second dimension (rows) when :math:`X` and :math:`Y` are matrices

- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector or matrix :math:`X`
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector or matrix :math:`Y`
- ``NormViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 or 1 containing the output :math:`dot`. dot product is accumulated is accumulated in the type of the elements of ``NormViewType``

.. note::
  
  This kernel supports both vector and matrix operations. When the input views :math:`X` and :math:`Y` are of rank 1, the kernel performs a vector operation (BLAS dot). ``Axis`` must be set to 0 for this case.
  When the input views :math:`X` and :math:`Y` are of rank 2, the kernel performs a vector operation along the specified axis (0 or 1), where each column or row is treated as a separate vector.
  The template argument ``Axis`` to specify the axis to perform the operation is required from 5.2.0.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_dot.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   dot works correctly!
