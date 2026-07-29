KokkosBatched::Nrm
##################

Defined in header: :code:`KokkosBatched_Nrm.hpp`

.. code:: c++

  template <typename NrmType>
  struct SerialNrm {
    template <typename XViewType, typename NormViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X, const NormViewType &norm);
  };

  template <typename MemberType, typename NrmType>
  struct TeamNrm {
    template <typename XViewType, typename NormViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const NormViewType &norm);
  };

  template <typename MemberType, typename NrmType>
  struct TeamVectorNrm {
    template <typename XViewType, typename NormViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const NormViewType &norm);
  };

Computes the :math:`L1`, :math:`L2` or :math:`L_\infty` norm of a vector :math:`X`.

.. math::

   \begin{align}
   norm &= ||x|| \: \text{(if NrmType == KokkosBatched::Norm::L1)} \\
   norm &= ||x||_2 \: \text{(if NrmType == KokkosBatched::Norm::L2 or NrmType == KokkosBatched::Norm::ScaledL2)} \\
   norm &= ||x||_\infty \: \text{(if NrmType == KokkosBatched::Norm::LInf)}
   \end{align}

1. If ``NrmType == KokkosBatched::Norm::L1``, this operation is equivalent to the BLAS routine `SASUM <https://www.netlib.org/blas/sasum.f>`_ (`SCASUM <https://www.netlib.org/blas/scasum.f>`_) or `DASUM <https://www.netlib.org/blas/dasum.f>`_ (`DZASUM <https://www.netlib.org/blas/dzasum.f>`_) for single or double precision for real (complex) vectors.
2. If ``NrmType == KokkosBatched::Norm::L2`` or ``NrmType == KokkosBatched::Norm::ScaledL2``, this operation is equivalent to the BLAS routine `SNRM2 <https://www.netlib.org/lapack/lapack-3.1.1/html/snrm2.f.html>`_ (`SCNRM2 <https://www.netlib.org/lapack/lapack-3.1.1/html/scnrm2.f.html>`_) or `DNRM2 <https://www.netlib.org/lapack/lapack-3.1.1/html/dnrm2.f.html>`_ (`DZNRM2 <https://www.netlib.org/lapack/lapack-3.1.1/html/dznrm2.f.html>`_) for single or double precision for real (complex) vectors.
3. If ``NrmType == KokkosBatched::Norm::LInf``, this operation is related to the BLAS routine `ISAMAX <https://www.netlib.org/blas/isamax.f>`_ (`ICAMAX <https://www.netlib.org/blas/icamax.f>`_) or `IDAMAX <https://www.netlib.org/blas/idamax.f>`_ (`IZAMAX <https://www.netlib.org/blas/izamax.f>`_) for single or double precision for real (complex) vectors, where the index of the maximum absolute value is returned in the output :math:`norm`. This routine returns the maximum absolute value instead.

.. note::
  
  Though ``NrmType == KokkosBatched::Norm::L2`` is more efficient, it may overflow for large vectors. For large vectors, ``NrmType == KokkosBatched::Norm::ScaledL2`` is recommended as it uses a numerically stable algorithm to compute the :math:`L2` norm that avoids overflow and underflow by scaling the input vector :math:`X` by the maximum absolute value of its elements that has been encountered.

Parameters
==========

:X: On input, :math:`X` is a length :math:`n` vector
:norm: On output, :math:`norm` is the computed norm of the vector :math:`X`.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamNrm`` and ``TeamVectorNrm``)

- ``NrmType`` must be one of the following:
   - ``KokkosBatched::Norm::L1`` for :math:`L1` norm
   - ``KokkosBatched::Norm::L2`` or ``KokkosBatched::Norm::ScaledL2`` for :math:`L2` norm
   - ``KokkosBatched::Norm::LInf`` for :math:`L_\infty` norm

- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector or matrix :math:`X`
- ``NormViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 containing the output :math:`norm`. The norm is accumulated in the type of the elements of ``NormViewType``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_nrm.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   nrm works correctly!
