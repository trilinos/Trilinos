KokkosLapack::gegqr
###################

Defined in header: :code:`KokkosLapack_gegqr.hpp`

.. code:: c++

  template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
  void gegqr(const ExecutionSpace& space, const int k, const AMatrix& A, const TauArray& Tau, const InfoArray& Info);

  template <class AMatrix, class TauArray, class InfoArray>
  void gegqr(const int k, const AMatrix& A, const TauArray& Tau, const InfoArray& Info);

Computes the matrix :math:`Q` from the QR factorization of matrix :math:`A` using the first :math:`k` reflectors

.. math::

   Q=H(0)H(1)...H(k-1)

where :math:`A` is, on input, a matrix previously factored using a call to ``geqrf`` and contains :math:`Q` on output and :math:`Tau` stores the associated scaling factors. 

1. Overwrites :math:`A` with the :math:`Q` using the resources of ``space``.
2. Same as 1. but uses the resources of the default execution space from ``AMatrix::execution_space``.

The function will throw a runtime exception if :math:`k < 0` or :math:`n < k` where :math:`n` is the number of columns in :math:`A`.

Parameters
==========

:space: execution space instance.

:k: The number of reflectors to use when computing :math:`Q`.

:A: On input, matrix that contains the :math:`QR` factors from a previous call to ``geqrf``, on output, it contains the first :math:`k` columns of :math:`Q`.

:Tau: rank-1 view of size min(M,N) that contains the scaling factors of the elementary reflectors.

:Info: rank-1 view of integers and of size 1: Info[0] = 0: successful exit; Info[0] < 0: if equal to '-i', the i-th argument had an illegal value.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``

- `Tau` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Tau::memory_space>::accessible``

- `Info` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Info::memory_space>::accessible``

Example
=======

.. literalinclude:: ../../../../example/docs/lapack/KokkosLapack_docs_gegqr.cpp
  :language: c++
