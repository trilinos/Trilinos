KokkosLapack::gemqr
###################

Defined in header: :code:`KokkosLapack_gemqr.hpp`

.. code:: c++

  template <class ExecutionSpace, class AMatrix, class TauArray, class CMatrix, class InfoArray>
  void gemqr(const ExecutionSpace& space, const char side[], const char trans[], const AMatrix& A, const TauArray& Tau,
             const CMatrix& C, const InfoArray& Info);


  template <class AMatrix, class TauArray, class CMatrix, class InfoArray>
  void gemqr(const char side[], const char trans[], const AMatrix& A, const TauArray& Tau, const CMatrix& C,
             const InfoArray& Info);

Applies the `Q` factor from the QR factorization of matrix :math:`A` to matrix :math:`C` using the prescribed side and operation

.. math::

   C=Q*C\quad\text{or}\quadC=C*Q

where :math:`A` is a matrix previously factored using a call to ``geqrf`` and :math:`Tau` stores the associated scaling factors. 

1. Overwrites :math:`C` with the :math:`Q*C` using the resources of ``space``.
2. Same as 1. but uses the resources of the default execution space from ``CMatrix::execution_space``.

The function will throw a runtime exception if the size of :math:`C` is incompatible with that of :math:`Q`.

Parameters
==========

:space: execution space instance.

:side: control parameter specifying on which side the solver is applied, supported values are ``L, l`` for left side and ``R, r`` for right side.

:trans: control parameter specifying what operation on the entries of :math:`Q` should be performed. Supported values are ``N, n`` for nothing, ``T, t`` for transpose mode and ``C, c`` for conjugate transpose mode.

:A: The input matrix that contains the :math:`QR` factors from a previous call to ``geqrf``.

:Tau: rank-1 view of size min(M,N) that contains the scaling factors of the elementary reflectors.

:C: The matrix to which we multiply the :math:`Q` factor.

:Info: rank-1 view of integers and of size 1: Info[0] = 0: successful exit; Info[0] < 0: if equal to '-i', the i-th argument had an illegal value.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``

- `Tau` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Tau::memory_space>::accessible``

- `CMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename CMatrix::memory_space>::accessible``

- `Info` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Info::memory_space>::accessible``

Example
=======

.. literalinclude:: ../../../../example/docs/lapack/KokkosLapack_docs_gemqr.cpp
  :language: c++

