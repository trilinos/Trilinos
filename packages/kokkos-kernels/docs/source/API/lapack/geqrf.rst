KokkosLapack::geqrf
###################

Defined in header: :code:`KokkosLapack_geqrf.hpp`

.. code:: c++

  template <class ExecutionSpace, class AMatrix, class TauArray, class InfoArray>
  void geqrf(const ExecutionSpace& space, const AMatrix& A, const TauArray& Tau, const InfoArray& Info);

  template <class AMatrix, class TauArray, class InfoArray>
  void geqrf(const AMatrix& A, const TauArray& Tau, const InfoArray& Info);

Performs the QR factorization of matrix :math:`A`

.. math::

   A=Q*R

where :math:`A` is the input coefficient matrix on entry and the resulting :math:`R` factor and scaled Householder vectors on exit. :math:`Tau` stores the scaling factors associated with the Householder vectors.

1. Overwrites :math:`A` with the :math:`QR` factors using the resources of ``space``.
2. Same as 1. but uses the resources of the default execution space from ``AMatrix::execution_space``.

The function will throw a runtime exception if ``Tau.extent(0) < Kokkos::min(A.extent(0), A.extent(1))`` or if ``Info.extent(0) < 1``.

Parameters
==========

:space: execution space instance.

:A: The input matrix on entry and the :math:`QR` factors on return.

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

TBD

