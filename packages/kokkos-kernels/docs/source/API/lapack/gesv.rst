KokkosLapack::gesv
##################

Defined in header: :code:`KokkosLapack_gesv.hpp`

.. code:: c++

  template <class ExecutionSpace, class AMatrix, class BXMV, class IPIVV>
  void gesv(const ExecutionSpace& space, const AMatrix& A, const BXMV& B, const IPIVV& IPIV);

  template <class AMatrix, class BXMV, class IPIVV>
  void gesv(const AMatrix& A, const BXMV& B, const IPIVV& IPIV);

Solve the dense system of linear equations

.. math::

   A*X=B

where :math:`A` is the input coefficient matrix, :math:`B` is the input right-hand-side vectors, and :math:`X` is the output solution vectors. Note that :math:`B` can contain multiple vectors to solve at a time. The solution is overwritten into `B`.

1. Overwrites the entries of :math:`B` with the solution of the linear system of equations using the resources of ``space``.
2. Same as 1. but use the resources of ``KernelHandle::HandleExecSpace{}``.

The function will throw a runtime exception if ``A.extent(0) < A.extent(1) || A.extent(0) != B.extent(0)``.

Parameters
==========

:space: execution space instance.

:A, B: The input matrix and right-hand-side vectors of the linear system. On return, B will hold the solution vectors.

:IPIV: Vector of pivots used to reorder the rows of :math:`A` while solving the system to improve numerical stability.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``

- `BXMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename BXMV::memory_space>::accessible``

- `IPIV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename IPIVV::memory_space>::accessible``

Example
=======

TBD

