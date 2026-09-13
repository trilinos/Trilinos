KokkosLapack::getrf
###################

Defined in header: :code:`KokkosLapack_getrf.hpp`

.. code:: c++

  template <class ExecutionSpace, class AMatrix, class IpivView, class InfoView>
  void getrf(const ExecutionSpace& space, const AMatrix& A, const IpivView& Ipiv, const InfoView& Info);

  template <class AMatrix, class IpivView, class InfoView>
  void getrf(const AMatrix& A, const IpivView& Ipiv, const InfoView& Info);

.. math::

   P*A=L*U

where :math:`A` is, on input, the matrix storing the linear system of equations and on output, contains the factors L and U in its lower and upper triangular parts. :math:`Ipiv` contains the pivots used to construct the row permutation matrix :math:`P`.

1. Overwrite :math:`A` with the :math:`L` and :math:`U` factors using the resources associated with ``space``.
2. Same as 1. but uses the resources of the default execution space from ``AMatrix::execution_space``.

Parameters
==========

:space: execution space instance.

:A: The input matrix that gets overwritten with the :math:`L` and :math:`U` factors on output.

:Ipiv: The vector holding the pivots selected during the matrix decomposition.

:Info: A scalar (store in a device view) that holds the return code from Lapack when that library gets called.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``

- `IpivView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Tau::memory_space>::accessible``

- `InfoView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Info::memory_space>::accessible``

Example
=======

.. literalinclude:: ../../../../example/docs/lapack/KokkosLapack_docs_geqrf.cpp
  :language: c++

output:

.. code::

   KokkosLapack::getrf() returned correct results!
