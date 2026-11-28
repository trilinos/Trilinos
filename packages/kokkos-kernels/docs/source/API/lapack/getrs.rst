KokkosLapack::getrs
###################

Defined in header: :code:`KokkosLapack_getrs.hpp`

.. code:: c++

  template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView>
  void getrs(const ExecutionSpace& space, const char trans[], const AMatrix& A, const IpivView& Ipiv, const BMatrix& B,
             const InfoView& Info);

  template <class AMatrix, class IpivView, class BMatrix, class InfoView>
  void getrs(const char trans[], const AMatrix& A, const IpivView& Ipiv, const BMatrix& B, const InfoView& Info);

.. math::

   A*X=B

where :math:`A` is a square matrix that store :math:`L` and :math:`U` factors obtained from `KokkosLapack::getrf`, :math:`B` is the right handside of the problem and :math:`X` is the left hand side solution that we are seeking.

1. Swap entries of :math:`B` using :math:`Ipiv` and apply two triangular solves for the lower and upper factors stored in :math:`A` to overwrite :math:`B` with the solution of to the system of equations. If a transpose or conjugate transpose mode is requested, the order of operations is reversed.
2. Same as 1. but uses the resources of the default execution space from ``AMatrix::execution_space``.

Parameters
==========

:space: execution space instance.

:trans: a transpose mode to apply to matrix A: "N" for no transpose, "T" for transpose and "C" for conjugate transpose.

:A: The input matrix that contains the :math:`L` and :math:`U` factors from a call to `getrf`.

:Ipiv: Pivots to apply to :math:`B`

:B: On input the right hand side vectors that gets overwritten on output by the solution vectors.

:Info: A scalar (stored in a device view) that holds the return code indicating if a zero pivot was detected.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``

- `IpivView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename IpivView::memory_space>::accessible``
  - ``std::is_integral_v<typename IpivView::non_const_value_type>``

- `BMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``
  - ``!std::is_const_v<typename BMatrix::value_type>``

- `InfoView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename Info::memory_space>::accessible``
  - ``std::is_same_v<typename InfoView::value_type, int>``

Example
=======

.. literalinclude:: ../../../../example/docs/lapack/KokkosLapack_docs_getrs.cpp
  :language: c++

output:
