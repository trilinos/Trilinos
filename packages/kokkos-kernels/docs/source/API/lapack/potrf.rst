KokkosLapack::potrf
###################

Defined in header: :code:`KokkosLapack_potrf.hpp`

.. code:: c++

  template <class execution_space, class AViewType>
  void potrf([[maybe_unused]] const execution_space& space, const char uplo[], AViewType& A) {

Computes the Cholesky factorization of a complex Hermitian positive definite matrix A :math:`A`

.. math::

    A = U**H * U,  if UPLO = 'U', or
    A = L  * L**H,  if UPLO = 'L', A=Q*R

where :math:`A` is the input matrix and is the factor U or L from the Cholesky factorization :math:`A = U**H*U or A = L*L**H` on exit.

1. Overwrites :math:`A` with the Cholesky factorization using the resources of ``space``.

Parameters
==========

:space: execution space instance.

:uplo: 'U':  Upper triangle of A is stored, else lower triangle

:A: The input matrix (lda,n) on entry and the Cholesky factorization on return.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AMatrix` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies:
  - ``AMatrix::rank == 2``, i.e., ``AMatrix`` represents a matrix,
  - ``AMatrix::value_type`` is a supported scalar type (for example, ``float``, ``double``, or a corresponding ``Kokkos::complex`` type),
  - ``AMatrix::array_layout`` is ``Kokkos::LayoutLeft``,
  - the memory space of ``AMatrix`` is accessible from ``ExecutionSpace``, i.e. ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible`` is ``true``.

Example
=======

.. literalinclude:: ../../../../example/docs/lapack/KokkosLapack_docs_potr.cpp
  :language: c++
