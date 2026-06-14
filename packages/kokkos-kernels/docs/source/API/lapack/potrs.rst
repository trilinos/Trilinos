KokkosLapack::potrs
###################

Defined in header: :code:`KokkosLapack_potrs.hpp`

.. code:: c++

  // Overload 1: with explicit execution space
  template <class ExecutionSpace, class AViewType, class BViewType>
  void potrs(const ExecutionSpace& space, const char uplo[], const AViewType& A, BViewType& B);

  // Overload 2: uses AViewType::execution_space
  template <class AViewType, class BViewType>
  void potrs(const char uplo[], const AViewType& A, BViewType& B);

Solves a system of linear equations :math:`A X = B` where :math:`A` is a
symmetric (or Hermitian) positive definite matrix whose Cholesky factorization
has already been computed by :code:`KokkosLapack::potrf`.

Given the factorization produced by ``potrf``:

.. math::

    A = U^H U \quad \text{if uplo = 'U'}, \qquad
    A = L L^H \quad \text{if uplo = 'L'},

``potrs`` solves for :math:`X` by two triangular solves, overwriting :math:`B`
with the solution :math:`X`.

Parameters
==========

:space: execution space instance (overload 1 only).

:uplo: ``'U'`` if the upper triangular factor is stored in ``A``; ``'L'`` if the lower triangular factor is stored.

:A: On entry, the triangular Cholesky factor of :math:`A` as returned by ``potrf`` — upper triangle if ``uplo = 'U'``, lower triangle if ``uplo = 'L'``.  Not modified.  Dimensions: ``(lda, n)``.

:B: On entry, the right-hand side matrix of dimensions ``(ldb, nrhs)``.  On exit, overwritten with the solution matrix :math:`X`.

Type Requirements
=================

- ``ExecutionSpace`` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ satisfying:

  - ``AViewType::rank == 2``,
  - ``AViewType::value_type`` is ``const Scalar`` for a supported scalar type (``float``, ``double``, ``Kokkos::complex<float>``, ``Kokkos::complex<double>``),
  - ``AViewType::array_layout`` is ``Kokkos::LayoutLeft``,
  - the memory space of ``AViewType`` is accessible from ``ExecutionSpace``.

- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ satisfying:

  - ``BViewType::rank == 2``,
  - ``BViewType::value_type`` is a (non-const) scalar type matching that of ``AViewType``,
  - ``BViewType::array_layout`` is ``Kokkos::LayoutLeft``,
  - the memory space of ``BViewType`` is accessible from ``ExecutionSpace``.

Example
=======

.. literalinclude:: ../../../../example/docs/lapack/KokkosLapack_docs_potr.cpp
  :language: c++
