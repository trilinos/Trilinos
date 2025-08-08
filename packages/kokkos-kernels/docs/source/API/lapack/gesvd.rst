KokkosLapack::svd
#################

Defined in header: :code:`KokkosLapack_svd.hpp`

.. code:: cppkokkos

  template <class ExecutionSpace, class AMatrix, class SVector, class UMatrix, class VMatrix>
  void svd(const ExecutionSpace& space, const char jobu[], const char jobvt[], const AMatrix& A,
           const SVector& S, const UMatrix& U, const VMatrix& Vt);

  template <class AMatrix, class SVector, class UMatrix, class VMatrix>
  void svd(const char jobu[], const char jobvt[], const AMatrix& A,
           const SVector& S, const UMatrix& U, const VMatrix& Vt);

Compute the singular value decomposition of matrix ``A`` which can be written as

.. math::

   A=U*S*transpose(V)

.. note::

   The function will compute the transpose of the right singular vectors ``Vt`` for computational efficiency.

1. Compute the singular value decomposition of ``A`` into ``S``, ``U`` and ``Vt`` using the resources associated with space.
2. Same as 1. but using the resources associated with ``AMatrix::execution_space()``.

Parameters
==========

:space: execution space instance.

:jobu, jobvt: characters used to control the calculation of left (jobu) and right (jobvt) singular vectors. ``A`` means all the singular vectors are computed, ``S`` means the first :math:`min(m, n)` vectors will be computed in ``U``, ``O`` means the first :math:`min(m, n)` vectors will be computed and overwritten in ``A``, ``N`` means no singular vectors will be computed.

:A: the matrix on which the singular value decomposition will be computed.

:S: the ``min(m, n)`` singular values of ``A`` in decreasing order.

:U, Vt: the left and right singular vectors of ``A`` computed according to the ``jobu`` and ``jobvt`` flags.

Type Requirements
=================

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- the three matrix types ``AMatrix``, ``UMatrix`` and ``VMatrix`` have the same requirements:

  - ``Kokkos::is_view_v<AMatrix> && AMatrix::rank() == 2 && Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible``
  - ``Kokkos::is_view_v<UMatrix> && UMatrix::rank() == 2 && Kokkos::SpaceAccessibility<ExecutionSpace, typename UMatrix::memory_space>::accessible``
  - ``Kokkos::is_view_v<VMatrix> && VMatrix::rank() == 2 && Kokkos::SpaceAccessibility<ExecutionSpace, typename VMatrix::memory_space>::accessible``

- the vector of singular values has the following requirements:

  - ``Kokkos::is_view_v<SVector> && SVector::rank() == 1 && Kokkos::SpaceAccessibility<ExecutionSpace, typename SVector::memory_space>::accessible``
  - ``std::is_same_v<typename SVector::array_layout, Kokkos::LayoutStride> && S.span_is_contiguous()`` only contiguous strided spans are allowed.

.. note::

  Should check that the ouput view ``SVector``, ``UMatrix`` and ``VMatrix`` are storing non-const data.

Example
=======
