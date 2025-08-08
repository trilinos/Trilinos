KokkosBlas::gemm
################

Defined in header: :code:`KokkosBlas3_gemm.hpp`

.. code:: cppkokkos

  template <class execution_space, class AViewType, class BViewType, class CViewType>
  void gemm(const execution_space& space, const char transA[], const char transB[],
            typename AViewType::const_value_type& alpha, const AViewType& A, const BViewType& B,
            typename CViewType::const_value_type& beta, const CViewType& C);

  template <class AViewType, class BViewType, class CViewType>
  void gemm(const char transA[], const char transB[], typename AViewType::const_value_type& alpha, const AViewType& A,
            const BViewType& B, typename CViewType::const_value_type& beta, const CViewType& C);

Perform a dense matrix-matrix multiplication

.. math::

   C = \beta * C + \alpha * A * B


1. Iterate over the entries of ``C``, scale them with ``beta``, compute the corresponding dot product of a row of ``A`` with a column of ``B``, scale the result with ``alpha`` and accumulate back into ``C`` using the resources from the ``space`` instance, or
2. Iterate over the entries of ``C``, scale them with ``beta``, compute the corresponding dot product of a row of ``A`` with a column of ``B``, scale the result with ``alpha`` and accumulate back into ``C`` using resources from the default instance of ``typename CViewType::execution_space``.

The function will throw a runtime exception if any of the following conditions are not met:

- ``transA`` stores a valid control value, see below.
- ``transB`` stores a valid control value, see below.
- ``A.extent(1) == B.extent(0) && A.extent(0) == C.extent(0) && B.extent(1) == C.extent(1)``

Parameters
==========

:space: execution space instance.

:transA, transB: control parameters to apply transformations to ``A`` and ``B`` respectively, the accepted values are "N" for no transpose, "T" for transpose and "C" for conjugate transpose, all characters after the first one are ignored.

:alpha, beta: scaling parameters for ``A*B`` and ``C``, respectively.

:A, B, C: the two matrices to multiply and the matrix where the scaled multiplication values will be accumulated.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename AViewType::memory_space>::accessible``

- `BViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename BViewType::memory_space>::accessible``

- `CViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename CViewType::memory_space>::accessible``

Example
=======

.. code:: cppkokkos

  #include<Kokkos_Core.hpp>
  #include<KokkosBlas3_gemm.hpp>

  int main(int argc, char* argv[]) {
    Kokkos::initialize();
    {
      int M = atoi(argv[1]);
      int N = atoi(argv[2]);

      Kokkos::View<double**> A("A",M,N);
      Kokkos::View<double**> B("B",N,M);
      Kokkos::View<double**> C("C",M,M);
   
      Kokkos::deep_copy(A,1.0);
      Kokkos::deep_copy(B,2.0);

      const double alpha = double(1.0);
      const double beta = double(0.0);
   
      KokkosBlas::gemm("N","N",alpha,A,B,beta,C);
    }
    Kokkos::finalize();
  }
