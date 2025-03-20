KokkosBlas::gemv
################

Defined in header: :code:`KokkosBlas2_gemv.hpp`

.. code:: cppkokkos

  template <class ExecutionSpace, class AViewType, class XViewType, class YViewType>
  void gemv(const ExecutionSpace& space, const char trans[], typename AViewType::const_value_type& alpha,
            const AViewType& A, const XViewType& x, typename YViewType::const_value_type& beta, const YViewType& y);

  template <class AViewType, class XViewType, class YViewType>
  void gemv(const char trans[], typename AViewType::const_value_type& alpha, const AViewType& A, const XViewType& x,
            typename YViewType::const_value_type& beta, const YViewType& y);

Perform a dense matrix-vector multiplication

.. math::

   y = \beta * y + \alpha * A * x

1. Iterate over the entries of ``y``, scale them with ``beta``, compute the dot product of the corresponding row of ``A`` with ``x``, scale it with ``alpha`` and accumulate to ``y`` using resources from the ``space`` instance, or
2. Iterate over the entries of ``y``, scale them with ``beta``, compute the dot product of the corresponding row of ``A`` with ``x``, scale it with ``alpha`` and accumulate to ``y`` using resources from the default instance of ``typename AViewType::execution_space``.

The function will throw a runtime exception if any of the following conditions are not met:

- ``trans`` stores a valid control value, see below
- ``trans == N && A.extent(1) == x.extent(0) && A.extent(0) == y.extent(0)``
- ``trans == T,C && A.extent(0) == x.extent(0) && A.extent(1) == y.extent(0)``

Parameters
==========

:space: execution space instance.

:trans: "N" for no transpose, "T" for transpose and "C" for conjugate transpose of ``A``, all characters after the first one are ignored.

:alpha, beta: scaling parameters for ``A*x`` and ``y``, respectively.

:A: matrix being multiplied with ``x``.

:x, y: input and output vectors for the matrix-vector multiplication, respectively.


Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `AViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AViewType::memory_space>::accessible``

- `XViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename XViewType::memory_space>::accessible``

- `YViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename YViewType::memory_space>::accessible``

Example
=======

.. code:: cppkokkos
  
  #include<Kokkos_Core.hpp>
  #include<KokkosBlas2_gemv.hpp>

  int main(int argc, char* argv[]) {
    Kokkos::initialize();
    {
      int M = atoi(argv[1]);
      int N = atoi(argv[2]);

      Kokkos::View<double**> A("A",M,N);
      Kokkos::View<double*> x("X",N);
      Kokkos::View<double*> y("Y",N);

      Kokkos::deep_copy(A,1.0);
      Kokkos::deep_copy(x,3.0);

      const double alpha = double(1.0);
      const double beta = double(0.0);
   
      KokkosBlas::gemv("N",alpha,A,x,beta,y);
    }
    Kokkos::finalize();
  }
