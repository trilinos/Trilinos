KokkosBlas::axpy
################

Defined in header: :code:`KokkosBlas1_axpby.hpp`

.. code:: c++

  template <class execution_space, class AV, class XMV, class YMV>
  void axpy(const execution_space& space, const AV& a, const XMV& X, const YMV& Y)

  template <class AV, class XMV, class YMV>
  void axpy(const AV& a, const XMV& X, const YMV& Y)

Add entries of :code:`X` scaled by coeffecient :code:`a` to entries of :code:`Y`: ``Y += a*X``

1. iterate over the entries of ``Y`` and add the corresponding entries of ``X`` scaled by ``a`` using the resources of ``space``
2. iterate over the entries of ``Y`` and add the corresponding entries of ``X`` scaled by ``a`` using the resources of the default instance of ``typename XMV::execution_space``

The function will throw a runtime exception if any of the following conditions are **not** met:
  - ``Y.extent(0) == X.extent(0) && Y.extent(1) == X.extent(1)``
  - ``Kokkos::is_view_v<AV> && (a.extent(0) == 1 || a.extent(0) == X.extent(1)``

Parameters
==========

:space: execution space instance

:a: scaling factor

:X, Y: vector(s) of input and output values respectively

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible == true``

- `YMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``YMV::rank == XMV::rank``
  - ``std::is_same_v<typename YMV::value_type, typename YMV::non_const_value_type> == true``
  - ``Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible == true``

- `AV` must be a scalar or a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``(AV::rank == XMV::rank - 1) || (AV::rank == 0)``
  - ``Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible == true``

Example
=======

.. code:: c++

  #include<Kokkos_Core.hpp>
  #include<KokkosBlas1_axpby.hpp>

  int main(int argc, char* argv[]) {
    Kokkos::initialize();
    {

      int N = atoi(argv[1]);

      Kokkos::View<double*> x("X",N);
      Kokkos::View<double*> y("Y",N);
      Kokkos::deep_copy(x,3.0);
      Kokkos::deep_copy(y,2.0);

      double alpha = 1.5;

      KokkosBlas::axpy(alpha,x,y);

      double sum = 0.0;
      Kokkos::parallel_reduce("CheckValue", N, KOKKOS_LAMBDA (const int& i, double& lsum) {
        lsum += y(i);
      },sum);

      printf("Sum: %lf Expected: %lf Diff: %e\n",sum,1.0*N*(2.0+1.5*3.0),sum-1.0*N);
    }

    Kokkos::finalize();
  }
