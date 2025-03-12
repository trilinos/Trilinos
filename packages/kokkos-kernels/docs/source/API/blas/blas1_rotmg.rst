KokkosBlas::rotmg
#################

Defined in header: :code:`KokkosBlas1_rotmg.hpp`

.. code:: c++

  template <class execution_space, class DXView, class YView, class PView>
  void rotmg(execution_space const& space, DXView const& d1, DXView const& d2,
             DXView const& x1, YView const& y1, PView const& param);
  
  template <class DXView, class YView, class PView>
  void rotmg(DXView const& d1, DXView const& d2, DXView const& x1,
             YView const& y1, PView const& param);

1. Computes the modified Givens rotation coefficients that zeroes out the second entry of vector :math:`[x1\sqrt{d1}, y1\sqrt{d2}]`, using the provided execution space instance.
2. Computes the modified Givens rotation coefficients that zeroes out the second entry of vector :math:`[x1\sqrt{d1}, y1\sqrt{d2}]`, using the default instance of :code:`PView::execution_space`.

The output flag parameter (stored in :code:`params(0)`) along with the other four elements of :code:`param` determine the matrix :math:`H` that applies the rotation:

flag = -1.0:

.. math::

   H=\begin{bmatrix} h_{11} && h_{12} \\ h_{21} && h_{22}\end{bmatrix}

flag = 0.0:

.. math::

   H=\begin{bmatrix} 1 && h_{12} \\ h_{21} && 1\end{bmatrix}

flag = 1.0:

.. math::

   H=\begin{bmatrix} h_{11} && 1 \\ -1 && h_{22}\end{bmatrix}

flag = -2.0:

.. math::

   H=\begin{bmatrix} 1 && 0\\ 0  && 1\end{bmatrix}

where elements :math:`1 \ldots 4` of :code:`params` contain :math:`h_{11}, h_{21}, h_{12}, h_{22}`.

Note that the values -1.0, 0.0 and 1.0 implied by the flag are not stored explicitly in :code:`params`.

Parameters
==========

:space: execution space instance
:d1, d2: Rank-0 views. On input, they contain the squares of the initial scaling factors for x,y respectively. On output, they contain the squares of the scaling factors to be applied to x,y after applying the rotation with :doc:`rotm <blas1_rotm>`.
:x1, y1: Rank-0 views containing the components of the vector to rotate.
:params: A 1D vector with compile-time dimension 5. On output, :code:`params(0)` is the control flag and :code:`params(1...4)` are the coefficients of the rotation matrix :math:`H` as described above.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `DXView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename DXView::memory_space>::accessible == true``
  - ``std::is_same_v<typename DXView::non_const_value_type, typename DXView::value_type> == true``
  - ``!Kokkos::ArithTraits<typename DXView::value_type>::is_complex``

- `YView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename YView::memory_space>::accessible == true``
  - ``!Kokkos::ArithTraits<typename YView::value_type>::is_complex``

- `PView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 and compile-time extent 5 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename PView::memory_space>::accessible == true``
  - ``std::is_same_v<typename PView::non_const_value_type, typename PView::value_type> == true``
  - ``!Kokkos::ArithTraits<typename PView::value_type>::is_complex``

Example
=======

This example shows how to eliminate an entry using a modified Givens rotation.
It uses :code:`rotmg` to compute the rotation parameters and :doc:`rotm <blas1_rotm>` to apply the rotation.

.. code:: c++

    #include <iostream>
    #include <Kokkos_Core.hpp>
    #include <Kokkos_Random.hpp>
    #include "KokkosBlas1_rotmg.hpp"
    #include "KokkosBlas1_rotm.hpp"
    #include "KokkosKernels_PrintUtils.hpp"
    
    using execution_space = Kokkos::DefaultExecutionSpace;
    using Scalar          = double;
    using Vector          = Kokkos::View<Scalar*, execution_space>;
    using ParamView       = Kokkos::View<Scalar[5], execution_space>;
    using ScalarView      = Kokkos::View<Scalar, execution_space>;
    
    int main(int argc, char* argv[]) {
      Kokkos::initialize();
      {
        const int N = 10;
        Vector x("x", N);
        Vector y("y", N);
        ScalarView d1("d1");
        ScalarView d2("d2");
        ParamView param("param");
    
        // Populate x,y with uniform random values between 0 and 10
        Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
        Kokkos::fill_random(x, rand_pool, Scalar(10));
        Kokkos::fill_random(y, rand_pool, Scalar(10));
    
        // Populate input vector scaling factors with 1
        Kokkos::deep_copy(d1, Scalar(1));
        Kokkos::deep_copy(d2, Scalar(1));
    
        std::cout << "x,y before applying modified Givens rotation:\n";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
        KokkosKernels::Impl::kk_print_1Dview(std::cout, y);
    
        // Calculate Givens rotation coefficients to eliminate y(0)
        KokkosBlas::rotmg<execution_space, ScalarView, ScalarView, ParamView>(
            execution_space(), d1, d2, Kokkos::subview(x, 0), Kokkos::subview(y, 0), param);
    
        auto paramHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), param);
    
        std::cout << "\nrotmg output (rotation parameters) to eliminate y(0):\n";
        std::cout << "d1 = ";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, d1);
        std::cout << "d2 = ";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, d2);
        std::cout << "flag = " << paramHost(0) << '\n';
        std::cout << "h components = ";
        for (int i = 0; i < 4; i++) std::cout << paramHost(1 + i) << " ";
        std::cout << '\n';
    
        // Zero out y(0), which was left unmodified by rotmg.
        Kokkos::deep_copy(Kokkos::subview(y, 0), Scalar(0));
    
        // Apply the rotation to the remaining entries of x and y
        KokkosBlas::rotm(execution_space(), Kokkos::subview(x, Kokkos::make_pair(1, N)),
                         Kokkos::subview(y, Kokkos::make_pair(1, N)), param);
    
        // Apply scaling factors: sqrt(d1) and sqrt(d2) to x and y respectively
        Kokkos::parallel_for(
            Kokkos::RangePolicy<execution_space>(0, N), KOKKOS_LAMBDA(int i) {
              x(i) *= Kokkos::sqrt(d1());
              y(i) *= Kokkos::sqrt(d2());
            });
    
        std::cout << "\nx,y after applying modified Givens rotation and scaling by [sqrt(d1), sqrt(d2)]):\n";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
        KokkosKernels::Impl::kk_print_1Dview(std::cout, y);
      }
      Kokkos::finalize();
    }
