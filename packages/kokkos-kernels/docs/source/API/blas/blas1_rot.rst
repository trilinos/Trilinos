KokkosBlas::rot
###############

Defined in header: :code:`KokkosBlas1_rot.hpp`

.. code:: c++

  template <class execution_space, class VectorView, class MagnitudeView, class ScalarView>
  void rot(
    execution_space const& space,
    VectorView const& X, VectorView const& Y,
    MagnitudeView const& c, ScalarView const& s);
  
  template <class VectorView, class MagnitudeView, class ScalarView>
  void rot(
    VectorView const& X, VectorView const& Y,
    MagnitudeView const& c, ScalarView const& s);

Applies plane rotation ``(c, s)`` to vector pair ``X, Y``. This can be used to apply a Givens rotation to ``X`` and ``Y``,
if the coefficients ``c`` and ``s`` were previously computed by :doc:`KokkosBlas::rotg(a, b, c, s) <blas1_rotg>`.

1. Replaces the values of ``X`` and ``Y`` with ``c*X + s*Y`` and ``c*Y - s*X`` respectively, using the provided ``space`` instance.
2. Replaces the values of ``X`` and ``Y`` with ``c*X + s*Y`` and ``c*Y - s*X`` respectively, using the default instance of type ``typename VectorView::execution_space``.

The function will throw a runtime exception if ``X.extent(0) != Y.extent(0)``

Parameters
==========

:space: execution space instance
:X, Y: Pair of vectors to rotate
:c, s: cosine and sine of the angle rotation.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_
- `VectorView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename VectorView::memory_space>::accessible == true``
  - ``std::is_same_v<typename VectorView::non_const_value_type, typename VectorView::value_type> == true``

- `MagnitudeView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename MagnitudeView::memory_space>::accessible == true``
  - ``!Kokkos::ArithTraits<typename MagnitudeView::value_type>::is_complex``

- `ScalarView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename ScalarView::memory_space>::accessible == true``

Example
=======

This example shows how to eliminate an entry using a Givens rotation.
It uses :doc:`rotg <blas1_rotg>` to compute the rotation coefficients and :code:`rot` to apply the rotation.

.. code:: c++

    #include <iostream>
    #include <Kokkos_Core.hpp>
    #include <Kokkos_Random.hpp>
    #include "KokkosBlas1_rotg.hpp"
    #include "KokkosBlas1_rot.hpp"
    #include "KokkosKernels_PrintUtils.hpp"
    
    using execution_space = Kokkos::DefaultExecutionSpace;
    using Scalar          = double;
    using Vector          = Kokkos::View<Scalar*, execution_space>;
    using ScalarView      = Kokkos::View<Scalar, execution_space>;
    
    int main(int argc, char* argv[]) {
      Kokkos::initialize();
      {
        const int N = 10;
        Vector x("x", N);
        Vector y("y", N);
    
        // Populate x,y with uniform random values between 0 and 10
        Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
        Kokkos::fill_random(x, rand_pool, Scalar(10));
        Kokkos::fill_random(y, rand_pool, Scalar(10));
    
        std::cout << "x,y before applying Givens rotation:\n";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
        KokkosKernels::Impl::kk_print_1Dview(std::cout, y);
    
        ScalarView c("c");
        ScalarView s("s");
    
        // Calculate Givens rotation coefficients to eliminate y(0)
        KokkosBlas::rotg<execution_space, ScalarView, ScalarView>(execution_space(), Kokkos::subview(x, 0),
                                                                  Kokkos::subview(y, 0), c, s);
    
        std::cout << "\nrotg output (rotation parameters) to eliminate y(0):\n";
        std::cout << "c = ";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, c);
        std::cout << "s = ";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, s);
        std::cout << "r = x(0) = ";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, Kokkos::subview(x, 0));
        std::cout << "z = ";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, Kokkos::subview(y, 0));
    
        // Zero out y(0), which now contains the output parameter z.
        // This completes the replacement of [x(0), y(0)] with [r, 0].
        Kokkos::deep_copy(Kokkos::subview(y, 0), Scalar(0));
    
        // Apply the rotation to the remaining entries of x and y
        KokkosBlas::rot(execution_space(), Kokkos::subview(x, Kokkos::make_pair(1, N)),
                        Kokkos::subview(y, Kokkos::make_pair(1, N)), c, s);
    
        std::cout << "\nx,y after applying Givens rotation:\n";
        KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
        KokkosKernels::Impl::kk_print_1Dview(std::cout, y);
      }
      Kokkos::finalize();
    }

