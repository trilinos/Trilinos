KokkosBlas::rotg
################

Defined in header: :code:`KokkosBlas1_rotg.hpp`

.. code:: c++

  template <class execution_space, class SViewType, class MViewType>
  void rotg(
    execution_space const& space,
    SViewType const& a, SViewType const& b,
    MViewType const& c, SViewType const& s);

Compute the Givens rotation coefficients ``c`` and ``s`` that align [a, b] along the direction of :math:`e_1`:

.. math::
   \begin{bmatrix} c && s\\ -s && c\end{bmatrix}\begin{bmatrix}a\\ b\end{bmatrix}=\begin{bmatrix}r\\ 0\end{bmatrix}

satisfying :math:`c^2+s^2=1` and :math:`r=\sqrt{a^2 + b^2}`

Parameters
==========

:space: execution space instance
:a, b: 0-D views. On input, the components of the 2D vector to rotate. On output, ``a`` contains the first component of the rotated vector and ``b`` contains the ``z`` parameter that provides an alternate way to define the rotation.
:c, s: cosine and sine of the rotation that rotates [a, b] onto :math:`e_1`

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `SViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename SViewType::memory_space>::accessible == true``

- `MViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename MViewType::memory_space>::accessible == true``
  - ``!Kokkos::ArithTraits<typename MViewType::value_type>::is_complex``

Example
=======

This example shows how to eliminate an entry using a Givens rotation.
It uses both :code:`rotg` to compute the rotation coefficients and :doc:`rot <blas1_rot>` to apply the rotation.

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

